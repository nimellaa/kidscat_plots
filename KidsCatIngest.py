from astro.util.TableConverter import TableConverter
    
from astro.main.RegriddedFrame import CoaddedRegriddedFrame
from astro.main.SourceList import SourceList
from astro.main.InspectFigure import InspectFigure
from astro.main.Config import SextractorConfig

from astro.external.STIFF import stiff

from common.database.Database import database
from common.database.Context import context
from common.database.typed_list import typed_list

from common.log.Comment import Comment
from common.log.Message import Message

from astro.experimental.kids.qcSourcelist_KiDS_ESO_DR1 import qcSourcelist
from astro.experimental.kids.ProductionControl import *

from pdb import pm
import numpy
import os
import shutil
import argparse
import csv
import traceback
import glob

class KidsCatIngest:
    """Ingest as a SourceList catalogs produced by KiDS-CAT.

    Based on ingest-kidscat-jdj.py.

    Public methods: ingest_and_attach_comment, ingest_inspectfigure.
    Slightly modified methods from ingest-kidscat-jdj.py:
    create_mask_thumbnail_if_needed, create_blending_thumbnail,
    invalidate_and_remove_comment.
    
    Instance variables required are path_kidscat, name_ob, and
    namepart_version_ingest. The string variable path_kidscat is the full 
    path where files (such as the catalog, SExtractor files) are locally 
    stored. The string variable name_ob follows KiDS-CAT format 
    convention, for example KIDS_174.0_0.5.u, and is the directory in the 
    full path where files specific to a KiDS tile are. The string 
    variable namepart_version_ingest will be appended to the ingested
    SourceList's name, for example INTDR3v4.
    
    Instance variables optional is sourcelist_comment, which is a string
    that will be set as a Comment of the SourceList.
    
    Instance variable sourcelist is set by the method 
    ingest_and_attach_comment.

    Ways of calling:

    Call with 1 OB name as argument using -o or --ob, path to files 
    given by argument -p, and string to append to SourceList name by 
    argument -v:

    > awe $AWEPIPE/astro/experimental/kids/KidsCatIngest.py -p /PATH/TO/FILES -o KIDS_174.0_0.5.u -v KCv1.2

    With optional argument to set a Comment for the SourceList, with the
    argument -c:

    > awe $AWEPIPE/astro/experimental/kids/KidsCatIngest.py -p /PATH/TO/FILES -o KIDS_174.0_0.5.u -v KCv1.2 -c 'INT-DR3 SL_Single'

    Call with input file with -f or --file (note: beware of memory 
    issues!):

    > awe $AWEPIPE/astro/experimental/kids/KidsCatIngest.py -f input.list -p /PATH/TO/FILES -v KCv1.2 -c 'INT-DR3 SL_Single'
    
    where input.list contains 1 column with OB names like:
    KIDS_130.0_0.5.i
    
    Make an object within awe:

    awe> from astro.experimental.kids.KidsCatIngest import KidsCatIngest
    awe> to_ingest = KidsCatIngest(path_kidscat=os.getcwd(), name_ob='KIDS_174.0_0.5.u', namepart_version_ingest="KCv1.2")
    awe> sourcelist = to_ingest.ingest_and_attach_comment()

    All exceptions handled by public methods by printing the traceback 
    and return none.
    """

    def __init__(self, path_kidscat, name_ob, namepart_version_ingest,
      sourcelist_comment=None, sourcelist=None
    ):
        """Class constructor."""
        self.path_kidscat = path_kidscat
        self.name_ob = name_ob
        self.namepart_version_ingest = namepart_version_ingest
        self.sourcelist_comment = sourcelist_comment
        self.sourcelist = sourcelist
        
    def __get_phot_apertures_from_file(self, filename):
        """Get photometric apertures included in the catalog.
        
        Argument string filename is the SExtractor configuration file.
        
        Return list.
        """
        # Read in the SextractorConfig and remove comments and extra whitespace.
        lines = [
          l.split("#")[0].strip()
          for l in open(filename).readlines()
        ]

        # Select non-empty lines and split in (key, value) pairs.
        lines = [l.split() for l in list(lines) if len(l) > 0]

        # The join is there for lines like
        # PHOT_AUTOPARAMS	2.5, 3.5
        lines = [(l[0], ' '.join(l[1:])) for l in list(lines)]

        for (k,v) in lines:

            if k == 'PHOT_APERTURES':

                return [(a.strip()) for a in v.split(",")]

    def __create_sexparam_from_file(self, filename):
        """Get the SExtractor parameters.
        
        Argument string filename is the SExtractor file default.param.
        
        Return list.
        """
        lines = [l.strip() for l in open(filename).readlines()]

        lines = [l for l in list(lines) if len(l) > 0]

        lines = [l for l in list(lines) if not l[0] == '#']

        return lines

    def __create_sexconf_from_file(self, filename):
        """Get the SExtractor configuration values.
        
        Argument string filename is the catalog's .sex file.
        
        Return SextractorConfig object.
        """
        # Read in the SextractorConfig and remove comments and extra whitespace.
        lines = [
          l.split("#")[0].strip()
          for l in open(filename).readlines()
        ]

        # Select non-empty lines and split in (key, value) pairs.
        lines = [l.split() for l in list(lines) if len(l) > 0]

        # The join is there for lines like
        # PHOT_AUTOPARAMS	2.5, 3.5
        lines = [(l[0], ' '.join(l[1:])) for l in list(lines)]

        # Convert the list of key, value pairs to SextractorConfig
        # Problematic cases are:
        # PHOT_APERTURES, for which AW allows only 1 value, while KiDSCAT
        #   uses a list.
        # PHOT_FLUXFRAC and WEIGHT_GAIN, which cannot be set in AW, however
        #   the values that KiDSCAT uses are (usually?) the default.
        sexconf = SextractorConfig()

        message = 'Method create_sexconf_from_file: %s'

        for (k,v) in lines:

            if hasattr(sexconf, k):
                class_new = sexconf.__getattribute__(k).__class__

                value_old = sexconf.__getattribute__(k)

                if class_new in [type('string'), type(3.3), type(4)]:
                    try:
                        v2 = class_new(v)
                        sexconf.__setattr__(k, v2)

                    except ValueError, e:
                        print message %\
                          "%s cannot be converted to %s: %s" %\
                          (k, class_new, v)

                        if k == 'PHOT_APERTURES':
                            sexconf.__setattr__(k, -1)
                            print message %\
                              "special case of PHOT_APERTURES, set to -1"

                elif class_new == typed_list:
                    item_types = value_old.item_type

                    try:
                        v2 = [item_types[0](a.strip()) for a in v.split(",")]

                        sexconf.__setattr__(k, v2)

                    except ValueError, e:
                        print message %\
                          "%s cannot be converted to %s of %s: %s" %\
                          (k, class_new, item_types, v)

                else:
                    print message %\
                      "unknown class for %s : %s" % (k, class_new)

            else:
                print message % "%s not in sexconf!" % (k)

        # Should these parameters be set to these defaults?
        #sexconf.WEIGHT_IMAGE = 'weight.fits'

        return sexconf

    def __get_columns(self, tableid):
        """Get the columns in a database table.
        
        Argument integer tableid is the table's number.
        
        Return list.
        """
        tablename = "SOURCELIST*SOURCES**%02i" % (tableid)

        query_table = "SELECT COLUMN_NAME FROM all_tab_columns WHERE TABLE_NAME = '%s' " % (tablename)

        c = database.cursor()

        c.execute(query_table)

        cols_table_result = c.fetchall()

        c.close()

        cols_table = [c[0] for c in cols_table_result]

        cols_table = list(set(cols_table))

        return cols_table

    def __ingest_kidscat_esodr1_sources(self, filename, apertures):
        """Ingest the sources of a KIDSCAT ESO DR1 catalog.
        
        Argument string filename is the catalog. Argument apertures
        is returned by the method __get_phot_apertures_from_file().
        
        Return SourceList object.
        """
        # Read the catalog and split up in header and table.
        data = open(filename).readlines()

        header = [l for l in data if l[0] == '#']

        sources = [l for l in data if l[0] != '#']

        # Parse the headers.
        header = [l.split() for l in list(header)]

        header = [(int(l[1]),l[2],' '.join(l[3:])) for l in list(header)]

        # Add _<number> to mimic vector columns.
        header_new, attributes_vector = [], []

        header_old = header[0]

        for i in header[1:]:

            if i[0] == header_old[0]+1:
                header_new.append(header_old)

            else:
                attributes_vector.append((header_old[1], i[0]-header_old[0]))

                for j in range(i[0] - header_old[0]):
                    header_new.append((
                      header_old[0]+j,
                      '%s_%i' % (header_old[1], j),
                      header_old[2],
                    ))

            header_old = i

        attributes_vector_expected = [('FLUX_APER', 27), ('FLUXERR_APER', 27)]

        assert attributes_vector == attributes_vector_expected, "Vector attributes different than expected %s" % (attributes_vector)

        # Do the same to the last column ('VIGNET')
        vignet = header[-1]

        header_vignet = [
          (vignet[0]+5*x+y, "%s_%s%s" % (vignet[1], x,y), vignet[2])
          for x in range(5) for y in range(5)
        ]

        header = [i for i in header_new] + header_vignet

        # Parse the table.
        sources = [l.split() for l in list(sources)]

        assert set(len(l) for l in sources) == set([len(header)]), "Not correct number of columns."

        sources = [
          [numpy.nan if (a=='nan' or a=='-nan') else eval(a) for a in l]
          for l in list(sources)
        ]

        sources = zip(*sources)

        sources = [numpy.array(l) for l in list(sources)]

        # Create a TableConverter to hold the data.
        tc = TableConverter()

        for h,s in zip(header, sources):
            tc.add_attribute(h[1], s.dtype)

            tc.data[h[1]] = s

            tc.check()

        # There are 27 apertures in the catalog:
        # Keep them, use _<aperturesize> as name
        # The INTDR2 catalogs have FLUX_APER_* instead of MAG_APER_*.
        #magconversion = {a:(a+1)*2 for a in range(25)}
        #magconversion[25] = 100
        #magconversion[26] = 200
        #magconversion
        #{0: 2, 1: 4, 2: 6, 3: 8, 4: 10, 5: 12, 6: 14, 7: 16, 8: 18, 9: 20, 10: 22, 11: 24, 12: 26, 13: 28, 14: 30, 15: 32, 16: 34, 17: 36, 18: 38, 19: 40, 20: 42, 21: 44, 22: 46, 23: 48, 24: 50, 25: 100, 26: 200}
        assert isinstance(apertures[0], str), "Give strings as apertures."

        assert len(apertures) == 27, "Exactly 27 apertures expected."

        magconversion = {
            i: aperture.replace(".", "p")
            for (i, aperture) in enumerate(apertures)
        }

        print magconversion

        # Remove the unnecessary magnitude columns.
        for i in range(27):
            assert i in magconversion, "Missing a magnitude! %s" % (i)

            if not i in magconversion:
                tc.remove_attribute("FLUX_APER_%i" % i)
                tc.remove_attribute("FLUXERR_APER_%i" % i)

        # Rename the necessary magnitude columns.
        # Is this necessary? 
        # Seems to rename one time, then rename again to original.
        for (mag_old,mag_new) in sorted(magconversion.items()):
            tc.rename_attribute("FLUX_APER_%i" %\
              mag_old, "FLUX_APER_R%s" % mag_new
            )
            
            tc.rename_attribute("FLUXERR_APER_%i" %\
              mag_old, "FLUXERR_APER_R%s" % mag_new
            )

        for (mag_old,mag_new) in sorted(magconversion.items()):
            tc.rename_attribute("FLUX_APER_R%s" %\
              mag_new, "FLUX_APER_%s" % mag_new)
            
            tc.rename_attribute("FLUXERR_APER_R%s" %\
              mag_new, "FLUXERR_APER_%s" % mag_new)

        # Remove the 'VIGNET' attributes, since it seems this can better
        # be solved with the cutout server.
        for a in [n for n in tc.attribute_order]:
            if 'VIGNET' in a:
                tc.remove_attribute(a)

        # Rename attributes that have a corresponding column in AW with
        # a different name.
        from astro.main.SourceList import FITSCOLUMN_NAMES

        attrs_rename = FITSCOLUMN_NAMES

        for (attribute_old,attribute_new) in attrs_rename.items():
            
            tc.rename_attribute(attribute_old, attribute_new)

        # Convert to the number format used in AW.
        keep64float = ['RA', 'B', 'Y_WORLD', 'DEC', 'X_WORLD', 'A', 'XM2', 
          'YM2', 'Corr']

        keep64int = ['FLAG']

        for i in tc.attribute_order:

            if (tc.attributes[i]['format'] == 'float64') and (i not in keep64float):

                tc.attributes[i]['format'] = 'float32'

                tc.data[i] = tc.data[i].astype('float32')

            if (tc.attributes[i]['format'] == 'int64') and (i not in keep64int):

                tc.attributes[i]['format'] = 'int32'

                tc.data[i] = tc.data[i].astype('int32')

        # Check whether all necessary columns are available in the 
        # database.
        columns = tc.attribute_order

        #tableid = 0
        #tableid = 9 # Table specific for KiDS
        tableid = 11 # New table specific for KiDS (introduced by Danny in April 2014)

        tablename = "SOURCELIST*SOURCES**%02i" % (tableid)

        columns_exist = [i.lower() for i in self.__get_columns(tableid)]

        columns_new = []

        for i in columns:
            if i.lower() in columns_exist or "sdss_%s" % (i.lower()) in columns_exist:
                pass

            else:
                columns_new.append(i)

        print "==== %-22s %4i %4i %4i" %\
          (tablename, len(columns_exist), len(columns_new), 
          len(columns_exist)+ len(columns_new))

        assert len(columns_new)==0, "New columns required."

        # Save the catalog as a SourceList-compatible FITS file
        sourcelist = tc.save_sourcelist(tablename_store=tablename,
          filename="%s.awcatalog.fits" % filename
        )

        # Set the pathname, wihch is replaced later with a hash
        sourcelist.pathname = "%s.awcatalog.fits" % filename

        return sourcelist

    def create_mask_thumbnail_if_needed(self, pixel_map, coadd, subtype):
        """Create mask InspectFigures.
        
        Return none.
        """
        thumbnail = InspectFigure.select(pixel_map, "PixelMap", subtype, 
          newest=True)

        if thumbnail:
            Message("Found %s for mask: %s" % (subtype,thumbnail.filename), 1)

        else:
            if os.path.exists(pixel_map.filename) == False:
                pixel_map.retrieve()

            if subtype == "Thumbnail__flag":
                thumbnail = InspectFigure(
                  pixel_map, "PixelMap",
                  subtype, extension="png"
                )

                thumbnail.set_attributes()

                thumbnail_plot = ThumbnailPlot(pixel_map, thumbnail.filename)

                thumbnail_plot.process_params.COLOR_MAP = "flag"

                thumbnail_plot.process_params.GAMMA = 1.0

                thumbnail_plot.process_params.CLIP_USING_QUANTILES = True

                thumbnail_plot.process_params.MAX_QUANTILE = 1.0

                thumbnail_plot.process_params.MIN_QUANTILE = 0.0

                thumbnail_plot.make()

            elif subtype == "Thumbnail_stiffb":
                thumbnail = InspectFigure(pixel_map, "PixelMap",
                  subtype, extension="png"
                )

                thumbnail.set_attributes()

                fn_stiff = "%s.tif" % thumbnail.filename

                if not os.path.exists(pixel_map.filename[:-3]):
                    os.system("gunzip %s" % (pixel_map.filename,))

                stiff(pixel_map.filename[:-3],
                  OUTFILE_NAME=fn_stiff, BINNING=12, NEGATIVE='Y',
                  MIN_TYPE='MANUAL', MAX_TYPE='MANUAL', MIN_LEVEL=0,
                  MAX_LEVEL=127, GAMMA=6.0
                )

                os.system("convert %s %s" % (fn_stiff, thumbnail.filename))

            elif subtype == "Thumbnail_coadd_mask_overlay":
                thumbnail = create_blending_thumbnail(coadd, pixel_map)

            if thumbnail != False:
                thumbnail.store()

                thumbnail.commit()

                Message("Created %s for mask: %s" % (subtype, thumbnail.filename), 1)

        return

    def ingest_inspectfigure(self, subtype, sourcelist=None):
        """Ingest plots as InspectFigures.
        
        Argument string subtype is for example, 'ABMAGSAT', 'COMPL',
        'FWHMSN', 'MLIM', 'SG'.
        
        Argument sourcelist is a SourceList for which the InspectFigure
        will be ingested.
        
        Uses the instance variables path_kidscat and name_ob.
        
        Return none.
        
        All exceptions handled by printing the traceback and return none.
        """
        try:
            if self.sourcelist is None and sourcelist is None:
                message = 'Please pass a SourceList object.'
                print message
                
                raise SystemExit
    
            if self.sourcelist is not None:
                sourcelist = self.sourcelist
    
            filename_png = glob.glob("%s/%s/PLOTS/*%s*.png" % (
              self.path_kidscat, self.name_ob, subtype
            ))
            
            assert len(filename_png) == 1, "Cannot find file %s" % ("%s/%s/PLOTS/*%s*.png" % (
              self.path_kidscat, self.name_ob, subtype
            ))
    
            filename_png = filename_png[0].split('/')[-1]
    
            for (filename, extension) in [
              (filename_png, 'png')
            ]:
                ifs = (InspectFigure.filename == filename) &\
                  (InspectFigure.subtype == subtype)
    
                if len(ifs):
                    #figure = None
                    #figure = ifs[0] # This line takes about one minute
                    Message("Found figure %s" % (filename), 1)
    
                else:
                    dir_main = os.getcwd()

                    # Chdir to where the plot is.
                    os.chdir("%s/%s/PLOTS/" % (self.path_kidscat, self.name_ob))
                    
                    figure = InspectFigure(sourcelist, 'SourceList',
                      subtype, extension)
    
                    figure.derive_timestamp()
    
                    figure.set_filename(filename)
    
                    figure.store()
    
                    figure.commit()
    
                    Message("Made inspect figure %s" % (filename), 1)

                    # Chdir to where the script was.
                    os.chdir(dir_main)
                    
            return

        except:
            traceback.print_exc()

            return

    def create_blending_thumbnail(self, coadd, mask):
        """Create blending thumbnail.
        
        Argument coadd is a CoaddedRegriddedFrame and mask is a 
        PixelMap.
        
        Return thumbnail or False.
        """
        thumbnails = coadd.get_inspect_figures()

        thumbnails = [i for i in list(thumbnails) if i.subtype == 'Thumbnail']

        if len(thumbnails):
            thumbnail = thumbnails[0]

            thumbnail.retrieve()

            obname = coadd.observing_block.name

            ob = '%s.%s' % (obname[:-2], obname[-1])

            print '\nCreating blending thumbnail for %s' % ob

            mask_filename = mask.filename[:-3]

            mask.retrieve()

            if not os.path.exists(mask.filename[:-3]):
                os.system("gunzip %s" % (mask.filename,))

            frame = BaseFrame(pathname=mask_filename)

            plot = ThumbnailPlot(frame)

            if os.path.exists(plot.pngname):
                print 'Mask thumbnail %s already exists!' % plot.pngname

            else:
                plot.process_params.MIN_TYPE = 'MANUAL'
                plot.process_params.MIN_LEVEL = 0
                plot.process_params.MAX_TYPE = 'MANUAL'
                plot.process_params.MAX_LEVEL = 1
                plot.process_params.NEGATIVE = True
                plot.process_params.COLOR_MAP = 'Greens'
                plot.make()

            thumbnail_image = thumbnail.filename

            thumbnail_mask = plot.jpgname

            composite_thumbnail = '%s_mask_blend.png' % obname

            os.system('convert %s %s -alpha on -compose dissolve -define compose:args=20 -composite %s' %\
              (thumbnail_image, thumbnail_mask, composite_thumbnail)
            )

            thumbnail = InspectFigure(target=mask, type=mask.__class__.__name__,
              subtype='Thumbnail_coadd_mask_overlay', extension='jpg'
            )

            thumbnail.set_attributes()

            shutil.copyfile(composite_thumbnail,thumbnail.filename)

            print 'Created figure with filename %s' % thumbnail.filename

            return thumbnail

        else:
            return False

    def __get_stats_for_ob(self):
        """Get values from the KiDS-CAT run.
        
        Looks at files in the directories OUTPUT_RES and LOG under the
        ob's path.
        
        Return dictionary.
        """
        stats = {}

        output_res = '%s/%s/OUTPUT_RES' % (self.path_kidscat, self.name_ob)

        filename = os.path.join(output_res, 'MLIM.dat')

        mlim = [[float(a) for a in l.strip().split()]
          for l in open(filename).readlines()[1:]
        ]

        stats['mlim'] = [{
          'SN_LIM': l[0],
          'MLIM': l[1],
          'p_i': l[2],
          }
          for l in list(mlim)
        ]

        filename = os.path.join(output_res, 'seeing.dat')
        
        fwhm = [
          float(a) for a in open(filename).readlines()[0].strip().split()
        ]

        stats['fwhm'] = {
          'fwhm': fwhm[0],
          'sigma': fwhm[1],
        }

        filename = os.path.join(
          '%s/%s/LOG' % (self.path_kidscat, self.name_ob),
          'completeness.log'
        )
        
        compl = float(open(filename).readlines()[0].strip().split()[1])

        stats['completeness'] = {'completeness': compl}

        filename = os.path.join(output_res, 'elliptic.dat')
        
        ell = [float(a) for a in open(filename).readlines()[0].strip().split()]

        stats['ellipticity'] = {
          'ellipticity_selected_stars': ell[0],
          'ellipticity_sure_stars': ell[1],
        }

        filename = os.path.join(output_res, 'abmagsat.dat')
        
        abmagsat = [
          float(a) for a in open(filename).readlines()[0].strip().split()
        ]

        stats['abmagsat'] = {
            'mag_saturation': abmagsat[0],
            'coeff1': abmagsat[1],
            'coeff2': abmagsat[2],
        }

        return stats

    def __set_path(self, file_type):
        """Set the path of a file.
        
        Argument string file_type is the part of the file's name to
        look for.

        Uses the instance variable name_ob.

        Return list.
        """
        for root, dirs, files in os.walk(os.getcwd()):
            if (self.name_ob in root) &\
              (len([i for i in files if file_type in i]) > 0):
                return [os.path.join(root, [i for i in files if file_type in i][0])]

    def ingest_and_attach_comment(self):
        """Ingest as a SourceList and comment a KiDS-CAT catalog.
        
        Ingestion is done only if there is no SourceList with the name
        as built with this script, using the instance variable
        namepart_version_ingest, for the CoaddedRegriddedFrame for which
        KiDS-CAT was run.
        
        The comment is attached only if given by the user and if not 
        already in the database.
        
        Return SourceList object.
        
        All exceptions handled by printing the traceback and return none.
        """
        try:
            context.set_project('KIDS')
            context.set_privileges(1)
        
            Message("Starting SourceList ingest for %s" % (self.name_ob), 1)

            # Determine all the necessary filenames.
            dir_ob = os.path.join(self.path_kidscat, self.name_ob)
            
            filename_catalog = glob.glob(
              os.path.join(dir_ob, '*.cat')
            )

            assert len(filename_catalog) == 1, "Cannot find catalog: %s" %\
              (filename_catalog)

            os.chdir(dir_ob)

            path_catalog = filename_catalog[0]

            filename_catalog = filename_catalog[0].split('/')[-1]

            path_sexconf = glob.glob(
              os.path.join(dir_ob, '%s*.sex' % filename_catalog[:-4])
            )[0]
            
            path_param = glob.glob(
              os.path.join(dir_ob, 'default.param')
            )[0]
            
            filename_coadd = "%s.fits" % (filename_catalog[:-11])

            # Check wether all the files are there.
            for i in [path_catalog, path_sexconf, path_param]:

                assert os.path.exists(i), "Cannot find %s" % (i)

            coadd = (CoaddedRegriddedFrame.filename == filename_coadd)
            
            assert len(coadd) == 1, "Cannot find coadd %s" % (filename_coadd)
            
            coadd = coadd[0]

            Message("All input for mask ingestion checked.", 1)

            # Construct AW mask filename.
            filename_mask_aw = glob.glob(
              os.path.join(dir_ob, '*.msk.v1.0.all*')
            )[0]
            
            # Construct sourcelist name
            name_sourcelist = "KiDS_INTDR3_%s_%s_%s_src_%s_%s" %\
              (self.name_ob[:-2].split("_")[1],
              self.name_ob[:-2].split("_")[2],
              self.name_ob[-1],
              filename_catalog.split("_")[-1][:-4],
              self.namepart_version_ingest)

            sls = (SourceList.name == name_sourcelist) &\
              (SourceList.frame == coadd)

            if len(sls):
                sl = sls[0]
                Message("Found the SourceList in database: %s" % (sl.filename), 1)

            else:
                Message("SourceList not found; ingesting SourceList %s" %\
                  (name_sourcelist), 1)

                # The SextractorConfiguration
                sexconf = self.__create_sexconf_from_file(path_sexconf)

                phot_apertures = self.__get_phot_apertures_from_file(
                  path_sexconf)

                phot_apertures_expected = [
                  '2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '22', 
                  '24', '26', '28.5', '30', '32', '34', '36', '38', '40', '42', 
                  '44', '46', '48', '50', '100', '200'
                ]

                assert phot_apertures == phot_apertures_expected, "phot_apertures unexpected %s" % (phot_apertures)

                if filename_mask_aw[-3:] == ".gz":
                    filename_sex_mask = filename_mask_aw[:-3]

                else:
                    filename_sex_mask = filename_mask_aw

                filename_sex_mask = filename_sex_mask.split('/')[-1]

                sexconf.FLAG_IMAGE = filename_sex_mask

                sexconf.WEIGHT_IMAGE = coadd.weight.filename

                sexconf.CATALOG_TYPE = 'FITS_LDAC'

                # The SextractorParameters
                sexparam = self.__create_sexparam_from_file(path_param)

                sexparam = [p for p in list(sexparam) if not 'VIGNET' in p]

                sourcelist = self.__ingest_kidscat_esodr1_sources(
                  filename_catalog, phot_apertures
                )

                sourcelist.frame = coadd

                sourcelist.copy_attributes()

                sourcelist.sexconf = sexconf

                sourcelist.sexparam = sexparam

                sourcelist.name = name_sourcelist

                # Store the corresponding FITS file for convenience
                sourcelist.store_with_hash_as_name(
                  prefix="%s-" % (sourcelist.name), suffix='fits'
                )

                sourcelist.commit()

                Message("SourceList created: %s" % (sourcelist.filename), 1)

            if self.sourcelist_comment is not None:
                production_control = ProductionControl()

                production_control.addCommentIfNeeded(
                  sourcelist, self.sourcelist_comment)

            # Source distribution inspect figure
            try:
                qcS = qcSourcelist(
                    sourcelist=sourcelist,
                    redo=False,
                    privileges=sourcelist._privileges,
                    commit=1
                )
                qcS.slSourceDistribution()

            except:
                traceback.print_exc()
                
                pass
                
            # SourceList statistics. Surround with try/except for the 
            # case there were files not produced by the KiDS-CAT run.
            try:
                stats = self.__get_stats_for_ob()
                
                stats1 = {k: stats[k] for k in ['mlim', 'completeness', 'fwhm']}
                
                stats2 = {k: stats[k] for k in ['ellipticity','abmagsat']}
                
                for i in [stats1, stats2]:
                    str_comment = repr(i)
                
                    assert len(str_comment) < 297, "Comment to long: %i %s" %\
                      (len(str_comment), str_comment)
                
                    production_control.addCommentIfNeeded(
                      sourcelist, str_comment)

            except:
                traceback.print_exc()
                
                pass

            Message("Finished SourceList ingest for %s SLID %s" %\
              (self.name_ob, sourcelist.SLID), 1)
            
            self.sourcelist = sourcelist
            
            return sourcelist

        except:
            traceback.print_exc()

            return

    def invalidate_and_remove_comment(self):
        """Invalidate a SourceList and remove given comment.
        
        Return none.
        
        All exceptions handled by printing the traceback and return none.
        """
        try:
            Message("Starting mask and sourcelist invalidation for %s" %\
              (self.name_ob),1)

            filename_catalog = glob.glob(
              os.path.join(dir_ob, '*.cat')
            )
            
            assert len(filename_catalog) == 1, "Cannot find catalog: %s" %\
              (filename_catalog)

            path_catalog = filename_catalog[0]

            filename_catalog = filename_catalog[0].split('/')[-1]

            path_sexconf = glob.glob(
              os.path.join(dir_ob, '%s*.sex' % filename_catalog[:-4])
            )[0]
            
            path_param = glob.glob(
              os.path.join(dir_ob, 'default.param')
            )[0]
            
            filename_coadd = "%s.fits" % (filename_catalog[:-11])

            # Check wether all the files are there.
            for fn in [path_catalog]:
                assert os.path.exists(fn), "Cannot find %s" % (fn)

            coadd = (CoaddedRegriddedFrame.filename == filename_coadd)
            
            assert len(coadd) == 1, "Cannot find coadd %s" % (filename_coadd)
            
            coadd = coadd[0]

            Message("All input checked.",1)

            # Construct sourcelist name
            name_sourcelist = "KiDS_INTDR3_%s_%s_%s_src_%s_%s" %\
              (self.name_ob[:-2].split("_")[1],
              self.name_ob[:-2].split("_")[2],
              self.name_ob[-1],
              filename_catalog.split("_")[-1][:-4],
              self.namepart_version_ingest)

            # Test if source list exists in database
            sourcelist = (SourceList.name == name_sourcelist) &\
              (SourceList.frame == coadd)

            if len(sourcelist):

                sourcelist = sourcelist[0]

                context.update_is_valid(sourcelist,0)

                production_control = ProductionControl()

                production_control.deleteCommentIfNeeded(
                  sourcelist,self.sourcelist_comment)

                Message("SourceList in database is now invalidated: %s" %\
                  (sourcelist.filename), 1)

            else:

                Message("SourceList not found: %s" % (name_sourcelist), 1)

            return

        except:
            traceback.print_exc()

            return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Provide input OB or input file')
    parser.add_argument('-f', '--file', help="Input file should contain 1 column with OB names (KIDS_130.0_0.5.i)")
    parser.add_argument('-o', '--ob', help="OB name should be of form KIDS_130.0_0.5.i")
    parser.add_argument('-p', '--path', help="Full path where files to ingest are stored.")
    parser.add_argument('-v', '--version', help="Version number which will be kept in SourceList name.")
    parser.add_argument('-c', '--comment', help="Comment for SourceList.")
    parser.add_argument('-i', '--invalidate', action='store_const', const=True, default=False, help="Invalidate masks and sourcelists")
    args = parser.parse_args()

    context.set_project('KIDS')
    context.set_privileges(1)

    # Check that required arguments were passed.
    if args.path is None:
        print KidsCatIngest.__doc__        

        raise SystemExit
    
    if args.version is None:
        print KidsCatIngest.__doc__        

        raise SystemExit

    if (args.invalidate):
        invalidate = True

    else:
        invalidate = False
        
    # If file is provided.
    if (args.file):
        # open input file
        inputfile = open(args.file, 'rb')

        inputreader = csv.reader(inputfile, delimiter=',', quotechar='"')

        # read input file
        obs = []

        for row in inputreader:
            obs.append(row[0])

        inputfile.close()

        # call ingest or invalidate method for each OB
        for ob in obs:
            to_ingest = KidsCatIngest(
              path_kidscat=args.path,
              name_ob=ob,
              namepart_version_ingest=args.version,
              sourcelist_comment=args.comment
            )

            if (invalidate):
                try:

                    to_ingest.invalidate_and_remove_comment()

                except:
                    traceback.print_exc()

                    print "ERROR FOR %s" % (ob.name)
            else:
                try:
                    to_ingest.ingest_and_attach_comment()

                except:
                    traceback.print_exc()

                    print "ERROR FOR %s" % (ob.name)

    # If OB name is provided.
    elif (args.ob):

        to_ingest = KidsCatIngest(
          path_kidscat=args.path,
          name_ob=args.ob,
          namepart_version_ingest=args.version,
          sourcelist_comment=args.comment
        )
        
        if (invalidate):
            to_ingest.invalidate_and_remove_comment()

        else:
            to_ingest.ingest_and_attach_comment()

    # If nothing is provided.
    else:

        print KidsCatIngest.__doc__

