from common.log.Comment import Comment
from astro.main.PixelMap import PixelMap
from astro.main.RegriddedFrame import CoaddedRegriddedFrame

import argparse
import string
import traceback
import os
import numpy as np
import subprocess
import glob
import fileinput

class KidsCatCreate:
    """Create from a CoaddedRegriddedFrame a KiDS-CAT catalog.

    Based on cvsroot/kidscat/kc_batch.py.

    Requires that the Coadd have a weight frame and a mask.

    Compilation by the private method install_kidscat requires that the
    environment variable AWEPIPE is set, as well as the AWE cvs password.

    Public methods: make and make_batch.

    Private method install_kidscat is used in the case compilation is 
    needed for the computer running the script.

	Instance variables required when running for a single tile are 
    comment_mask, path_kidscat, name_object, and name_filter. The 
    variable comment_mask is a Comment attached to a PixelMap of a mask. 
    The variables name_object and name_filter are an OBJECT and 
    filter.name attributes, respectively.

    The required variable path_kidscat is the full path to where the 
    kidscat shell and executable files are. But also the instance 
    variable install_kidscat can be set to True, so that a local 
    checkout of the KiDS-CAT files is done followed by compilation. In 
    this case, set the variable path_kidscat to where the files should 
    be downloaded.

    Instance variable name_ob is set by method get_name_ob; coadd and 
    mask by methods get_by_comment and get_coadd; weight is that 
    attribute from the coadd; log_errors and log_done by method do.

    The variables log_errors and log_done are text files stored in the 
    kidscat directory. They are used to keep track of tiles processed 
    successfully or unsuccessfully in make_batch. This is done to follow 
    kc_batch.py.

    Ways of calling:

	Call for one tile by passing the tile OBJECT as argument -o, the 
    filter as -f, the comment attached to the mask as -c, and the path 
    where the KiDS-CAT shell files and executables are (or will be) as 
    -p (also doing installation with the flag -i):

	> awe $AWEPIPE/astro/experimental/kids/KidsCatCreate.py -c 'KiDS-INT-DR3_Mask' -p /PATH/TO/DIRECTORY -o KIDS_140.0_-1.5 -f OCAM_u_SDSS -i

    Call same as above but without the installation flag:

    > awe $AWEPIPE/astro/experimental/kids/KidsCatCreate.py -c 'KiDS-INT-DR3_Mask' -p /PATH/TO/DIRECTORY/kidscat -o KIDS_140.0_-1.5 -f OCAM_u_SDSS

	Call for a batch of masks with a given comment as argument -c, by 
    also passing the flag -b:

    > awe $AWEPIPE/astro/experimental/kids/KidsCatCreate.py -b -c 'KiDS-INT-DR3_Mask' -p /PATH/TO/DIRECTORY -i

    Call within AWE for a tile including installation:

    awe> from astro.experimental.kids.KidsCatCreate import KidsCatCreate
    awe> comment_mask = 'KiDS-INT-DR3_Mask'
    awe> path_kidscat = '/PATH/TO/DIRECTORY'
    awe> name_object = 'KIDS_174.0_0.5'
    awe> name_filter = 'OCAM_u_SDSS'
    awe> to_create = KidsCatCreate(comment_mask=comment_mask, path_kidscat=path_kidscat, name_object=name_object, name_filter=name_filter)
    awe> to_create.install_kidscat = True
    awe> newpath, name_ob = to_create.make()

    To continue and ingest the catalog as a SourceList:

    awe> from astro.experimental.kids.KidsCatIngest import KidsCatIngest
    awe> to_ingest = KidsCatIngest(path_kidscat=newpath, name_ob=name_ob, namepart_version_ingest='TEST200')
    awe> sourcelist = to_ingest.ingest_and_attach_comment()

    And then get the statistics:

    awe> from astro.experimental.kids.KidsCat import KidsCat
    awe> to_analyze = KidsCat(sl_slid=sourcelist.SLID)
    awe> statistics = to_analyze.make_statistics_all()
    
    All exceptions handled by public methods by printing the traceback 
    and return none.    
    """

    def __init__(self, comment_mask, path_kidscat, install_kidscat=False,
      name_object=None,
      name_filter=None, name_ob=None, coadd=None, mask=None, weight=None,
      log_errors=None, log_done=None):
        """Class constructor."""
        self.install_kidscat = install_kidscat
        self.path_kidscat = path_kidscat
        self.name_object = name_object
        self.name_filter = name_filter
        self.name_ob = name_ob
        self.comment_mask = comment_mask
        self.coadd = coadd
        self.mask = mask
        self.weight = weight
        self.log_errors = log_errors
        self.log_done = log_done

    def __install_kidscat(self):
        """Checkout kidscat from cvsroot and compile.
        
        Requires that the environment variable AWEPIPE is set, as well 
        as the AWE cvs password.
        
        Changes the value of path_kidscat to the path where the files 
        were downloaded.
        
        Return none.
        """
        os.chdir(self.path_kidscat)

        # Get the checkout of Gert Sikkema's kidscat code.
        if not os.path.isdir(os.path.join(os.getcwd(), 'kidscat')):
            subprocess.call('cvs -d cvs.astro-wise.org:/cvsroot checkout kidscat', 
              shell=True)
              
        else:
            message = '\nPlease pass a path_kidscat that does not contain a directory named kidscat.'
            print message
            
            raise SystemExit

        # Set the kidscat directory.
        dir_kidscat = os.path.join(self.path_kidscat, 'kidscat')

        # Enter the checked out port3 directory.
        dir_port3 = os.path.join(dir_kidscat, 'src/port3')
        os.chdir(dir_port3)

        try:
            subprocess.call('make clean', shell=True)
            subprocess.call('make F77=gfortran', shell=True)

        except:
            subprocess.call('make clean', shell=True)

            raise

        # Enter the checked out src directory.
        dir_src = os.path.join(dir_kidscat, 'src')
        os.chdir(dir_src)

        try:
            subprocess.call('make clean', shell=True)
            subprocess.call('make FC=gfortran CC=gfortran', shell=True)

        except:
            subprocess.call('make clean', shell=True)

            raise

        # Change text in the shell file that will be used, to correspond 
        # to user's variables and to location of dfits executable.
        filename_kidscat = os.path.join(
          dir_kidscat, 'kidscat_v1.6_awe.sh'
        )

        to_replace = {
          'export LD_LIBRARY_PATH':
          'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:%s\n' % dir_src,
          'export PATH':
          'export PATH=%s:${PATH}\n' % dir_src,
          'export SEXHOME':
          'export SEXHOME=%s\n' % dir_kidscat
        }

        with open(filename_kidscat, 'r+') as fp:

            for line in fileinput.input(filename_kidscat):

                if 'dfits' in line:
                    fp.write(
                      """sat=`../dfits ../${FITS}  | grep SATUR | awk """
                    )
                    fp.write(
                      repr('{printf("%e\n",$3)}')
                    )
                    fp.write("`\n")
                    continue

                curr = to_replace.get(line.split('=')[0])

                if curr is not None:
                    fp.write(curr)
                else:
                    fp.write(line)

        fileinput.close()

		# Set instance variable to that expected by other methods.
        self.path_kidscat = dir_kidscat

        message = '\nSuccess downloaded checkout and compiled.'
        print message

    def __rename_files(self):
        """Rename Coadd, compressed and uncompressed mask, and weight.

        Return none.
        """
        new_coadd = '%s.fits' % self.name_ob

        new_mask = '%s.msk.v1.0.all.fits' % self.name_ob

        new_weight = '%s.weight.fits' % self.name_ob

        shutil.move(self.coadd.filename, new_coadd)

        shutil.move(self.mask.filename, new_mask)

        shutil.move('%s.gz' % self.mask.filename, '%s.gz' % new_mask)

        shutil.move(self.weight.filename, new_weight)

        return new_coadd, new_mask, new_weight

    def __get_name_ob(self, coadd):
        """Format KiDS tile name following KiDS-CAT convention.
        
        Return string.
        """
        obj = coadd.OBJECT

        band = string.split(coadd.filter.name, '_')[1]

        return '%s.%s' % (obj, band)

    def __get_effective_gain(self):
        """Get the effective gain.
        
        Use the coadd's RegriddedFrames.
        
        Return float.
        """
        l, t, d2 = [], [], []

        for i in self.coadd.regridded_frames:
            g = i.get_gain()

            d2.append(g.gain)

            l.append(i.FLXSCALE)

            if i.DATE_OBS not in t:
                t.append(i.DATE_OBS)

        effective_gain = np.median(d2) * len(t)/np.median(l)

        return effective_gain

    def __set_filenames(self, rename):
        """Get file names and retrieve, optionally renaming them.
        
        Return 3 strings.
        """
        # Set to always retrieve, thinking that KiDS-CAT cannot run 
        # without the files being stored locally.
        #if retrieve == 1:

        for i in self.mask, self.coadd, self.weight:
            if not os.path.isfile(i.filename):
                i.retrieve()

        # Decompress the mask, keeping the compressed file.
        if not os.path.isfile(self.mask.filename.replace('.gz', '')):
            command='gunzip < %s > %s' %\
              (self.mask.filename, self.mask.filename.replace('.gz', ''))

            os.system(command)

        if rename == 1:

            new_coadd, new_mask, new_weight = self.__rename_files()

            return new_coadd, new_mask, new_weight

        return self.coadd.filename, self.mask.filename.replace('.gz', ''), self.weight.filename

    def __write_done(self, error):
        """Write to log file the name_ob when successfully processed.

        Used by method make_batch.
        
        Return none.
        """
        f = open(self.log_done, 'a').write('%s\n' % self.name_ob)

        return

    def __write_error(self, error):
        """Write to log file the name_ob when unsuccessfully processed.

        Used by method make_batch.
        
        Return none.
        """
        if self.name_ob is not None:
            f = open(self.log_errors, 'a').write('%s\n' % self.name_ob)

        else:
            f = open(self.log_errors, 'a').write(error)

        return

    def __do(self, rename=0, remove=0):
        """Run KiDS-CAT.
        
        First set some variables. Then make call to KiDS-CAT shell 
        script. If remove is 1, remove downloaded files; otherwise move 
        them to the specific tile's directory, named name_ob. Finally 
        prepare the ouput to structure and file names expected by 
        ingestion scripts.
        
        Return none.
        """
        os.chdir(self.path_kidscat)

        self.log_errors = os.path.join(os.getcwd(), 'OBs_error.txt')

        self.log_done = os.path.join(os.getcwd(), 'OBs_done.txt')

        # Set the file names, and optionally rename them by setting 
        # rename to 0 or 1.
        filename_coadd, filename_mask, filename_weight =\
          self.__set_filenames(rename)

        effective_gain = self.__get_effective_gain()

        command = os.path.abspath('./') +\
          '/kidscat_v1.6_awe.sh %s default.sex %s default.param %s %s' %\
          (filename_coadd, effective_gain, filename_mask, filename_weight)

        message = '\nStarting KiDS-CAT\n%s' % command
        print message

        os.system(command)

        if remove == 1:
            for i in [filename_coadd, filename_mask, filename_weight]:
                os.remove(i)

            os.remove(filename_mask.replace('.fits', '.fits.gz'))

        # Move the files to the directory expected by the ingest script.
        else:
            for i in [filename_coadd, filename_mask, 
              filename_mask.replace('.fits', '.fits.gz'), filename_weight]:
                subprocess.call(
                  'mv %s %s/' % (i, filename_coadd.replace('.fits', '')),
                  shell=True
                )

        # Rename the directory holding files to KIDS_<ra>_<dec>.<filter>,
        # as expected by the ingestion script.
        subprocess.call(
          'mv %s %s' % (filename_coadd.replace('.fits', ''), self.name_ob),
          shell=True
        )

        # Rename the .cat and .sex to name expected by ingest script.
        catalog = glob.glob('%s/*.cat' % self.name_ob)[0]

        sextractor = glob.glob('%s/Sci*.sex' % self.name_ob)[0]

        post = catalog.split('_')[-1].replace('.cat', '')

        subprocess.call(
          'mv %s %s/%s_%s.cat' %\
          (catalog, self.name_ob, filename_coadd.replace('.fits', ''), post),
          shell=True
        )

        subprocess.call(
          'mv %s %s/%s_%s.sex' %\
          (sextractor, self.name_ob, filename_coadd.replace('.fits', ''), post),
          shell=True
        )

    def __get_by_comment(self, batch=False):
        """Get PixelMap by given Comment.
        
        If batch is True, return all found. Otherwise, return PixelMap 
        # and Coadd corresponding to name_object and name_filter.
        """
        query_mask = Comment.content == self.comment_mask

        masks = [
          PixelMap(object_id=c.db_object_id) for c in query_mask
        ]

        if batch:
            return masks

        filenames = [i.filename for i in masks]

        filenames = [string.split(i, '.') for i in list(filenames)]

        filenames = ['%s.%s.fits' % (i[0], i[1]) for i in list(filenames)]

        coadd = (CoaddedRegriddedFrame.OBJECT == self.name_object) &\
          (CoaddedRegriddedFrame.filter.name == self.name_filter)

        coadd = [i for i in list(coadd) if i.filename in filenames][0]

        indx = filenames.index(coadd.filename)

        return coadd, masks[indx]

    def __get_coadd(self):
        """Get a CoaddedRegriddedFrame using a PixelMap's filename.
        
        Return CoaddedRegriddedFrame object.
        """
        filename_coadd = '%s.%s.fits' %\
          (self.mask.filename.split('.')[0],
          self.mask.filename.split('.')[1])

        return (CoaddedRegriddedFrame.filename  == filename_coadd)[0]

    def make_batch(self):
        """Run KiDS-CAT for a batch of masks with a given comment.
        
        All exceptions handled by printing the traceback and return none.
        
        Return none.
        """
        try:
            if (self.install_kidscat):
                self.__install_kidscat()

            masks = self.__get_by_comment(batch=True)

            for mask in masks:
                try:
                    self.mask = mask

                    self.coadd = self.__get_coadd()

                    self.weight = self.coadd.weight

                    self.name_ob = self.__get_name_ob(self.coadd)

                    self.__do()

                    self.__write_done()

                except:
                    self.__write_error(traceback.format_exc())

                    continue
            
            return
        
        except:
            traceback.print_exc()

            raise SystemExit

    def make(self):
        """Run KiDS-CAT for a KiDS tile, which has a mask with a given 
        comment.
        
        All exceptions handled by printing the traceback and return none.
        
        Return none.
        """
        try:
            if (self.install_kidscat):
                self.__install_kidscat()

            self.coadd, self.mask = self.__get_by_comment()

            self.weight = self.coadd.weight

            self.name_ob = self.__get_name_ob(self.coadd)

            message = 'Found Coadd with observing_block %s, filename %s.' %\
              (self.coadd.observing_block.name, self.coadd.filename)
            print message

            self.__do()

            return self.path_kidscat, self.name_ob

        except:
            traceback.print_exc()

            raise SystemExit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Provide input OB or input file')
    parser.add_argument('-o', '--ob', help="OB name should be of form KIDS_130.0_0.5")
    parser.add_argument('-f', '--filter', help="Filter name should be of form OCAM_i_SDSS")
    parser.add_argument('-c', '--comment', help="Comment for Mask")
    parser.add_argument('-p', '--path', help="Path where KiDS-CAT shell and executable files are stored, or where they will be stored")
    parser.add_argument('-i', '--install', action='store_const', const=True, default=False, help="Install KiDS-CAT from CVS checkout")
    parser.add_argument('-b', '--batch', action='store_const', const=True, default=False, help="Run script in batch mode")
    args = parser.parse_args()

    # Check that necessary argument was passed.
    if args.comment is None or args.path is None:
        print KidsCatCreate.__doc__

        raise SystemExit

    # If batch is True, run in batch mode.
    if (args.batch):
        to_create = KidsCatCreate(
          comment_mask=args.comment,
          path_kidscat=args.path
        )

        if (args.install):
            to_create.install_kidscat = True

        to_create.make_batch()

    # Elif check that ob and filter were passed.
    elif args.ob is None or args.filter is None:
        print KidsCatCreate.__doc__

        raise SystemExit

    # Else run for the passed comment, ob, and filter.
    else:
        to_create = KidsCatCreate(
          comment_mask=args.comment,
          path_kidscat=args.path,
          name_object=args.ob,
          name_filter=args.filter,
        )

        if (args.install):
            to_create.install_kidscat = True

        to_create.make()

