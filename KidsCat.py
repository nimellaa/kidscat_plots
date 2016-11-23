"""Contains classes related to implementation of KiDS-CAT
(La Barbera et al).

Use class KidsCat to calculate individual statistics, calculate all
statistics, make all plots, or calculate all statistics while also
making all plots.

Use class KidsCatBase as base for child classes KidsCatAbmagsat,
KidsCatCompl, KidsCatFwhmsn, KidsCatMlim, KidsCatSg.
"""

from astro.main.SourceList import SourceList
from astro.util.TableConverter import TableConverter

import os
import glob
import subprocess
from datetime import datetime
import traceback
import sys
import numpy as np

# Used here to ignore numpy warning:
# 'Warning: invalid value encountered in log10'
np.seterr(all='ignore')

class KidsCat:
    """Useful methods related to KiDS-CAT for calculating statistics and
    making plots.

    Instance variables string sl_object, string sl_filter, string
    sl_filename, and integer sl_slid correspond to attributes of
    SourceList OBJECT, filter.name, filename or filename.like, and SLID,
    respectively. Construct with either SLID, or OBJECT and filter.name,
    or OBJECT, filter.name, and filename. They will be used to query AWE
    for a SourceList. If a SLID is given, it will be used. Otherwise the
    other variables will be used and giving all three will result in the
    most specific query.

    Examples:

    > sl_object = 'KIDS_140.0_-1.5'
    > sl_filter = 'OCAM_r_SDSS'

    Useful attribute to specify further the SourceList that will be
    used:

    > sl_filename = '*KCv1.6_INTDR3v4*'

    Or pass a SourceList.SLID:

    > sl_slid = 11395661

    Public methods are get_ellip_median_all, get_ellip_median_stars,
    get_magsat, get_mcompl, get_fwhm_mean, get_fwhm_sigma, get_mlim_sn5,
    get_mlim_sn10, get_mlim_sn15, make_statistics_all, make_plot_all,
    and make_all.

    There are two main ways to get the statistics calculated by KiDS-CAT:
    get a specific statistic or get all the statistics.

    To get a specific statistic use any of the methods:

    > get_ellip_median_all
    > get_ellip_median_stars
    > get_magsat
    > get_mcompl
    > get_fwhm_mean
    > get_fwhm_sigma
    > get_mlim_sn5
    > get_mlim_sn10
    > get_mlim_sn15

    which take no arguments. Return a float.

    The get all of the statistics at once use the method:

    > make_statistics_all

    which takes no arguments. Return dictionary with string keys
    ellip_median_all, ellip_median_stars, magsat, fwhm_mean, fwhm_sigma,
    fwhm_mean_stars, fwhm_sigma_stars, mcompl, mlim_sn5, mlim_sn10, and
    mlim_sn15.

    To make all of the KiDS-CAT plots at once use the method:

    > make_plot_all

    which takes no arguments. Stores the plots in current working
    directory.

    To get all of the statistics and make all of the plots use the
    method:

    > make_all

    which takes no arguments. Return same dictionary as the method
    make_statistics_all.

    Each of the public methods handle all exceptions by printing the
    traceback and return none.
    """

    def __init__(self, sl_object=None, sl_filter=None, sl_filename=None,
      sl_slid=None
    ):
        """Class constructor."""
        self.sl_object = sl_object
        self.sl_filter = sl_filter
        self.sl_filename = sl_filename
        self.sl_slid = sl_slid

    def get_ellip_median_all(self):
        """Get the median of the ellipticity of all sources.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_abmagsat(['ellip_median_all'])

            return dic['ellip_median_all']

        except:
            traceback.print_exc()

            return

    def get_ellip_median_stars(self):
        """Get the median of the ellipticity of sure stars.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_abmagsat(['ellip_median_stars'])

            return dic['ellip_median_stars']

        except:
            traceback.print_exc()

            return

    def get_magsat(self):
        """Get the saturation magnitude.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_abmagsat(['magsat'])

            return dic['magsat']

        except:
            traceback.print_exc()

            return

    def get_mcompl(self):
        """Get the completeness magnitude.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_compl()

            return dic['mcompl']

        except:
            traceback.print_exc()

            return

    def get_fwhm_mean(self):
        """Get the biweight location of the FWHM, in pixels.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_abmagsat(['fwhm_mean'])

            message = '\nFinished get_fwhm_mean. Value returned is in pixels.'
            print message

            return dic['fwhm_mean']

        except:
            traceback.print_exc()

            return

    def get_fwhm_sigma(self):
        """Get the biweight scale of the FWHM, in pixels.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_abmagsat(['fwhm_sigma'])

            message = '\nFinished get_fwhm_sigma. Value returned is in pixels.'
            print message

            return dic['fwhm_sigma']

        except:
            traceback.print_exc()

            return

    def get_mlim_sn5(self):
        """Get the limiting magnitude at a S/N of 5.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_mlim(['mlim_sn5'])

            return dic['mlim_sn5']

        except:
            traceback.print_exc()

            return

    def get_mlim_sn10(self):
        """Get the limiting magnitude at a S/N of 10.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_mlim(['mlim_sn10'])

            return dic['mlim_sn10']

        except:
            traceback.print_exc()

            return

    def get_mlim_sn15(self):
        """Get the limiting magnitude at a S/N of 15.

        Returns a float.
        """
        try:
            dic = self.__get_statistics_mlim(['mlim_sn15'])

            return dic['mlim_sn15']

        except:
            traceback.print_exc()

            return

    def __get_statistics_abmagsat(self, in_statistics):
        """Called by get_ellip_median_all, get_ellip_median_stars,
        get_magsat.

        Get one or all of the statistics calculated by KidsCatAbmagsat,
        specified by the argument in_statistics, which is a list.

        Returns a dictionary with the statistic(s) as string key(s).
        """
        from astro.experimental.kids.KidsCatAbmagsat import KidsCatAbmagsat

        plot = KidsCatAbmagsat(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )
        plot.make_statistics()

        statistics = {
          'ellip_median_all': plot.ellip_median_all,
          'ellip_median_stars': plot.ellip_median_stars,
          'magsat': plot.magsat,
          'fwhm_mean': plot.fwhm_mean,
          'fwhm_sigma': plot.fwhm_sigma,
          'fwhm_mean_stars': plot.fwhm_mean_stars,
          'fwhm_sigma_stars': plot.fwhm_sigma_stars,
        }

        return_statistics = {}

        for i in in_statistics:
            return_statistics[i] = statistics[i]

        return return_statistics

    def __get_statistics_compl(self):
        """Called by get_mcompl.

        Get the statistic calculated by KidsCatCompl.

        Returns a dictionary with the statistic as string key.
        """
        from astro.experimental.kids.KidsCatCompl import KidsCatCompl

        plot = KidsCatCompl(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )
        plot.make_statistics()

        statistics = {}

        statistics['mcompl'] = plot.mcompl

        return statistics

    def __get_statistics_mlim(self, in_statistics):
        """Called by get_mlim_sn5, get_mlim_sn10, get_mlim_sn15.

        Get one or all of the statistics calculated by KidsCatMlim,
        specified by the argument in_statistics, which is a list.

        Returns a dictionary with the statistic(s) as string key(s).
        """
        from astro.experimental.kids.KidsCatMlim import KidsCatMlim

        plot = KidsCatMlim(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )
        plot.make_statistics()

        statistics = {
          'mlim_sn5': plot.mlim_sn5,
          'mlim_sn10': plot.mlim_sn10,
          'mlim_sn15': plot.mlim_sn15
        }

        return_statistics = {}

        for i in in_statistics:
            return_statistics[i] = statistics[i]

        return return_statistics

    def make_statistics_all(self):
        """Get all statistics produced by KiDS-CAT.

        Returns a dictionary, where the string keys are:
        ellip_median_all, ellip_median_stars, magsat, fwhm_mean,
        fwhm_sigma, fwhm_mean_stars, fwhm_sigma_stars, mcompl, mlim_sn5,
        mlim_sn10, mlim_sn15.
        """
        try:
            statistics = self.__do_all(level='statistics')

            return statistics

        except:
            traceback.print_exc()

            return

    def make_plot_all(self):
        """Use to make all the KiDS-CAT plots, without calculating
        unnecessary statistics.

        Each of the plots is saved locally.
        """
        try:
            self.__do_all(level='plot')

            return

        except:
            traceback.print_exc()

            return

    def make_all(self):
        """Make all of the KiDS-CAT plots, including the statistics
        corresponding to each of the plots.

        Each of the plots is saved locally.

        Returns a dictionary, where the string keys are:
        ellip_median_all, ellip_median_stars, magsat, fwhm_mean,
        fwhm_sigma, fwhm_mean_stars, fwhm_sigma_stars, mcompl, mlim_sn5,
        mlim_sn10, mlim_sn15.
        """
        try:
            statistics = self.__do_all(level='make')

            return statistics

        except:
            traceback.print_exc()

            return

    def __do_all(self, level):
        """Method called by make_plot_all, make_statistics_all, and
        make_all.

        The method make_plot_all passes the argument 'plot'. The method
        make_statistics_all passes the argument 'statistics'. The method
        make_all passes the argument 'make'.

        Use correspondingly to make all plots, or all statistics of
        KiDS-CAT, or all plots and corresponding statistics. It actually
        runs the methods of the child classes, creating an object of
        each of the classes.

        Each of the plots is saved locally.

        If the argument is 'make' or 'statistics', returns a dictionary
        where the string keys are: ellip_median_all, ellip_median_stars,
        magsat, fwhm_mean, fwhm_sigma, fwhm_mean_stars, fwhm_sigma_stars,
        mcompl, mlim_sn5, mlim_sn10, and mlim_sn15.
        """
        from astro.experimental.kids.KidsCatAbmagsat import KidsCatAbmagsat
        from astro.experimental.kids.KidsCatCompl import KidsCatCompl
        from astro.experimental.kids.KidsCatFwhmsn import KidsCatFwhmsn
        from astro.experimental.kids.KidsCatMlim import KidsCatMlim
        from astro.experimental.kids.KidsCatSg import KidsCatSg

        if level == 'make':
            get_stats = True
            write_to_file = True
            do_plot = True

        elif level == 'statistics':
            get_stats = True
            write_to_file = False
            do_plot = False

        elif level == 'plot':
            get_stats = False
            write_to_file = False
            do_plot = True

        plot = KidsCatAbmagsat(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )

        statistics = {}

        stats = plot.make(get_stats, write_to_file, do_plot)

        # Set the path so that it is not asked for every KidsCat subclass.
        path = plot.path

        # Set the attribute sourcelist.
        self.sourcelist = plot.sourcelist

        if get_stats:
            statistics.update(stats)

        plot = KidsCatCompl(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )

        plot.path = path

        stats = plot.make(get_stats, write_to_file, do_plot)

        if get_stats:
            statistics.update(stats)

        plot = KidsCatFwhmsn(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )

        plot.path = path

        stats = plot.make(get_stats, write_to_file, do_plot)

        if get_stats:
            statistics.update(stats)

        plot = KidsCatMlim(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )

        plot.path = path

        stats = plot.make(get_stats, write_to_file, do_plot)

        if get_stats:
            statistics.update(stats)

        plot = KidsCatSg(
          sl_object=self.sl_object, sl_filter=self.sl_filter,
          sl_filename=self.sl_filename, sl_slid=self.sl_slid
        )

        plot.path = path

        stats = plot.make(get_stats, write_to_file, do_plot)

        if get_stats:
            statistics.update(stats)

        if get_stats:
            return statistics

        return

class KidsCatBase:
    """Base for classes related to KiDS-CAT.

    Public methods: make, make_statistics, and make_plot. Public methods
    meant to be called only by subclasses: find_sourcelist, load_data,
    print_message, bin_mat, find_stars, get_attribute.

    Instance variables string sl_object, string sl_filter, string
    sl_filename, and integer sl_slid correspond to attributes of
    SourceList OBJECT, filter.name, filename or filename.like, and SLID,
    respectively. Construct with either SLID, or OBJECT and filter.name,
    or OBJECT, filter.name, and filename. They will be used to query AWE
    for a SourceList. If a SLID is given, it will be used. Otherwise the
    other variables will be used and giving all three will result in the
    most specific query.

    Examples:

    > sl_object = 'KIDS_140.0_-1.5'
    > sl_filter = 'OCAM_r_SDSS'

    Useful attribute to specify further the SourceList that will be
    used:

    > sl_filename = '*KCv1.6_INTDR3v4*'

    Or pass a SourceList.SLID:

    > sl_slid = 11395661

    The attribute data contains the columns that are necessary for the
    plot. It is set by in each of the plots calling the method
    load_data, so it is run before making each of the plots.

    The attributes fwhm_mean, fwhm_sigma, fwhm_mean_stars, and
    fwhm_sigma_stars are returned by the method find_stars. It is
    called by KidsCatFwhmsn and called by KidsCatAbmagsat, which also
    uses a file created by find_stars.

    The attributes ellip_median_all and ellip_median_stars are set after
    the make methods of KidsCatABMAGSAt.

    The attribute mcompl is set by the make methods of KidsCatCompl.

    The attributes mlim_sn5, mlim_sn10, and mlim_sn15 are set by the
    make methods of KidsCatMlim.

    The attributes magsat and msat are set by the make methods of
    KidsCatAbmagsat.

    The attribute path is set by the method set_path. It is the main path 
    where files will be saved. It contains the SourceList OBJECT and the 
    filter following the KiDS-CAT format, e.g. KIDS_140.0_-1.5.u. It is 
    joined with OUTPUT_RES or TEMP for some files, and with PLOTS for
    the plots. If in the current path there is no directory named in
    this format, the main path will be set to a structure following the
    AWE file structure. That is, path including environment variable
    $AWEPIPE (or otherwise to path including the current working
    directory) joined with astro/experimental/kids/kidscat_temp/ or with 
    astro/experimental/kids/kidscat_output_res/.

    Make methods handle all exceptions by printing the traceback and
    return none.
    """

    def __init__(self, sl_object=None, sl_filter=None, sl_filename=None,
      sl_slid=None, sourcelist=None, data=None, fwhm_mean=None, fwhm_sigma=None,
      fwhm_mean_stars=None, fwhm_sigma_stars=None, ellip_median_all=None,
      ellip_median_stars=None, mcompl=None, mlim_sn5=None, mlim_sn10=None,
      mlim_sn15=None, msat=None, magsat=None, path=None
    ):
        """Class constructor."""
        self.sl_object = sl_object
        self.sl_filter = sl_filter
        self.sl_filename = sl_filename
        self.sl_slid = sl_slid
        self.sourcelist = sourcelist
        self.data = data
        self.fwhm_mean = fwhm_mean
        self.fwhm_sigma = fwhm_sigma
        self.fwhm_mean_stars = fwhm_mean_stars
        self.fwhm_sigma_stars = fwhm_sigma_stars
        self.ellip_median_all = ellip_median_all
        self.ellip_median_stars = ellip_median_stars
        self.mcompl = mcompl
        self.mlim_sn5 = mlim_sn5
        self.mlim_sn10 = mlim_sn10
        self.mlim_sn15 = mlim_sn15
        self.msat = msat
        self.magsat = magsat
        self.path = path

    def __mad(self, list_in):
        """Calculates the Median Absolute Deviation of a list of numbers.

        The argument list_in is the list of numbers.

        Return float median absolute deviation.

        Should be called by method __biweight_location.
        """
        list_median = np.median(list_in)

        # Calculate |x-M| values
        diff_moduli = []

        for i in list_in:
            diff_moduli.append(abs(i - list_median))

        mad = np.median(diff_moduli)

        return mad

    def __biweight_location(self, list_in, num_sigma):
        """Calculates the biweight location estimator (like a robust
        average) of a list of numbers.

        Implemented from the function biweightLocation from the module
        astStats, see:
        http://astlib.sourceforge.net/docs/astLib/astLib.astStats-module.html

        The argument list_in is the list of numbers. The argument
        num_sigma is a float to indicate the number of standard 
        deviations from the MAD that we want to keep.

        Return float biweight location.

        Should be called by method find_stars.
        """
        tuning_constant = num_sigma / 0.6745
        
        list_median = np.median(list_in)

        list_mad = self.__mad(list_in)

        if list_mad != 0:
            u_values = []

            for i in list_in:
                u_values.append(
                  (i - list_median) / (tuning_constant * list_mad)
                )

            top = 0     # numerator equation (5) Beers 1990AJ....100...32B
            bottom = 0  # denominator equation (5) Beers 1990AJ....100...32B

            for i in range(len(u_values)):

                if abs(u_values[i]) <= 1.0:
                    
                    top += (
                      (list_in[i] - list_median) \
                      * (1.0 - (u_values[i]**2))**2 \
                    )

                    bottom += (
                      (1.0 - (u_values[i]**2))**2
                    )

            cbi = list_median + (top / bottom)

        return cbi

    def __biweight_scale(self, list_in, num_sigma):
        """Calculates the biweight scale estimator (like a robust
        standard deviation) of a list of numbers.

        Implemented from the function biweightScale from the module
        astStats, see:
        http://astlib.sourceforge.net/docs/astLib/astLib.astStats-module.html

        The argument list_in is the list of numbers.The argument
        num_sigma is a float to indicate the number of standard 
        deviations from the MAD that we want to keep.

        Return float biweight scale.

        Should be called by method find_stars.
        """
        tuning_constant = num_sigma / 0.6745
        
        # Calculate |x-M| values and u values
        list_median = np.median(list_in)

        list_mad = self.__mad(list_in)

        diff_moduli = []

        for i in list_in:
            diff_moduli.append(abs(i - list_median))

        u_values = []

        for i in list_in:
            u_values.append(
              (i - list_median) / (tuning_constant * list_mad)
            )

        top = 0         # numerator equation (9) Beers 1990AJ....100...32B
        bottom = 0      # denomination equation (9) Beers 1990AJ....100...32B
        ctr = 0         # Count values where u<1 only

        for i in range(len(u_values)):

            # Skip u values >1
            if abs(u_values[i]) <= 1.0:
                
                u2Term = 1.0 - (u_values[i]**2)

                u4Term = u2Term**4

                top += (
                  (diff_moduli[i]**2) * u4Term
                )

                bottom += (
                  u2Term * (1.0 - (5.0 * (u_values[i]**2)))
                )

                ctr += 1

        top = np.sqrt(top)
        
        bottom = abs(bottom)

        sbi = float(ctr)**0.5 * (top / bottom)

        return sbi

    def __chi(self, diff, num_sigma):
        """Calculates part of the biweight scale.

        Implemented from chi subroutine by Naples KiDS-CAT
        (La Barbera et al.).

        The argument diff is the value to compare. The argument
        num_sigma is a float to indicate the number of standard 
        deviations from the MAD that we want to keep.
        
        Return float.

        Should be called by method bwsm_loop.
        """
        if abs(diff) <= num_sigma:
            chi = diff**2 / 2.
        
        else:
            chi = num_sigma**2 / 2.
        
        return chi

    def __psi(self, diff):
        """Calculates part of the biweight location.

        Implemented from psi subroutine by Naples KiDS-CAT
        (La Barbera et al.).

        The argument diff is the value to compare.

        Return float.

        Should be called by method bwsm_loop.
        """
        if abs(diff) <= 1.:
            psi = diff * (1. - diff**2)**2
        
        else:
            psi = 0.
        
        return psi

    def __bwsm_loop(self, n, x, num_sigma, par_2, t1, s1):
        """Loop used by method bwsm.

        Implemented from bwsm subroutine by Naples KiDS-CAT
        (La Barbera et al.).

        Return list of floats with biweight location and scale.

        Should be called by method bwsm.
        """
        s2 = 0.
        
        for i in x:
            s2 += self.__chi(
              (i - t1) / s1, num_sigma
            )
        
        s2 = np.sqrt(s2) * s1 / np.sqrt(par_2 * float(n-1))
        
        t2 = 0.
        
        for i in x:
            t2 += self.__psi((i - t1) / s2)
        
        t2 = t1 + t2 * s2/float(n)
        
        return [t2, s2]

    def __bwsm(self, xx, num_sigma):
        """Calculates the biweight location and scale estimator.

        Implemented from bwsm subroutine by Naples KiDS-CAT
        (La Barbera et al.).

        The argument xx is the list of numbers. The argument
        num_sigma is a float to indicate the number of standard 
        deviations from the location that we want to keep.

        Return list of floats with biweight location and scale.

        Should be called by method bin_mat.
        """
        n = len(xx)
        
        x = []
      
        # Some safety precautions in case there are very few sources.
        if n == 0:
            theta = -1.
            sigma = -1.
    
        if n == 1:
            theta = x[0]
            sigma = 0.
    
        kc = 0
    
        for i in xx:

            if i == xx[0]:
                kc += 1
    
        if kc == n:
            theta = xx[0]
            sigma = 0.1 * abs(theta)
    
        if (n <= 1) or (kc == n):
            return [theta, sigma]
    
        for i in xx:
            x.append(i)
    
        # Comes from a double integral of a function. Value was 
        # calculated using scipy.
        par_2 = 0.488780198154

        # Calculate iteratively the biweight location and scale.
        t1 = np.median(x)
    
        s1 = np.std(x)
            
        k = 0
        
        t2, s2 = self.__bwsm_loop(n, x, num_sigma, par_2, t1, s1)

        while (abs((t1-t2)/t1) > 0.0001) or (abs((s1-s2)/s2) > 0.0001):
            t1 = t2
            s1 = s2
            
            if k >= 200:
                return [t1, s1]
                
            else:
                t2, s2 = self.__bwsm_loop(n, x, num_sigma, par_2, t1, s1)
            
            k += 1
            
        theta = t1
          
        sigma = s1
        
        return [theta, sigma]

    def set_path(self):
        """Set path where files will be saved.

        Uses the instance variables sl_object and sl_filter to find a
        directory that contains them in the current path.

        Return string.
        """
        name_ob = '%s.%s' % (
          self.sl_object, 
          self.sl_filter.replace('OCAM', '').replace('SDSS', '').replace('_', '')
        )
        
        candidates = []
        
        for root, dirs, files in os.walk(os.getcwd()):
            if (name_ob in root.split('/')[-1]):#) & (dir_type in root):
                #return root
                candidates.append(root)

        # No directory found. Use path that follows AWE file structure.
        if len(candidates) == 0:
            path = os.path.join(
              os.getenv('AWEPIPE', os.getcwd()),
              'astro/experimental/kids/%s/' %\
              name_ob
            )
            
            message = 'No %s directory found. Saving files to %s.' %\
              (name_ob, path)
              
            return path
            
        # Exactly one directory found, so return it.
        elif len(candidates) == 1:
            return candidates[0]
        
        # Multiple directories found. Let the user choose from those
        # found.
        else:
            
            path_index = -9999
            
            while path_index-1 not in range(len(candidates)):
                print '\nChoose a path:'
                
                for i in range(len(candidates)):
                    print '%s. %s' % (i+1, candidates[i])
                    
                path_index = raw_input('Enter a number from the list: ')
                
                try:
                    path_index = int(path_index)
                    
                except:
                    path_index = -9999
            
            return candidates[path_index-1]
            
    def find_sourcelist(self):
        """Query for SourceList from which the plots and statistics will
        be made.

        Return SourceList object.

        Should be called by child classes.
        """
        message = '\nObject has\nsl_object: %s\nsl_filter: %s\nsl_filename: %s\nsl_slid: %s\nConstruct with at least strings SourceList OBJECT and filter.name, or an integer SourceList SLID.\n'

        # Check if there is a SLID, or a OBJECT and filter.name.
        if (self.sl_object is None or self.sl_filter is None) and\
          (self.sl_slid is None):
            print message %\
              (self.sl_object, self.sl_filter, self.sl_filename,
              self.sl_slid)

            raise SystemExit

        # If a SourceList SLID is given, use it.
        # Query for a SourceList and if not found, print message and
        # exit; otherwise return the found SourceList.
        if self.sl_slid is not None:
            sourcelist = SourceList.SLID == self.sl_slid

            if not len(sourcelist):
                message = '\nNo SourceList found with SLID %s.\n' %\
                  self.sl_slid
                print message

                raise SystemExit

            else:
                sourcelist = sourcelist[0]

        elif self.sl_object is not None and self.sl_filter is not None and\
          self.sl_filename is not None:
            sourcelist = (
              (SourceList.OBJECT == self.sl_object) &\
              (SourceList.filter.name == self.sl_filter) &\
              (SourceList.filename.like(self.sl_filename))
            ).max('creation_date')

        elif self.sl_object is not None and self.sl_filter is not None:
            sourcelist = (
              (SourceList.OBJECT == self.sl_object) &\
              (SourceList.filter.name == self.sl_filter)
            ).max('creation_date')

        else:
            print message %\
              (self.sl_object, self.sl_filter, self.sl_filename,
              self.sl_slid)

            raise SystemExit

        message = '\nFound SourceList with name %s and SLID %s.' %\
          (sourcelist.name, str(sourcelist.SLID))
        print message

        # Try to set object attributes here, so that we can use them in 
        # other methods.
        try:
            self.sl_object = sourcelist.OBJECT
            self.sl_filter = sourcelist.filter.name
            self.sl_filename = sourcelist.filename
            self.sl_slid = sourcelist.SLID
        
        except:
            pass
            
        return sourcelist

    def load_data(self, atts):
        """Method to get the necessary columns for a plot or for the
        method find_stars.

        The columns should be in the list argument atts.

        Returns a list of numpy arrays, were each array is for each of
        the necessary columns.

        Should be called by child classes.
        """
        tc = TableConverter()

        tc.load_sourcelist(self.sourcelist, attr_list=atts)

        # Copy only the attributes needed for plotting, in order
        # to save on processing load to the computer. Then store
        # in a list of numpy arrays.
        data = [
            np.array(tc.copy(attributes=atts).data[i]) for i in atts
        ]

        return data

    def print_message(self, status, plot_type, get_stats, do_plot,
      start_time=None):
        """Used twice per plot to print to the terminal."""
        message = "\n%s %s make" % (status, plot_type)

        if get_stats and not do_plot:
            message += " statistics"

        if not get_stats and do_plot:
            message += " plot"

        if status == "Starting":
            return "%s." % message

        try:
            return "%s, took %s." % (message, (datetime.now()-start_time))
        except:
            return

    def bin_mat(self, nx, num_per_bin_min, num_per_bin, x, y, plot=None):
        """Do binned trend, in each bin calculating average of edges,
        median, and standard deviation.

        Implemented from bin_mat subroutine by Naples KiDS-CAT
        (La Barbera et al.).

        Argument plot is set by KidsCatCompl to run additional steps
        for that plot.

        Loop through x, which has length nx, and make bins with
        num_per_bin or at least num_per_bin_min.

        Keep in the xb list the average of the edges of the bin. Keep in
        the yb list the numpy medians of the bins. Keep in the eyb list
        the numpy standard deviation of the median. Count the number of
        sources and store in nb.

        Return xb, yb, nb, eyb.

        Should be called by child classes.
        """
        # Sort the parameters x and y.
        x_sorted = sorted(x)
        y_sorted = [j for (i,j) in sorted(zip(x, y))]

        # Lists to hold values, which will be returned.
        xb, yb, eyb = [], [], []

        # Index of lower edge of the bin.
        i1 = -1 * num_per_bin

        # Counter.
        i = 0

        # Loop over all sources in the parameter x.
        while i1 <= nx:
            i += 1

            # Modify the index of the lower edge to correspond to
            # the desired number per bin.
            i1 += num_per_bin

            # Index of the higher edge of the bin.
            i2 = i1 + num_per_bin - 1

            # Return if looped over all sources or if remaining
            # number of sources is less that the desired minimum
            # number per bin.
            if (i1 >= nx) or ((nx - i1) < num_per_bin_min):

                nb = i

                return xb, yb, nb, eyb

            # In case the while is at the last possible bin.
            i2 = min(i2, nx)

            # Get the value of the lower edge of the bin.
            x1 = x_sorted[i1]

            # Get the value of the higher edge of the bin.
            try:
                x2 = x_sorted[i2]

            except:
                i2 -= 1
                x2 = x_sorted[i2]

            # Append to list that will be returned, the average.
            xb.append(
              (x1 + x2) / 2.
            )

            # Holds values for this bin.
            v = []

            for ii in range(i1, i2+1):

                v.append(
                  y_sorted[ii]
                )

            # Alternative step for the SG plot. Append to eyb the 
            # sigma value calculated by the method bwsm.
            if plot == 'sg':
                yb.append(
                  np.median(v)
                )
                
                theta, sigma = self.__bwsm(v, num_sigma=2.5)
                
                eyb.append(
                  sigma
                )
            
            # Alternative step for the COMPL plot. Append to yb the 
            # 98th-percentile value of the bin.
            elif plot == 'compl':

                v_sorted = sorted(v)

                yb.append(
                  v_sorted[:int(len(v_sorted) * 0.98)][-1]
                )

                eyb.append(np.std(v))
                
            # Append to lists that will be returned, the standard
            # deviation of the bin and the median of the bin.
            else:

                yb.append(np.median(v))

                eyb.append(np.std(v))

    def find_stars(self, plot_type=None):
        """Find sure stars using criteria.

        Criteria is described in the OmegaCen wiki, step 2 at:
        http://wiki.astro-wise.org/projects:kids:data_deliveries:kids-eso-dr1:release_notes:star-galaxy_separation

        The argument plot_type is only set by KidsCatSg as 'SG'.

        Files will be saved to path including environment variable
        $AWEPIPE, or otherwise to path including the current working
        directory. File stars_cand_SL_name.dat will be saved in this
        path joined with astro/experimental/kids/kidscat_temp/. File
        seeing_SL_name.dat will be saved in this path joined with
        astro/experimental/kids/kidscat_output_res/. In original KiDS-CAT
        files were saved in directories TEMP and OUTPUT_RES.

        Return, if plot_type='SG', the indices of the selected sources
        and of the sure stars.

        Return otherwise a dictionary where string keys are fwhm_mean,
        fwhm_sigma, fwhm_mean_stars, and fwhm_sigma_stars.

        Should be called by child classes.
        """
        # SourceList columns needed.
        atts = [
          'MAGERR_AUTO', 'IMAFLAGS_ISO', 'FWHM_IMAGE', 'ELLIPTICITY',
          'CLASS_STAR', 'NPIX', 'Flag', 'SeqNr'
        ]

        # Load the columns.
        data = self.load_data(atts)

        # Keep the indices, which makes writing and reading the code
        # easier.
        indx_magerr_auto = atts.index('MAGERR_AUTO')
        indx_imaflags_iso = atts.index('IMAFLAGS_ISO')
        indx_fwhm_image = atts.index('FWHM_IMAGE')
        indx_ellipticity = atts.index('ELLIPTICITY')
        indx_class_star = atts.index('CLASS_STAR')
        indx_isoarea_image = atts.index('NPIX')
        indx_flag = atts.index('Flag')
        indx_nums = atts.index('SeqNr')

        # Hold values used to calculate the statistics.
        fwhm = []

        # Hold values to return.
        indices_stars_pre = []

        # Value following the description in the wiki.
        ellip_value = 0.1

        if plot_type != 'SG':
            # Prepare for writing new stars_cand_SL_filename.dat
            # file by deleting any previously made.
            dirname = os.path.join(
              self.path,
              'TEMP' 
            )
            
            if self.sourcelist.filename:
                filename = os.path.join(
                  dirname,
                  'stars_cand_%s.dat' %\
                  self.sourcelist.filename.replace('.fits', '')
                )

            else:
                filename = os.path.join(
                  dirname,
                  'stars_cand_%s.dat' %\
                  self.sourcelist.SLID
                )

            if not os.path.isdir(dirname):
                os.makedirs(dirname)

            else:
                if glob.glob(filename):
                    subprocess.call('rm -f %s' % filename, shell=True)

                else:
                    message = '\nFile %s was not found for deletion. It will be appended to.' %\
                      filename
                    print message

        # Find candidate stars.
        # 20150804: had to introduce taking sources with
        # IMAFLAGS_ISO == 128 since some sources kept by FlB have
        # this number == 0 but in AWE it is == 128.
        #
        # TO DO: check if number of sources in file is same as
        # number of sources used by FlB.
        for i in range(len(data[0])):
            if (1/data[indx_magerr_auto][i] >= 50) &\
              (data[indx_isoarea_image][i] > 7) &\
              (data[indx_imaflags_iso][i] == 0 or\
              data[indx_imaflags_iso][i] == 128) &\
              (data[indx_fwhm_image][i] > 0.) &\
              (data[indx_ellipticity][i] < ellip_value) &\
              (data[indx_flag][i] == 0 or data[indx_flag][i] == 16):

                fwhm.append(data[indx_fwhm_image][i])

                if plot_type=='SG':
                    indices_stars_pre.append(i)

                # Write to file created above.
                else:
                    line = '%i %1.4f %1.2f %i\n' % (
                      data[indx_nums][i],
                      data[indx_magerr_auto][i],
                      data[indx_fwhm_image][i],
                      data[indx_flag][i]
                    )

                    with open(filename, 'ab') as fp:
                        fp.write(line)

        # If too few stars with ELLIPTICITY < 0.1, increase to < 0.2.
        # Follows description of plot from the wiki.
        if len(fwhm) < 10:
            message = '\nToo few sure stars with ELLIPTICITY < 0.1; will use ELLIPTICITY < 0.2.'
            print message

            subprocess.call('rm -f %s' % filename, shell=True)

            fwhm = []
            indices_stars_pre = []
            ellip_value = 0.2

            for i in range(len(data[0])):

                if (1/data[indx_magerr_auto][i] >= 50) &\
                  (data[indx_isoarea_image][i] > 7) &\
                  (data[indx_imaflags_iso][i] == 0 or\
                  data[indx_imaflags_iso][i] == 128) &\
                  (data[indx_fwhm_image][i] > 0.) &\
                  (data[indx_ellipticity][i] < ellip_value) &\
                  (data[indx_flag][i] == 0 or data[indx_flag][i] == 16):

                    fwhm.append(data[indx_fwhm_image][i])

                    if plot_type=='SG':
                        indices_stars_pre.append(i)

                    else:
                        line = '%i %1.4f %1.2f %i\n' % (
                          data[indx_nums][i],
                          data[indx_magerr_auto][i],
                          data[indx_fwhm_image][i],
                          data[indx_flag][i]
                        )

                        with open(filename, 'ab') as fp:
                            fp.write(line)

        fwhm_mean = self.__biweight_location(
          fwhm, num_sigma=4
        )

        fwhm_sigma = self.__biweight_scale(
          fwhm, num_sigma=4
        )

        if plot_type != 'SG':
            # Prepare for writing new seeing_SL_filename.dat file by
            # deleting any previously made.
            dirname = os.path.join(
              self.path,
              'OUTPUT_RES' 
            )
            
            if self.sourcelist.filename:
                filename = os.path.join(
                  dirname,
                  'seeing_%s.dat' %\
                  self.sourcelist.filename.replace('.fits', '')
                )

            else:
                filename = os.path.join(
                  dirname,
                  'seeing_%s.dat' %\
                  self.sourcelist.SLID
                )

            if not os.path.isdir(dirname):
                os.makedirs(dirname)

            else:
                if glob.glob(filename):
                    subprocess.call('rm -f %s' % filename, shell=True)

                else:
                    message = '\nFile %s was not found for deletion. It will be appended to.' %\
                      filename
                    print message

            # Keep average FWHM and its width sFWHM. Use scale used
            # by FlB, which is AWE value multiplied by 0.2
            # (conversion between pixels and arcsec).
            line = '%1.7f %1.7f 0.2\n' % (
              fwhm_mean * 0.2, fwhm_sigma * 0.2
            )

            with open(filename, 'wb') as fp:
                fp.write(line)

        # Determine sources within + or - 2.5 times the sFWHM to be
        # sure stars. Follows description of plot from the wiki.
        width_low = -2.5 * fwhm_sigma
        width_high = 2.5 * fwhm_sigma

        indices_stars = []
        fwhm_stars = []

        for i in range(len(data[0])):

            if (width_low <\
              data[indx_fwhm_image][i] - fwhm_mean <\
              width_high) &\
              (1/data[indx_magerr_auto][i] >= 50) &\
              (data[indx_isoarea_image][i] > 7) &\
              (data[indx_imaflags_iso][i] == 0 or\
              data[indx_imaflags_iso][i] == 128) &\
              (data[indx_fwhm_image][i] > 0.) &\
              (data[indx_ellipticity][i] < ellip_value) &\
              (data[indx_flag][i] == 0 or data[indx_flag][i] == 16):

                if plot_type=='SG':
                    indices_stars.append(i)

                fwhm_stars.append(data[indx_fwhm_image][i])

        fwhm_mean_stars = self.__biweight_location(
          fwhm_stars, num_sigma=4
        )

        fwhm_sigma_stars = self.__biweight_scale(
          fwhm_stars, num_sigma=4
        )

        if plot_type=='SG':
            return [indices_stars_pre, indices_stars]

        statistics = {
          'fwhm_mean': fwhm_mean*0.2,
          'fwhm_sigma': fwhm_sigma*0.2,
          'fwhm_mean_stars': fwhm_mean_stars*0.2,
          'fwhm_sigma_stars': fwhm_sigma_stars*0.2
        }

        return statistics

    def get_attribute(self, sl_attribute):
        """Check for attribute. Return it or a replacement"""
        if sl_attribute == 'filename':
            if self.sourcelist.filename:
                return self.sourcelist.filename.replace('.fits', '')

        if sl_attribute == 'name':
            try:
                if self.sourcelist.name:
                    return self.sourcelist.name
            
            except:
                return 'fromCatFile'

        return self.sourcelist.SLID

    def make(self):
        """Overriden in each of the child classes.

        Should only be called by child classes.
        """
        pass

    def make_statistics(self, write_to_file=True):
        """Called by subclasses to only make statistics calculated while
        making a specific plot.

        Boolean variable write_to_file sets if the statistics are also
        written to a file or not.

        The child classes (so each of the plots) extend this method.
        Following are the return statements when called by each of the
        plots.

        Return dictionary with string keys ellip_median_all,
        ellip_median_stars, magsat when called by KidsCatAbmagsat.

        Return dictionary with string keys mcompl when called by
        KidsCatCompl.

        Return dictionary with string keys fwhm_mean, fwhm_sigma,
        fwhm_mean_stars, and fwhm_sigma_stars when called by
        KidsCatFwhmsn.

        Return dictionary with string keys mlim_sn5, mlim_sn10, and
        mlim_sn15 when called by KidsCatMlim.

        Return dictionary with string keys the file names where the
        statistics were saved when called by KidsCatSg.

        Should not be called by KidsCat.
        """
        try:
            statistics = self.make(
              get_stats=True, write_to_file=write_to_file, do_plot=False
            )

            return statistics

        except:
            traceback.print_exc()

            return

    def make_plot(self):
        """Called by subclasses to make only the plot.

        Plot is saved in current working directory.

        Return none.

        Should not be called by KidsCat.
        """
        try:
            self.make(
              get_stats=False, write_to_file=False, do_plot=True
            )

            return

        except:
            traceback.print_exc()

            return
