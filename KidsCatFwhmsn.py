from KidsCat import KidsCatBase

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import glob

class KidsCatFwhmsn(KidsCatBase):
    """Make the KiDS-CAT FWHMSN figure and statistics.

    Implemented from plot_ss.f by Naples KiDS-CAT (La Barbera et al.)

    Instance variables are inherited from KidsCatBase.

    Public method make overrides KidsCatBase.make.

    Extends methods KidsCatBase.make_statistics and
    KidsCatBase.make_plot.

    Extends from KidsCatBase the methods find_sourcelist, load_data,
    print_message, get_attribute, find_stars, make_statistics, and
    make_plot.

    Uses the file stars_cand_SL_name.dat, which is made by the function
    find_stars.

    All exceptions are handled by printing the traceback and return none.
    """

    atts = [
      'MAGERR_AUTO', 'FWHM_IMAGE', '2DPHOT', 'Flag', 'SeqNr'
    ]

    def __flb_flag(self, ifll):
        """Implemented from subroutine flag in plot_ss.f by Naples
        KiDS-CAT (La Barbera et al.).

        If flag value is not one, in the dictionary set the value
        to 1 corresponding to the key equal to the flag value. Append
        the values in reverse order to a list.

        Return list of integers.

        20150824: Not used, but keep in case we find out it is better
        to use it.
        
        Should be called only by KidsCatFwhmsn make methods.
        """
        # Original method. Keep here for reference.
        #
        # Check the Flag value in powers of 2 from 128 to 1. If
        # value is larger than current power of 2, append 1 to list
        # iflag and substract current power of 2 from the value;
        # otherwise append 0 to iflag. Repeat until 2**0.
        #iflag = []
        #
        #for j in range(8, 0, -1):

            #if ifll >= 2**(j-1):

                #ifll = ifll - 2**(j-1)

                #iflag.append(1)
            #else:
                #iflag.append(0)

        #return iflag

        # My implementation.
        dic = {1: 0, 2: 0, 4: 0, 8:0, 16: 0, 32: 0, 64:0, 128:0}

        if ifll != 0:
            dic[ifll] = 1

        flags = []

        for key in sorted(dic.keys(), reverse=True):
            flags.append(dic[key])

        return flags

    def make(self, get_stats=True, write_to_file=True, do_plot=True):
        """Calculate statistics, write them to a file, and plot.

        Arguments determine if statistics are calculated, if files are
        written, and if plots are made, respectively.

        Return none if not get_stats.

        Return otherwise dictionary with statistics. String keys are
        fwhm_mean, fwhm_sigma, fwhm_mean_stars, and fwhm_sigma_stars.
        """
        try:
            print self.print_message('Starting', 'FWHMSN', get_stats, do_plot)

            start_time = datetime.now()

            self.sourcelist = self.find_sourcelist()

            if self.path is None:
                self.path = self.set_path()

            self.data = self.load_data(self.atts)

            # Keep statistics that will be returned and set instance
            # variables for convenience.
            statistics = self.find_stars()
            self.fwhm_mean = statistics['fwhm_mean']
            self.fwhm_sigma = statistics['fwhm_sigma']
            self.fwhm_mean_stars = statistics['fwhm_mean_stars']
            self.fwhm_sigma_stars = statistics['fwhm_sigma_stars']

            # Get seeing values, used in the sFWHM calculation and to
            # set limits for the plot.
            self.find_stars()

            sl_attribute = self.get_attribute('filename')

            filename = os.path.join(
              self.path,
              'OUTPUT_RES',
              'seeing_%s.dat' %\
              sl_attribute
            )

            with open(filename, 'r') as fp:
                seeing = fp.readlines()[0].split()

            seeing_a1 = round(float(seeing[0]), 4)
            seeing_a2 = round(float(seeing[1]), 4)
            seeing_a3 = float(seeing[2])

            aa1 = seeing_a1
            y1 = (seeing_a1 - 1.5*seeing_a2) / seeing_a3
            y2 = (seeing_a1 + 1.5*seeing_a2) / seeing_a3

            ymi = max((seeing_a1 - 5.*seeing_a2) / seeing_a3, 0.)

            xma = -100.
            yma = -100.

            # Get info from stars from file. Put into a1 list.
            #
            # 20150804: Should be read from file KIDS_OBJECT_stars.lis,
            # but cannot find how that file is made. Tried making the
            # same file from the final stars chosen by the function
            # find_stars, but get additional sources with SeqNr ['1922',
            # '4648', '11797', '13699', '14280', '15263', '16564',
            # '17745', '18220']
            a1 = []

            sl_attribute = self.get_attribute('filename')

            filename = os.path.join(
              self.path,
              'TEMP',
              'stars_cand_%s.dat' %\
              sl_attribute
            )

            with open(filename, 'rb') as fp:
                while 1:
                    curr = fp.readline().split()
                    if not curr:
                        break

                    a1.append(int(curr[0]))

            # Hold values that will be plotted
            x, y, xv_final, yv_final = [], [], [], []

            self.atts = [
              'MAGERR_AUTO', 'FWHM_IMAGE', '2DPHOT', 'Flag', 'SeqNr'
            ]
            for i in range(len(self.data[0])):

                num = self.data[self.atts.index('SeqNr')][i]

                mag_auto_err = self.data[self.atts.index('MAGERR_AUTO')][i]

                x.append(mag_auto_err)

                fwhm = self.data[self.atts.index('FWHM_IMAGE')][i]

                y.append(fwhm)

                if (num in a1) &\
                  (self.data[self.atts.index('Flag')][i] != 32) &\
                  (self.data[self.atts.index('2DPHOT')][i] == 5):

                    xv_final.append(mag_auto_err)
                    yv_final.append(fwhm)

                # Previous implementation. Keep for reference.
                #
                # The list a1 contains the NUMBER column for a set of
                # sources, which should be the candidate stars. Find
                # that number in the same column in AWE, which is num
                # created above. Use the number to call flb_flag, check
                # the result and the 2DPHOT column, and if the source
                # passes then append it to be plotted as a sure star.
                #if num in a1:
                    #ifll = flag
                    #iflag = self.__flb_flag(ifll)
                    #if iflag[2] == 0:
                        #if self.data[self.atts.index('2DPHOT')][i] == 5:
                            #xv_final.append(mag_auto_err)
                            #yv_final.append(fwhm)

            if (1. / min(x)) > xma:
                xma = 1. / min(x)

            if max(y) > yma:
                yma = max(y)

            x = -np.log10(list(x))
            xv_final = -np.log10(list(xv_final))

            # Set limits for the plot.
            xmi = np.log10(5.)
            xma = np.log10(min(xma+100., 5000.))
            yma = min((aa1+15.*seeing_a2)/seeing_a3, yma)

            if do_plot:
                plt.ioff()

                fig, ax = plt.subplots(nrows=1, ncols=1)

                # Plot not sure stars.
                ax.scatter(
                  x, y,
                  facecolor='k', edgecolor='k', s=10
                )

                # Plot sure stars.
                ax.scatter(
                  xv_final, yv_final,
                  facecolor='red', edgecolor='red', s=5
                )

                # Plot + and - 1.5 times the sFWHM.
                ax.axhline(y=y1, color='b', linestyle='--')
                ax.axhline(y=y2, color='b', linestyle='--')

                plt.xscale('log')
                plt.xlabel('S/N ratio')
                plt.ylabel('FWHM [pixels]')
                plt.xlim(xmi, xma)
                plt.ylim(ymi, yma)

                sl_attribute = self.get_attribute('name')

                plt.title('%s' % sl_attribute)

                filename_plot = os.path.join(
                  self.path,
                  'PLOTS',
                  '%s-%s.png' %\
                  (sl_attribute, 'FWHMSN')
                )
                
                # Does not work in AWE.
                #plt.tight_layout()

                fig.savefig(filename_plot)

                if glob.glob(filename_plot):
                    message = '\nPlot saved as %s.' % filename_plot
                    print message

            print self.print_message(
              'Finished', 'FWHMSN', get_stats, do_plot, start_time
            )

            if not get_stats:
                return

            message = '\nall FWHM mean: %s\nall FWHM sigma: %s\nstars FWHM mean: %s\nstars FWHM sigma: %s' %\
              (self.fwhm_mean, self.fwhm_sigma, self.fwhm_mean_stars,
              self.fwhm_sigma_stars)
            print message

            return statistics

        except:
            import traceback

            traceback.print_exc()

            return
