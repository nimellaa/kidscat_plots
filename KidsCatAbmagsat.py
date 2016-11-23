from KidsCat import KidsCatBase

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os
import glob

# Used here to ignore numpy warning:
# 'Warning: invalid value encountered in log10'
np.seterr(all='ignore')

class KidsCatAbmagsat(KidsCatBase):
    """Make the KiDS-CAT ABMAGSAT figure and statistics.

    Implemented from satur_mag.f by Naples KiDS-CAT (La Barbera et al.)

    Instance variables are inherited from KidsCatBase.

    Public method make overrides KidsCatBase.make.

    Extends from KidsCatBase the methods find_sourcelist, load_data,
    print_message, find_stars, get_attribute, make_statistics, and
    make_plot.

    Uses the file seeing_SL_name.dat, which is made by the method
    find_stars.

    All exceptions are handled by printing the traceback and return none.
    """

    atts = [
      'MAG_AUTO', 'FWHM_IMAGE', 'FLUX_APER_10', 'IMAFLAGS_ISO', '2DPHOT',
      'ELLIPTICITY', 'MaxVal'
    ]

    def __set_fwhm(self, get_stats):
        """Call function to make file, then read it and calculate FWHM.

        The argument get_stats is a boolean to do or not do statistics.

        Return float FWHM.

        Should be called only by KidsCatAbmagsat make methods.
        """
        sl_attribute = self.get_attribute('filename')

        self.find_stars()

        dirname = os.path.join(
          self.path,
          'OUTPUT_RES' 
        )
        
        filename = os.path.join(
          dirname,
          'seeing_%s.dat' %\
          sl_attribute
        )

        with open(filename, 'r') as fp:
            seeing = fp.readlines()[0].split()

        fwhm = float(seeing[0]) / float(seeing[2])

        # It is not the intention that this class be used to get
        # these statistics, but include this step for convenience.
        if get_stats:
            statistics = self.find_stars()
            self.fwhm_mean = statistics['fwhm_mean']
            self.fwhm_sigma = statistics['fwhm_sigma']
            self.fwhm_mean_stars = statistics['fwhm_mean_stars']
            self.fwhm_sigma_stars = statistics['fwhm_sigma_stars']

        return fwhm

    def make(self, get_stats=True, write_to_file=True, do_plot=True):
        """Calculate statistics, write them to a file, and plot.

        Arguments determine if statistics are calculated, if files are
        written, and if plots are made, respectively.

        Return none if not get_stats.

        Return otherwise dictionary with statistics. String keys are
        ellip_median_all, ellip_median_stars, and magsat.
        """
        try:
            print self.print_message(
              'Starting', 'ABMAGSAT', get_stats, do_plot
            )

            # Statistics can only be written to file if they are
            # calculated, so prevent a discrepancy with the following.
            if not get_stats:
                write_to_file = False

            start_time = datetime.now()

            self.sourcelist = self.find_sourcelist()

            if self.path is None:
                self.path = self.set_path()

            self.fwhm = self.__set_fwhm(get_stats)

            self.data = self.load_data(self.atts)

            # Set limits for the plots. Values are changed in loop below.
            xmi = 100000.
            ymi = 100000.
            xma = -100000.
            yma = -100000.

            # If get_stats is True, these are lists to hold values for
            # calculating the medians of the ellipticities of all
            # sources, elall, and sure stars, elss.
            if get_stats:
                elall, elss = [], []

            # Hold values to plot.
            xv1, xv2 = [], []

            # Hold points that should be fitted to, later fitted with
            # numpy's polyfit.
            xv1_fit, xv2_fit = [], []

            for i in range(len(self.data[0])):

                flux_aper_10 = self.data[self.atts.index('FLUX_APER_10')][i]
                mag_aper_10 = -2.5 * np.log10(flux_aper_10)

                flux_max = self.data[self.atts.index('MaxVal')][i]
                surf_mag = -2.5 * np.log10(flux_max)

                if (self.data[self.atts.index('MAG_AUTO')][i] < 23.) &\
                  (self.data[self.atts.index('FWHM_IMAGE')][i] < (self.fwhm * 1.5)) &\
                  (flux_aper_10 > 0.) & (flux_max > 0.):

                    xv1.append(mag_aper_10)

                    xv2.append(surf_mag)

                    xmi = min(xmi, xv1[-1])
                    ymi = min(ymi, xv2[-1])
                    xma = max(xma, xv1[-1])
                    yma = max(yma, xv2[-1])

                if get_stats:

                    ima_flag = int(
                      self.data[self.atts.index('IMAFLAGS_ISO')][i]
                    )

                    twod_phot = int(
                      self.data[self.atts.index('2DPHOT')][i]
                    )

                    if (ima_flag == 0) &\
                      ((twod_phot == 1) or (twod_phot == 4) or\
                      (twod_phot == 5)):

                        elall.append(
                          self.data[self.atts.index('ELLIPTICITY')][i]
                        )

                    if (ima_flag == 0) &\
                      ((twod_phot == 1) or (twod_phot == 5)):

                        elss.append(
                          self.data[self.atts.index('ELLIPTICITY')][i]
                        )
            
            dirname = os.path.join(
              self.path,
              'OUTPUT_RES' 
            )
            
            # Set file names to save statistics.
            sl_attribute = self.get_attribute('filename')

            filename_elliptic = os.path.join(
              dirname,
              'elliptic_%s.dat' %\
              sl_attribute
            )
            
            filename_abmagsat = os.path.join(
              dirname,
              'abmagsat_%s.dat' %\
              sl_attribute
            )
        
            # Write the medians of elall and elss to 
            # elliptic_SL_name.dat, file with SourceList's filename.
            if get_stats:
                self.ellip_median_all = np.median(elall)
                self.ellip_median_stars = np.median(elss)

            if write_to_file:
                with open(filename_elliptic, 'wb') as fp:
                    fp.write(
                      '%1.4f %1.4f' %\
                      (self.ellip_median_all, self.ellip_median_stars)
                    )

            # Set the msat value from the SourceList's SExtractor
            # configuration value SATUR_LEVEL.
            self.msat = -2.5 * np.log10(self.sourcelist.sexconf.SATUR_LEVEL)

            # Update limits for the plot.
            xmi -= 2.
            ymi = self.msat - 1.

            # Get points that should be fitted to, then fit with numpy's
            # polyfit.
            xv1_fit = [i for i in xv1 if (xmi < i < 19)]
            xv2_fit = [xv2[i] for i,j in enumerate(xv1) if (xmi < j < 19)]

            # Fit the points.
            var_polyfit = np.polyfit(xv1_fit, xv2_fit, 1)

            my_func = np.poly1d(var_polyfit)

            x_polyfit = np.linspace(xmi, xma, num=1000)
            y_polyfit = my_func(x_polyfit)

            # Determine magsat at msat value with coefficients A, B
            # given by fit in equation magsat = (msat - B) / A.
            self.magsat = (self.msat - var_polyfit[1]) / var_polyfit[0]

            # Keep self.magsat value in abmagsat_SL_name.dat.
            if write_to_file:
                with open(filename_abmagsat, 'w') as fp:
                    fp.write('%2.2f' % self.magsat)

            if do_plot:
                plt.ioff()

                fig, ax = plt.subplots(nrows=1, ncols=1)

                # Plot the sources.
                ax.scatter(
                  xv1, xv2,
                  facecolor='k', edgecolor='k', s=1
                )

                # Plot the fit.
                ax.plot(x_polyfit, y_polyfit, 'm-')

                # Plot the msat and annotate it.
                ax.axhline(y=self.msat, color='k', linestyle='--')
                ax.annotate(
                  r'$\mu$ saturation', xy=(xma-0.5*(xma-xmi), self.msat+0.1),
                  color='k'
                )

                # Plot the magsat (ABMAGSAT) and annotate it.
                ax.axvline(x=self.magsat, color='r', linestyle='--')
                ax.annotate(
                  'ABMAGSAT',
                  xy=(self.magsat+0.01*(xma-xmi),yma-0.1*(yma-ymi)),
                  color='r'
                )

                plt.xlabel('mag (2"diameter)')
                plt.ylabel(r'$\mu_{max}$ [mag/arcsec^2]')
                plt.xlim(xmi, xma)
                plt.ylim(ymi, yma)

                sl_attribute = self.get_attribute('name')

                plt.title('%s' % sl_attribute)

                filename_plot = os.path.join(
                  self.path,
                  'PLOTS',
                  '%s-%s.png' %\
                  (sl_attribute, 'ABMAGSAT')
                )

                # Does not work in AWE.
                #plt.tight_layout()
                
                fig.savefig(filename_plot)

                if glob.glob(filename_plot):
                    message = '\nPlot saved as %s.' % filename_plot
                    print message

            print self.print_message(
              'Finished', 'ABMAGSAT', get_stats, do_plot, start_time
            )

            if not get_stats:
                return

            message = '\nEllipticity median of\nall sources: %s\nsure stars: %s' %\
              (self.ellip_median_all, self.ellip_median_stars)
            print message

            message = '\nABMAGSAT: %s' % self.magsat
            print message

            statistics = {
              'ellip_median_all': self.ellip_median_all,
              'ellip_median_stars': self.ellip_median_stars,
              'magsat': self.magsat
            }

            return statistics

        except:
            import traceback

            traceback.print_exc()

            return
