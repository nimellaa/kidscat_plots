from KidsCat import KidsCatBase

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import glob
import os

class KidsCatCompl(KidsCatBase):
    """Make the KiDS-CAT COMPL figure and statistics.

    Implemented from completeness.f by Naples KiDS-CAT
    (La Barbera et al.)

    Instance variables are inherited from KidsCatBase.

    Public method make overrides KidsCatBase.make.

    Extends from KidsCatBase the methods find_sourcelist, load_data,
    print_message, get_attribute, bin_mat, make_statistics, and
    make_plot.

    All exceptions are handled by printing the traceback and return none.
    """

    atts = [
      'FLUX_APER_2', 'FLUXERR_APER_2', 'MAG_AUTO', 'MAGERR_AUTO', '2DPHOT',
      'IMAFLAGS_ISO', 'Flag'
    ]

    def make(self, get_stats=True, write_to_file=True, do_plot=True):
        """Calculate statistics, write them to a file, and plot.

        Arguments determine if statistics are calculated, if files are
        written, and if plots are made, respectively.

        Return none if not get_stats.

        Return otherwise dictionary with statistic. String key is mcompl.
        """
        try:
            print self.print_message('Starting', 'COMPL', get_stats, do_plot)

            start_time = datetime.now()

            self.sourcelist = self.find_sourcelist()

            if self.path is None:
                self.path = self.set_path()

            self.data = self.load_data(self.atts)

            # Keep track of number of sources for sure stars, sources
            # within the criteria, and not stars outside the criteria.
            nx, nss, ni = 0, 0, 0

            # Hold values used in plotting.
            xv1, xv2, xss, yss, xi, yi = [], [], [], [], [], []

            self.atts = [
              'FLUX_APER_2', 'FLUXERR_APER_2', 'MAG_AUTO', 'MAGERR_AUTO', '2DPHOT',
              'IMAFLAGS_ISO', 'Flag'
            ]

            for i in range(len(self.data[0])):

                flux_aper_2 = self.data[self.atts.index('FLUX_APER_2')][i]
                mag_aper_2 = -2.5 * np.log10(flux_aper_2)

                mag_auto = self.data[self.atts.index('MAG_AUTO')][i]

                if (flux_aper_2 > 0.) &\
                  (int(self.data[self.atts.index('IMAFLAGS_ISO')][i]) == 0.) &\
                  (int(self.data[self.atts.index('Flag')][i]) == 0.) &\
                  (mag_aper_2 < 30.) &\
                  (int(self.data[self.atts.index('2DPHOT')][i]) >= 0.):

                    if ((1.0857/self.data[self.atts.index('MAGERR_AUTO')][i]) > 2.) &\
                      ((flux_aper_2/self.data[self.atts.index('FLUXERR_APER_2')][i]) > 2.7):

                        # Store the black points, i.e. sources within
                        # the desired criteria and not sure stars.
                        nx += 1
                        xv1.append(mag_aper_2)
                        xv2.append(mag_auto)

                        # Store the sure stars.
                        if int(self.data[self.atts.index('2DPHOT')][i]) >= 4:

                           nss += 1
                           xss.append(mag_aper_2)
                           yss.append(mag_auto)

                    # Store all other sources.
                    else:

                        ni += 1
                        xi.append(mag_aper_2)
                        yi.append(mag_auto)

            # Set the limits for the plot, stored in variables xma, xmi,
            # yma, and ymi.
            xv1_sorted = sorted(xv1)

            xma = xv1_sorted[:int(len(xv1_sorted) * .98)][-1] + 0.5

            xmi = xma - 5.

            xv2_sorted = sorted(xv2)

            yma = xv2_sorted[:int(len(xv2_sorted) * .98)][-1] + 0.5

            ymi = yma - 5.

            # Binning the class/star vector for sure stars according to
            # S/N, each bin including the same number of objects
            # num_per_bin, with minimum allowed number num_per_bin_min.
            num_per_bin = int(0.03 * float(nss))
            num_per_bin = max(num_per_bin, 30)

            num_per_bin_min = max(1, num_per_bin/2.)

            xb, yb, nb, eyb = self.bin_mat(
              nx, num_per_bin_min, num_per_bin, xv2, xv1, plot='compl'
            )

            yb_sorted = sorted(yb)

            yb_limit = yb_sorted[:int(len(yb_sorted) * 0.8)][-1]

            # Calculate the MAG_AUTO median of upper 80% of sources.
            v = []

            for i in range(len(yb)):

                if yb[i] > yb_limit:
                    v.append(yb[i])

            t = np.median(v)
            s = np.std(v)

            xv, yv = [], []

            for i in range(len(xb)):

                if yb[i] < (t - 5.*s - 0.3):

                    if yb[i] > ((t - 5.*s) - 2.3):

                       xv.append(yb[i])

                       yv.append(xb[i])

            var_polyfit = np.polyfit(xv, yv, 1)

            my_func = np.poly1d(var_polyfit)

            x_polyfit = np.linspace(xmi, xma, num=1000)
            y_polyfit = my_func(x_polyfit)

            # Get the coefficients A, B of the fit and calculate mcompl
            # at the median t with equation y = A * x + B.
            a, b = var_polyfit[0], var_polyfit[1]

            self.mcompl = a*t + b

            if do_plot:
                plt.ioff()

                fig, ax = plt.subplots(nrows=1, ncols=1)

                # Plot not sure stars and sources outside criteria.
                ax.scatter(
                  xi, yi,
                  facecolor='0.7', edgecolor='0.7', s=5
                )

                # Plot not sure stars and sources within criteria.
                ax.scatter(
                  xv1, xv2,
                  facecolor='k', edgecolor='k', s=1
                )

                # Plot sure stars.
                ax.scatter(
                  xss, yss,
                  facecolor='b', edgecolor='b', s=5,
                )

                # Plot binned means.
                ax.plot(
                  yb, xb, 'r', linewidth=1.5
                )

                # Plot fit to binned means.
                ax.plot(x_polyfit, y_polyfit, 'c', linewidth=1.5)

                # Plot mcompl and annotate it.
                ax.axhline(y=self.mcompl, color='c', linestyle='--')
                ax.annotate(
                  '98% COMPLETENESS', xy=(xmi+0.05*(xma-xmi),
                  self.mcompl+0.01*(yma-ymi)),
                  color='c'
                )
                ax.annotate(
                  'MAG=%2.2f' % self.mcompl, xy=(xmi+0.05*(xma-xmi),
                  self.mcompl-0.04*(yma-ymi)),
                  color='c'
                )

                # Plot the median of upper 80% of sources.
                ax.axvline(x=t, color='r', linestyle='--')

                # Annotate legend explaining what the blue points are.
                ax.annotate(
                  'STARS', xy=(xmi+0.05*(xma-xmi), ymi+0.2*(yma-ymi)),
                  color='b'
                )

                plt.xlabel('DETECTION MAG')
                plt.ylabel('MAG AUTO')
                plt.xlim(xmi, xma)
                plt.ylim(ymi, yma)

                sl_attribute = self.get_attribute('name')

                plt.title('%s' % sl_attribute)

                filename_plot = os.path.join(
                  self.path,
                  'PLOTS',
                  '%s-%s.png' %\
                  (sl_attribute, 'COMPL')
                )

                # Does not work in AWE.
                #plt.tight_layout()

                fig.savefig(filename_plot)

                if glob.glob(filename_plot):
                    message = '\nPlot saved as %s.' % filename_plot
                    print message

            print self.print_message(
              'Finished', 'COMPL', get_stats, do_plot, start_time
            )

            if not get_stats:
                return

            message = '\nMCOMPL: %s' % self.mcompl
            print message

            statistics = {'mcompl': self.mcompl}

            return statistics

        except:
            import traceback

            traceback.print_exc()

            return

