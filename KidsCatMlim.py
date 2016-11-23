from KidsCat import KidsCatBase

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import glob

# Used here to ignore numpy warning:
# 'Warning: invalid value encountered in log10'
np.seterr(all='ignore')

class KidsCatMlim(KidsCatBase):
    """Make the KiDS-CAT MLIM figure and statistics.

    Implemented from maglim.f by Naples KiDS-CAT (La Barbera et al.)

    Instance variables are inherited from KidsCatBase.

    Public method make overrides KidsCatBase.make.

    Extends from KidsCatBase the methods find_sourcelist, load_data,
    print_message, get_attribute, bin_mat, make_statistics, and
    make_plot.

    All exceptions are handled by printing the traceback and return none.
    """

    atts = [
      'FLUX_APER_10', 'FLUXERR_APER_10', 'NIMAFLAGS_ISO', 'IMAFLAGS_ISO',
      'Flag'
    ]

    def make(self, get_stats=True, write_to_file=True, do_plot=True):
        """Calculate statistics, write them to a file, and plot.

        Arguments determine if statistics are calculated, if files are
        written, and if plots are made, respectively.

        Return none if not get_stats.

        Return otherwise dictionary with statistics. String keys are
        mlim_sn5, mlim_sn10, and mlim_sn15.
        """
        try:
            print self.print_message('Starting', 'MLIM', get_stats, do_plot)

            start_time = datetime.now()

            self.sourcelist = self.find_sourcelist()

            if self.path is None:
                self.path = self.set_path()

            self.data = self.load_data(self.atts)

            # Holds number of sources.
            nx = 0

            # Holds values to be plotted.
            xv1, xv2 = [], []

            for i in range(len(self.data[0])):

                flux_aper = self.data[self.atts.index('FLUX_APER_10')][i]
                mag_aper = -2.5 * np.log10(flux_aper)

                if (flux_aper > 0.) &\
                  (int(self.data[self.atts.index('IMAFLAGS_ISO')][i]) == 0) &\
                  (int(self.data[self.atts.index('Flag')][i]) == 0) &\
                  (mag_aper < 30.):

                    xv1.append(mag_aper)

                    xv2.append(
                      flux_aper / self.data[self.atts.index('FLUXERR_APER_10')][i]
                    )

                    nx += 1

            # Binning the S/N vector according to mag, each bin
            # including the same number of objects num_per_bin, with
            # minimum allowed num_per_bin_min.
            num_per_bin = 30
            num_per_bin_min = 10

            xb, yb, nb, eyb = self.bin_mat(
              nx, num_per_bin_min, num_per_bin, xv1, xv2
            )

            # Use KiDS-CAT method for initial calculation of the
            # limiting magnitude at S/N = 5, which is in turn used to
            # calculate the limits of the plot xmi and xma.
            par_200000, par_10, par_100000 = [], [], []

            self.mlim_sn5, b, c, par_1, par_2, par_3 = 0., 0., 0., 0., 0., 0.

            while self.mlim_sn5 == 0.:

                acc = 1.0e-5

                nmax = 10000

                if self.mlim_sn5 == 0.:
                    mcut = 22.

                else:
                    mcut = self.mlim_sn5 - 4.

                for i in range(len(xb)):

                    if xb[i] > mcut:
                        par_200000.append(1.)

                    else:
                        par_200000.append(0.)

                    par_10.append(xb[i])

                    par_100000.append(yb[i])

                if self.mlim_sn5 == 0.:

                    par_1 = 1.
                    par_2 = 10.**(-0.4 * (par_10[0] - 25.)) / par_100000[0]**2
                    par_3 = par_1 * 10.**(-0.4 * (par_10[-1] - 25.)) / par_100000[-1]
                    par_3 = par_3**2

                snlim = 5.

                b = -(snlim / par_1)**2*par_2
                c = -(snlim / par_1)**2*par_3

                self.mlim_sn5 = (-b + np.sqrt(b**2 - 4.*c)) / 2.
                self.mlim_sn5 = -2.5 * np.log10(self.mlim_sn5) + 25.

            # Declare variables to hold limits for plot.
            # Are also used in current implementation of mlim calculation.
            xmi = self.mlim_sn5 - 3.5
            xma = self.mlim_sn5 + 1.
            ymi = 0.1
            yma = 99.

            # Use numpy's polyfit to fit the binned trend.
            #
            # TO DO:
            # Take a look or implement KiDS-CAT comment below.
            #c only faint objects (>=21) are used in the fit
            xb_fit, yb_fit = [], []

            for i in range(len(xb)):

                if xmi < xb[i] < xma+0.5:

                    xb_fit.append(xb[i])
                    yb_fit.append(yb[i])

            var_polyfit = np.polyfit(xb_fit, yb_fit, 3)

            my_func = np.poly1d(var_polyfit)

            x_polyfit = np.linspace(xmi, xma+0.5, num=1000)

            y_polyfit = my_func(x_polyfit)

            # Create a polynomial using the found coefficients.
            a, b, c, d =\
              var_polyfit[0], var_polyfit[1], var_polyfit[2], var_polyfit[3]

            p = np.poly1d([a, b, c, d])

            # Use the roots and the value of y we want to find x at, i.e.
            # y=5, to find our mlim. See how-to answer from Charles
            # Harris at:
            # http://stackoverflow.com/questions/16827053/solving-for-x-values-of-polynomial-with-known-y-in-scipy-numpy
            #
            # Round to 3 places to match KiDS-CAT.
            self.mlim_sn5 = round(np.real(min((p - 5).r)), 3)

            self.mlim_sn10 = round(np.real(min((p - 10).r)), 3)

            self.mlim_sn15 = round(np.real(min((p - 15).r)), 3)

            # Write to file MLIM_SL_name.dat the limiting magnitude at
            # S/N = 5, 10, and 15.
            dirname = os.path.join(
              self.path,
              'OUTPUT_RES' 
            )

            sl_attribute = self.get_attribute('filename')

            filename = os.path.join(
              dirname, 'MLIM_%s.dat' %\
              sl_attribute
            )

            if not os.path.isdir(dirname):
                os.makedirs(dirname)

            line = '#SN_LIM MLIM\n5.0 %s\n10.0 %s\n15.0 %s\n' % (
              self.mlim_sn5, self.mlim_sn10, self.mlim_sn15
            )

            with open(filename, 'wb') as fp:
                fp.write(line)

            if do_plot:
                plt.ioff()

                fig, ax = plt.subplots(nrows=1, ncols=1)

                # Plot the sources.
                ax.scatter(
                  xv1, xv2,
                  facecolor='k', edgecolor='k', s=1
                )

                # Plot the binned trend.
                ax.plot(
                  xb, yb, 'r', linewidth=1.5
                )

                # Plot the fit to the binned trend.
                ax.plot(
                  x_polyfit, y_polyfit, 'b', linewidth=1.5
                )

                # Plot line at S/N = 5.
                ax.axhline(y=5, color='k', linestyle='--', linewidth=1)

                # Plot the magnitude at S/N = 5.
                ax.axvline(x=self.mlim_sn5, color='k', linewidth=1)

                # Annotate as necessary.
                ax.annotate(
                  'binned trend', xy=(xma-0.6*(xma-xmi), yma-0.1*(yma-ymi)),
                  color='r'
                )

                ax.annotate(
                  'best-fit', xy=(xma-0.6*(xma-xmi),yma-0.17*(yma-ymi)),
                  color='b'
                )

                ax.annotate(
                  r'mag$_{lim}$=%s' % str(round(self.mlim_sn5, 2)),
                  xy=(xma-0.6*(xma-xmi),yma-0.24*(yma-ymi)),
                  color='k'
                )

                plt.xlabel('mag (2" diameter)')
                plt.ylabel('S/N')
                plt.xlim(xmi, xma)
                plt.ylim(ymi, yma)

                sl_attribute = self.get_attribute('name')

                plt.title('%s' % sl_attribute)

                filename_plot = os.path.join(
                  self.path,
                  'PLOTS',
                  '%s-%s.png' %\
                  (sl_attribute, 'MLIM')
                )

                # Does not work in AWE.
                #plt.tight_layout()

                fig.savefig(filename_plot)

                if glob.glob(filename_plot):
                    message = '\nPlot saved as %s.' % filename_plot
                    print message

            print self.print_message(
              'Finished', 'MLIM', get_stats, do_plot, start_time
            )

            if not get_stats:
                return

            message = '\nMLIM at S/N=5: %s\nMLIM at S/N=10: %s\nMLIM at S/N=15: %s' %\
              (self.mlim_sn5, self.mlim_sn10, self.mlim_sn15)
            print message

            statistics = {
              'mlim_sn5': self.mlim_sn5,
              'mlim_sn10': self.mlim_sn10,
              'mlim_sn15': self.mlim_sn15,
            }

            return statistics

        except:
            import traceback

            traceback.print_exc()

            return
