from KidsCat import KidsCatBase

import numpy as np
from datetime import datetime
import os
import glob
import subprocess
import matplotlib.pyplot as plt
import sys

# Used here to ignore numpy warning:
# 'Warning: invalid value encountered in log10'
np.seterr(all='ignore')

class KidsCatSg(KidsCatBase):
    """Make the KiDS-CAT SG figure and statistics.

    Implemented from stargal.f by Naples KiDS-CAT (La Barbera et al.)

    Instance variables are inherited from KidsCatBase.

    Public method make overrides KidsCatBase.make.

    Extends from KidsCatBase the methods find_sourcelist, load_data,
    print_message, find_stars, bin_mat, get_attribute, make_statistics,
    and make_plot.

    All exceptions are handled by printing the traceback and return none.
    """

    atts = [
      'FLUX_APER_10', 'IMAFLAGS_ISO', 'MAGERR_AUTO', 'CLASS_STAR', '2DPHOT',
      'Flag'
    ]

    def __do_get_stats(self, x, y):
        """Write statistics to file.

        SG_CURVE_SL_name.dat contains the S/N and corresponding
        CLASS_STAR. SG_CLASS_SL.name.dat contains the CLASS_STAR of
        some sources.

        Files will be saved to path including environment variable
        $AWEPIPE, or otherwise to path including the current working
        directory. This is joined with
        astro/experimental/kids/kidscat_output_res/. In original
        KiDS-CAT files were saved in directory OUTPUT_RES.

        Return list with file names.

        Should be called only by KidsCatSg make methods.
        """
        indices_stars_pre, indices_stars =\
          self.find_stars(plot_type='SG')

        dirname = os.path.join(
          self.path,
          'OUTPUT_RES' 
        )

        sl_attribute = self.get_attribute('filename')

        filename_curve = os.path.join(
          dirname, 'SG_CURVE_%s.dat' %\
          sl_attribute
        )

        filename_class = os.path.join(
          dirname, 'SG_CLASS_%s.dat' %\
          sl_attribute
        )

        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        if glob.glob(filename_class):
            subprocess.call('rm -f %s' % filename_class, shell=True)

        else:
            message = '\nFile %s was not found for deletion. It will be appended to.' %\
              filename_class
            print message

        # Write to SG_CURVE_SL_name.dat.
        with open(filename_curve, 'wb') as fp:
            fp.write('#S/N CLASS_STAR\n')

            for i in range(len(x)):

                fp.write(
                  '%4.4f %1.3f\n' % (10.**(x[i]), y[i])
                )

        # Append to SG_CLASS_SL_name.dat.
        for i in range(len(self.data[0])):

            dmin = 1.

            lmin = 0

            for j in range(len(x)):

                dist = np.abs(
                  np.log10(
                    1.0857 / self.data[self.atts.index('MAGERR_AUTO')][i]
                  ) - x[j]
                )

                if dist < dmin:
                   dmin = dist
                   lmin = j

            if self.data[self.atts.index('CLASS_STAR')][i] >= y[lmin]:
                class_flag = 4

            else:
                class_flag = 0

            twod_phot = int(self.data[self.atts.index('2DPHOT')][i])

            if int(class_flag + twod_phot) == 8:
                class_flag = 4

            elif int(class_flag + twod_phot) == 9:
                class_flag = 5

            if i in indices_stars_pre:
                class_flag = 4

            if i in indices_stars:
                class_flag = 5

            else:
                class_flag = twod_phot

            with open(filename_class, 'ab') as fp:
                fp.write(
                  '%i\n' % class_flag
                )

        return [filename_curve, filename_class]

    def make(self, get_stats=True, write_to_file=True, do_plot=True):
        """Calculate statistics, write them to a file, and plot.

        Arguments determine if statistics are calculated, if files are
        written, and if plots are made, respectively.

        Return none if not get_stats.

        Return otherwise dictionary with statistics. String keys are
        'filename SG curve' and 'filename SG class'.
        """
        try:
            print self.print_message('Starting', 'SG', get_stats, do_plot)

            start_time = datetime.now()

            self.sourcelist = self.find_sourcelist()

            if self.path is None:
                self.path = self.set_path()
            
            self.data = self.load_data(
                self.atts
            )
                
            # Hold number of sources.
            nss, nss2 = 0, 0

            # Hold values that will be plotted.
            xv1, xv2, xss, yss, xss2, yss2 = [], [], [], [], [], []

            #import numpy as np

            for i in range(len(self.data[0])):

                flux_aper = self.data[self.atts.index('FLUX_APER_10')][i]

                magerr_auto = self.data[self.atts.index('MAGERR_AUTO')][i]

                class_star = self.data[self.atts.index('CLASS_STAR')][i]

                twod_phot = self.data[self.atts.index('2DPHOT')][i]

                try:
                    indx_flag = self.atts.index('Flag')
                except:
                    indx_flag = self.atts.index('FLAGS')
                
                #import numpy as np
                
                if (flux_aper > 0.) &\
                  (int(self.data[self.atts.index('IMAFLAGS_ISO')][i]) == 0) &\
                  (int(self.data[indx_flag][i]) == 0) &\
                  ((-2.5 * np.log10(flux_aper)) < 30.):

                    xv1.append(
                      10.**(np.log10(1.0857 / magerr_auto) +\
                      np.random.normal(0., 0.01))
                    )

                    xv2.append(
                      class_star + np.random.normal(0., 0.005)
                    )

                    if (int(twod_phot) == 5) or (int(twod_phot) == 1):

                        nss += 1

                        xss.append(xv1[-1])

                        yss.append(xv2[-1])

                    if class_star > 0.8:
                        if (1. / magerr_auto) < 50.:
                            nss2 += 1

                            xss2.append(
                              10.**(np.log10(1.0857 / magerr_auto) +\
                              np.random.normal(0., 0.01))
                            )

                            yss2.append(
                              class_star +\
                              np.random.normal(0., 0.005)
                            )

            num_per_bin = int(0.02 * float(nss))

            num_per_bin_min = max(1, num_per_bin/2.)

            xb, yb, nb, eyb = self.bin_mat(
              nss, num_per_bin_min, num_per_bin, xss, yss, plot='sg'
            )

            xbs, ybs, nbs, eybs = self.bin_mat(
              nss2, num_per_bin_min, num_per_bin, xss2, yss2, plot='sg'
            )

            xv1 = np.log10(np.array(xv1))

            xss = np.log10(np.array(xss))

            if do_plot:
                plt.ioff()

                fig, ax = plt.subplots(nrows=1, ncols=1)

                # Plot not sure stars.
                ax.scatter(
                  xv1, xv2,
                  facecolor='k', edgecolor='k', s=1
                )

                # Plot sure stars.
                ax.scatter(
                  xss, yss,
                  facecolor='r', edgecolor='r', s=3
                )

            # Use info from binned trend to plot it.

            xx = np.log10(np.array(xb))
            
            yy = [j - 3.*eyb[i] - 0.002 for i, j in enumerate(yb)]
            yy = [0.75 if i <= 0.75 else i for i in list(yy)]

            xx2 = np.log10(np.array(xbs))
            
            yy2 = [
              j - 2.*eybs[i] - (ybs[-1] - 2.*eybs[i]) + yy[0] for i, j in enumerate(ybs)
            ]
            yy2 = [0.8 if i <= 0.8 else i for i in list(yy2)]

            # Plot the binned trend.
            xx_final = np.concatenate((xx2, xx))
            yy_final = yy2 + yy

            if do_plot:
                ax.plot(
                  xx_final, yy_final, 'b', linewidth=1.5
                )

            # Run at this point, since following KiDS-CAT method the
            # lists xx_final and yy_final that we need are the ones created just
            # above.
            if get_stats:

                filename_curve, filename_class = self.__do_get_stats(xx_final, yy_final)

            # Set limits for the plot.
            xmi = np.log10(4000.)
            xma = np.log10(5.)
            ymi = -0.01
            yma = 1.02

            if do_plot:
                plt.xlabel('S/N')
                plt.ylabel('CLASS STAR')
                plt.xscale('log')
                plt.xlim(xmi, xma)
                plt.ylim(ymi, yma)

                sl_attribute = self.get_attribute('name')

                plt.title('%s' % sl_attribute)

                filename_plot = os.path.join(
                  self.path,
                  'PLOTS',
                  '%s-%s.png' %\
                  (sl_attribute, 'SG')
                )

                # Does not work in AWE.
                #plt.tight_layout()

                fig.savefig(filename_plot)

                if glob.glob(filename_plot):
                    message = '\nPlot saved as %s.' % filename_plot
                    print message

            print self.print_message(
              'Finished', 'SG', get_stats, do_plot, start_time
            )

            if not get_stats:
                return

            message = '\nStatistics written to files\n%s\n%s' %\
              (filename_curve, filename_class)
            print message

            filenames = {
              'filename SG curve': filename_curve,
              'filename SG class': filename_class
            }

            return filenames

        except:
            import traceback

            traceback.print_exc()

            return
