"""
A simple class to read in the reduced polarised-differential VAMPIRES data from the IDL pipeline output,
and do some plots.

Get a data object with
    diffData = vampDiffdata(dataFilename, cubeInfoFilename)
Which has the following data attributes:
    .vhhv       - polarised differential visibilities for instrument Stokes Q
    .vhhverr    - 1-sigma errors for polarised differential visibilities for instrument Stokes Q
    .vhhvu      - polarised differential visibilities for instrument Stokes U
    .vhhvuerr   - 1-sigma errors for polarised differential visibilities for instrument Stokes U
    .diffCP     - polarised differential closure-phases for instrument Stokes Q
    .diffCPerr  - 1-sigma errors for polarised differential closure-phases for instrument Stokes Q
    .diffCPu    - polarised differential closure-phases for instrument Stokes U
    .diffCPuerr - 1-sigma errors for polarised differential closure-phases for instrument Stokes U
    .blengths   - baseline lengths
    .bazims     - baseline azimuths (in instrument coords)
    .u_coords   - u coordinates (in instrument coords)
    .v_coords   - v coordinates (in instrument coords)

It also has attributes containing geometry data, as follows:
    .BL2H_IX     - The 'baseline-to-hole indices'. A (numBLs x 2) array which specifies the pair of hole
                        numbers corresponding to each baseline.
                        I.e., given a baseline, bl2h_ix gives the 2 holes that go to make it up
    .H2BL_IX     - The 'hole-to-baseline indices'. An (numHoles x numHoles) array which baselines correspond
                        to any pair of holes.
                        I.e., given a pair of holes i,j H2BL_IX[i,j] gives the number of the baseline
    .BS2BL_IX    - The 'bispectrum-to-baseline' indices. Relates bispectrum/closure phase triangle indices
                        to baselines. I.e., given a point in the bispectrum, BS2BL_IX gives the 3 baselines
                        which make the triangle.
    .BL2BS_IX    - The 'baseline-to-bispectrum indexes'. BL2BS_IX gives the index of all points in the
                        bispectrum containing a given baseline
    .FILTER      - 2-element array giving the centre wavelength and bandpass of the filter

and metadata:
    .UTCs        - UTC times for each observation
    .ras, .decs  - telescope RA and dec for each observation
    .mask        - name of aperture mask used
    .emgains     - EM gain of camera for each observation
    .mffile      - the name of the matched-filter file used in data reduction. Thie file contains
                   information on u,v sampling, correspondence of closure phase triangles to baselines, etc.
    .pkflux      - peak flux in frame for each observation
    .totflux     - total flux in frame for each observation


Make a plot of the differential visibilities with
    diffData.plotVHVVdata()
or of the differential closure phases with
    diffdata.plotDiffCPdata()

"""

import numpy as np
from scipy import io
import matplotlib.pyplot as plt

class vampDiffdata:
    def __init__(self, filename, cubeInfoFilename):
        """
        Reads in the polarised-differential calibrated VAMPIRES data produced by the IDL pipeline
        :param filename: The full filename of the data file. E.g. diffdata_vega_.......idlvar
        """

        dObj = io.readsav(filename, python_dict=False, verbose=False)
        self.vhvv = dObj.vhvv
        self.vhvverr = dObj.vhvverr
        self.vhvvu = dObj.vhvvu
        self.vhvvuerr = dObj.vhvvuerr
        self.blengths = dObj.blengths
        self.bazims = dObj.bazims
        self.inFilename = filename
        try:
            self.diffCP = dObj.cp
            self.diffCPerr = dObj.cperr
            self.diffCPu = dObj.cpu
            self.diffCPuerr = dObj.cpuerr
            self.BL2H_IX = dObj.BL2H_IX
            self.H2BL_IX = dObj.H2BL_IX
            self.BL2BS_IX = dObj.BL2BS_IX
            self.BS2BL_IX = dObj.BS2BL_IX
        except:
            print "Couldn't find diff CP data for "+filename
        try:
            self.u_coords = dObj.u_coords
            self.v_coords = dObj.v_coords
        except:
            # Some older files didn't have these saved
            self.u_coords = ()
            self.v_coords = ()
        del (dObj)

        # Get useful metadata from cubeinfo file
        cubeinfoObj = io.readsav(cubeInfoFilename, python_dict=False, verbose=False)
        self.UTCs = cubeinfoObj.olog.utc[0]
        self.filters = cubeinfoObj.olog.filter[0]
        self.ras = cubeinfoObj.olog.ra[0]
        self.decs = cubeinfoObj.olog.dec[0]
        self.mask = cubeinfoObj.olog.mask[0]
        self.adate = cubeinfoObj.olog.adate[0]
        self.emgains = cubeinfoObj.olog.emgain[0]
        self.mffile = cubeinfoObj.plog.mf_file[0]
        self.pkflux = cubeinfoObj.framestats.pkflx[0]
        self.totflux = cubeinfoObj.framestats.totflx[0]
        self.cubename = cubeinfoObj.olog.cube_fname[0][0]
        del (cubeinfoObj)
        print "Read "+filename


    def plotVHVVdata(self, maxBL = 8., figNum=1):
        az = self.bazims
        bl = self.blengths
        blCols = bl #/ np.max(bl)

        az = az[bl <= maxBL]
        vhvv = self.vhvv[bl <= maxBL]
        vhvverr = self.vhvverr[bl <= maxBL]
        vhvvu = self.vhvvu[bl <= maxBL]
        vhvvuerr = self.vhvvuerr[bl <= maxBL]
        blCols = blCols[bl <= maxBL]
        bl = bl[bl <= maxBL]

        plt.figure(figNum)
        # This is a horrible hack to get the error bar colors to match...
        self.ax = plt.subplot(211)
        scatterPlt = self.ax.scatter(az, vhvv, c=blCols, marker='x', alpha=0.8)
        clb = plt.colorbar(scatterPlt)
        barColor = clb.to_rgba(blCols)

        plt.clf()
        #plt.suptitle(self.inFilename)
        self.ax = plt.subplot(211)
        scatterPlt = self.ax.scatter(az, vhvv, c=blCols, marker='x', alpha=0.8)
        a,b,c = self.ax.errorbar(az, vhvv, yerr=vhvverr, marker='', ls='',
                             alpha=0.8, capsize=0, zorder=0)
        c[0].set_color(barColor)
        sigString = "$\sigma$ = %.4f" % np.std(vhvv)
        chi2String = "$\chi^2_{null}$ = %.2f" % self.reducedChi2(vhvv, vhvverr)
        infoString = sigString + '      ' + chi2String
        self.ax.text(0.5, 0.9, infoString, horizontalalignment='center', verticalalignment='center',
                     transform = self.ax.transAxes)
        self.ax.text(0.5, 0.1, self.inFilename, horizontalalignment='center', verticalalignment='center',
                    transform=self.ax.transAxes, fontsize=8)
        plt.colorbar(scatterPlt, label='Baseline length (m)')
        plt.tight_layout()

        self.ax = plt.subplot(212)
        scatterPlt = self.ax.scatter(az, vhvvu, c=blCols, marker='x', alpha=0.8)
        a,b,c = self.ax.errorbar(az, vhvvu, yerr=vhvvuerr, marker='', ls='',
                             alpha=0.8, capsize=0, zorder=0)
        c[0].set_color(barColor)
        sigString = "$\sigma$ = %.4f" % np.std(vhvvu)
        chi2String = "$\chi^2_{null}$ = %.2f" % self.reducedChi2(vhvvu, vhvvuerr)
        infoString = sigString + '      ' + chi2String
        self.ax.text(0.5, 0.9, infoString, horizontalalignment='center', verticalalignment='center',
                     transform = self.ax.transAxes)
        self.ax.text(0.5, 0.1, self.inFilename, horizontalalignment='center', verticalalignment='center',
                    transform=self.ax.transAxes, fontsize=8)
        cb = plt.colorbar(scatterPlt, label='Baseline length (m)')
        plt.tight_layout()


    def reducedChi2(self, vhvv, vhvverr):
        # Reduced chi^2 for null hypothesis (vhvv=1)
        chi2 = np.sum( (vhvv-1.)**2 / vhvverr**2 )
        dof = len(vhvv) - 1
        return chi2 / dof


    def plotDiffCPdata(self, figNum=1, nBins=25):
        try:
            self.diffCP
        except:
            print "No closure phase data present for this file"
        else:
            cp = self.diffCP/np.pi * 180
            cpU = self.diffCPu/np.pi * 180
            plt.figure(figNum)
            plt.clf()
            #self.figureObj.suptitle(vhvvData[0])
            self.ax = plt.subplot(211)
            histPlt = self.ax.hist(cp, nBins)
            self.ax.set_xlabel('Differential closure phase (deg)')
            self.ax.set_ylabel('Frequency')
            sigString = "$\sigma$ = %.3f" % np.std(cp)
            self.ax.text(0.99, 0.9, sigString, horizontalalignment='right', verticalalignment='center',
                         transform=self.ax.transAxes)
            self.ax.text(0.99, 0.75, self.inFilename, horizontalalignment='right', verticalalignment='center',
                         transform=self.ax.transAxes, fontsize=8)

            self.ax = plt.subplot(212)
            histPlt = self.ax.hist(cpU, nBins)
            self.ax.set_xlabel('Differential closure phase (deg)')
            self.ax.set_ylabel('Frequency')
            sigString = "$\sigma$ = %.3f" % np.std(cpU)
            self.ax.text(0.99, 0.9, sigString, horizontalalignment='right', verticalalignment='center',
                         transform=self.ax.transAxes)
            self.ax.text(0.99, 0.75, self.inFilename, horizontalalignment='right', verticalalignment='center',
                         transform=self.ax.transAxes, fontsize=8)
            plt.tight_layout()



#######################################################################################################


srcPath = "/Users/bnorris/Don't Backup/etaCrv_Multi_18hN_redoMF/"
srcCubeInfoFilename = 'cubeinfoAug2016.idlvar'
srcFilename = 'diffdata_etaCrv_18hNCombined_20160320_750-50_18holeNudged_0_0.idlvar'


diffData = vampDiffdata(srcPath+srcFilename, srcPath+srcCubeInfoFilename)

diffData.plotVHVVdata()
diffData.plotDiffCPdata(figNum=2)
plt.show()
