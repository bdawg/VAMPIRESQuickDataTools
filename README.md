# VAMPIRES Data and Tools Description
*Version 0.1, October 2016*

Contact: Barnaby Norris - <mailto:barnaby.norris@sydney.edu.au>

## Contents
* [Polarisation-differential calibration](#sec_PolzDiffCal)
	* [File description](#sec_PolzDiffCal_FileDesc)
	* [Data format description & python tools](#sec_PolzDiffCal_FormatDesPython)
	* [Position angle offsets](#sec_PA)
* [Conventional non-polarised calibration](#sec_NonpolzCal)
	* [File description](#sec_NonpolzCal_FileDesc)
	* [Data format description & python tools](#sec_NonpolzCal_FormatDesPython)

## Overview
VAMPIRES data products are still in the form of IDL save files, as its aperture-masking data is initially reduced using the IDL aperture-masking pipeline. IDL is becoming uncommon, so a set of Python tools to read and perform basic tasks with this data is being developed. The Python tools are however in a very early stage of development!!

VAMPIRES allows two distinct calibration modes, leading to two distinct data product formats:

* **Polarisation-differential calibration** - This is the main mode VAMPIRES is designed for. It calibrates visibilities and closure phases between orthogonal polarisation states using a three-tiered polarisation switching scheme. This gives very good calibration precision (it is quasi-simultaneous).
* **Conventional non-polarised calibration** - This is just the traditional interferometry approach, where the visibilities and closure phases of the science target are calibrated with those observed from a separate calibrator star. This relies on the PSF being constant between targets, and currently this works less well as the SCExAO PSF seems to vary with time.

A description of the calibration scheme can be found in [Norris, et al. (2015), MNRAS 447 3](http://adsabs.harvard.edu/abs/2015MNRAS.447.2894N)

<a name="sec_PolzDiffCal"></a>
##  Polarisation-differential calibration

<a name="sec_PolzDiffCal_FileDesc"></a>
### File description
The data products of the differential calibration are stored in IDL save files with the naming format

	diffdata_[target name]_[description string]_[calmode].idlvar

or, if it has been split into PA bins,

	diffdata_[target name]_[description string]_PABinNum[num]_[PA Range]_[calmode].idlvar
	
``[calmode]`` refers to the number and combination of tiers used in the differential calibration.	By default, data products with several different calibration modes are produced for diagnostic purposes. **Mode 0 refers to the full 3-tiered calibration, and you should always use this unless you're doing something special!** The calmode number is defined as follows:

* 0 - Triple calibration. **Use this unless you have some special reason to do otherwise!**
* 1 - Double calibration, just use beamsplitter and LCVR
* 2 - Double calibration, just use beamsplitter and half-wave plate
* 3 - Double calibration, just use LCVR and half-wave plate
* 4 - Single calibration, just use beamsplitter
* 5 - Single calibration, just use LCVR
* 6 - Single calibration, just use half-wave plate

Following the number (except in the case of 0) is a letter, which specifies which of the possible combinations of that calmode are used. E.g. ``2a`` would calibrate using the beamsplitter and LCVR, using only data from half-wave plate position 0, while ``2b`` would calibrate using the beamsplitter and LCVR, using only data from half-wave plate position 1. Usually only the ``a`` file is saved.

Since the differential visibility data cannot simply be derotated and added (since this would also rotate the axis of the polarisation analyser with respect to the target), calibrated data is also saved into bins of position-angle (PA), to allow simultaneous model fitting to all of them. In these cases a bin number and range of PA in degrees ``[PA Range]`` is included in the filename.

In addition to the ``diffdata*`` files, a large collection of metadata is saved in an IDL save file with the naming format
	
	cubeinfo[id string].idlvar
	
Both sets of files can be parsed using the Python tools described subsequently...


<a name="sec_PA"></a>	
### Position angle offsets
The VAMPIRES software queries the Gen2 telescope control system upon acquiring each datacube, and records the relevant keyword/value pairs associated with the telescope pointing and image rotator. The VAMPIRES data reduction software then uses these values to derive a sky position angle. The simplest method to determine this angle is by using the IMR.PAD keyword (saved as *PA* in the Python tools), which is calculated by Gen2 as a function of the parallactic angle and the image-rotator angle.

However there remains a constant offset due to the various out-of-plane optics in SCExAO. Then PA_Sky, the position angle with respect to sky north, measured in degrees east of north, is given by
	PA_Sky = PA_Obs + IMR.PAD + offset
	where PA_Obs is the position angle observed with respect to VAMPIRES’ detector-up direction (where pixel (0,0) is at the bottom-left).

As of Q1 2016 the offset is measured to be

	offset = 54 +/- 4 degrees

	
<a name="sec_PolzDiffCal_FormatDesPython"></a>
### Data format description and Python tools

The simplest way to access the differential data is using ``readDiffdata.py`` and the ``vampDiffdata`` class therein. 

A data object can be returned by doing

	diffData = vampDiffdata(dataFilename, cubeInfoFilename)
	
where 	``dataFilename`` is the ``diffdata_*.idlvar`` file. This will return an object with the following attributes:

~~~
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
~~~

It also contains geometry data...

~~~
.BL2H_IX     - The 'baseline-to-hole indices'. A (numBLs x 2) array which specifies the pair of 
				hole numbers corresponding to each baseline. I.e., given a baseline, bl2h_ix 
				gives the 2 holes that go to make it up
.H2BL_IX     - The 'hole-to-baseline indices'. An (numHoles x numHoles) array which baselines 
				correspond to any pair of holes. I.e., given a pair of holes i,j H2BL_IX[i,j] 
				gives the number of the baseline
.BS2BL_IX    - The 'bispectrum-to-baseline' indices. Relates bispectrum/closure phase triangle 
				indices to baselines. I.e., given a point in the bispectrum, BS2BL_IX gives the 
				3 baselines which make the triangle.
.BL2BS_IX    - The 'baseline-to-bispectrum indexes'. BL2BS_IX gives the index of all points in 
				the bispectrum containing a given baseline
.FILTER      - 2-element array giving the centre wavelength and bandpass of the filter

~~~

and some metadata:

~~~
.UTCs        - UTC times for each observation
.ras, .decs  - telescope RA and dec for each observation
.mask        - name of aperture mask used
.emgains     - EM gain of camera for each observation
.mffile      - the name of the matched-filter file used in data reduction. This file contains
              	information on u,v sampling, correspondence of closure phase triangles to 
               	baselines, etc.
.pkflux      - peak flux in frame for each observation
.totflux     - total flux in frame for each observation
~~~

The class also includes methods to make simple plots of differential visibilities - ``.plotVHVVdata()`` - and closure phases - ``.plotDiffCPdata()``. An example usage would be

~~~	python
diffData = vampDiffdata(srcFilename, srcCubeInfoFilename)
diffData.plotVHVVdata()
diffData.plotDiffCPdata(figNum=2)
plt.show()
~~~
	
Additionally, a GUI tool (*vampCrawler*) is being developed to allow quick browsing and examination of reduced data, as described later in these notes.


<a name="sec_NonpolzCal"></a>
## Conventional non-polarised calibration
While the polarisation-differential data has already been calibrated, this is not the case for conventional calibration since there is not necessarily one 'correct' way to do it. To this end, the *vampCalNP.py* Python tool is being developed. This will perform visibility and closure-phase calibration between specified data sets and sub-sets. It will also provide various plotting functions. It is still in development but useful basic features are available. 

Performing conventional calibration of VAMPIRES data is complicated by the fact that there are 16 polarisation-switching states (4 halfwave-plate positions x 2 beam-splitter channels x 2 LCVR states). Each of these stats has a slightly different PSF due to spatially-dependent instrumental polarisation (such as curved mirrors or inconsistent coatings on optical surfaces) and, in the case of the beamsplitter channel, non-common path errors. Therefore best performance should be found by **calibrating data of each of these states against the corresponding state of the calibrator data**.


<a name="sec_NonpolzCal_FileDesc"></a>
### File description
The raw visibilities and closure-phases produced by the data-reduction pipeline are stored in a set of IDL save files with the naming format

	bs_[target name]_[description string]_[serial number]_[chan]_[lcvr].idlvar
	
where ``[serial number]`` corresponds to the file number of the original data (usually starting at 0 and counting up), ``[chan]`` is the beamsplitter channel, labelled as '1' or '2' and ``[lcvr]`` is the LCVR state, labelled as 'A' or 'B'. Data must be in the standard VAMPIRES acquisition sequence, i.e. the first file is HWP=0, the second file is HWP=22.5, then 45, then 67.5, then back to 0. Each file must have ``_1_A``, ``_1_B``, ``_2_A`` and ``_2_B`` variants.

The ``cubeinfo`` file described above is also used.


<a name="sec_NonpolzCal_FormatDesPython"></a>
### Data format description and Python tools
``vampCalNP.py`` consists of a set of functions used to calibrate the data and plot it in various ways. The accompanying script ``doCalScript.py`` gives some how-to examples.

To read in a set of data, several parameters must be specified:

* *srcPath* - path to the data files
* *srcPrefix* - the part of the date filename before the serial number, e.g. ``bs_vega_18hNCombined_20160101_750-50_18holeNudged_``
* *srcStartNum* - first file number to use. = 0 if you want the whole set.
* *numSrc* - the number of data files to read in. This must be a multiple of 4 (i.e. includes equal numbers of all HWP positions)
* *srcExtn* - usually this is ``.idlvar``
* *rootDir* - the directory containing the *templates* directory if using the full IDL pipeline. Otherwise, this should be ``./`` and the *templates* directory placed in the working folder.
* *srcCubeInfoFilename* - the filename of the *cubeinfo* file, as described earlier.

A set of data can then be read into a data object as follows:

~~~python
import vampCalNP as vc
...
sciData = vc.readDataSet(srcPath, srcPrefix, srcStartNum, numSrc, srcExtn, rootDir, 
				srcCubeInfoFilename)
~~~

This returns an object containing the entire data set. It is a list object with one entry for each dataset (so 4x the number of input files due to 4 polarisation states per file)
Each dataset has the following attributes:

~~~
.u_coords, .v_coords    - arrays of uv coordinates for that baseline
.vis2                   - array of squared visibilities
.vis2Err                - array of 1-sigma visibilty errors (NB covariance is not yet taken 
							into account!)
.cp                     - array of closure phases
.cpErr                  - array of 1-sigma closure phase errors (NB covariance is not yet 
							taken into account!)
.pa                     - *relative* position-angle for this data (see Position Angle 
							Offsets section)
.mfFilename             - the name of the matched-filter file used in data reduction. 
							This file contains information on (u,v) sampling, correspondence 
							of closure phase triangles to baselines, etc.
.hwpMode                - Integer specifying the half-wave plate state for this data, as 
							per the following:
							0 = 0 degrees; 1 = 22.5 degrees; 2 = 45 degrees; 3 = 67.5 degrees
.chanMode               - Integer specifying the polarising beamsplitter channel for this 
							data (0 or 1)
.lcvrMode               - Integer specifying the LCVR state for this data (0 or 1)
~~~


It also contains an object called ```mfFileData```, which contains metadata regarding the physical geometry of the mask, used, for example, to work out which closure phases correspond to which baselines, etc. It has these attributes:

~~~
.mfFileData.u
.mfFileData.v           - The u,v coordinates for each baseline
.mfFileData.BL2H_IX     - The 'baseline-to-hole indices'. A (numBLs x 2) array which 
							specifies the pair of hole numbers corresponding to each baseline.
							I.e., given a baseline, BL2H_IX gives the 2 holes that go to 
							make it up
.mfFileData.H2BL_IX     - The 'hole-to-baseline indices'. An (numHoles x numHoles) array 
							which specifies baselines correspond to any pair of holes.
							I.e., given a pair of holes i,j H2BL_IX[i,j] gives the number 
							of the baseline
.mfFileData.BS2BL_IX    - The 'bispectrum-to-baseline' indices. Relates bispectrum / 
							closure phase triangle indices to baselines. I.e., given a point 
							in the bispectrum, BS2BL_IX gives the 3 baselines which make the 
							triangle.
.mfFileData.BL2BS_IX    - The 'baseline-to-bispectrum indexes'. BL2BS_IX gives the index 
							of all points in the bispectrum containing a given baseline
.mfFileData.FILTER      - 2-element array giving the centre wavelength and bandpass of 
							the filter
~~~

See the various examples in ```doCalScript.py`` for some details on the various tools available. But for a quick starting point, a simple calibration might go like this:

~~~python
import vampCalNP as vc

# Read in data
sciData = vc.readDataSet(srcPath, srcPrefix, srcStartNum, numSrc, srcExtn, rootDir, srcCubeInfoFilename)
calData = vc.readDataSet(calPath, calPrefix, calStartNum, numCal, calExtn, rootDir, calCubeInfoFilename)

# Do basic calibration:
# Combine all calibrator data but keep it in its 16 polarisation states
calDataCombinedSepStates = vc.combineDataSeparated(calData)

# Make set of all science data calibrated against the matching calibrator polz states
calibratedDataSepStatesAll = vc.calibrateSeparately(sciData, calDataCombinedSepStates)

# De-rotate the calibrated data as per its PA
calibratedDataSepStatesAllDerot = vc.derotate(calibratedDataSepStatesAll)

# Make a plot of the calibrated data (only do this if you have ~10 input files, or graph
# gets too slow and cluttered)
vc.plotData.plotByPolzstate(calibratedDataSepStatesAllDerot, xlims=[0, 1.2e7], 
	ylims=[0,2], alpha = 0.1, figNum=1)
	
# Or make a 2D plot of the visibilities from just the first file:
vc.plotData.plot2DSingleVis2(calibratedDataSepStatesAllDerot, interpolate=True, figNum=2)
	
# Plot histogram of closure phases of calibrated data,
# keeping 16 polarisation states separate
vc.plotData.plotCPByPolzstate(calibratedDataSepStatesAllDerot, alpha=0.1, figNum=3)
~~~

***vampCalNP does not yet have OIFITS support*** - this is the next thing to be implemented. For now, once data is read into Python objects it can be saved into a generic format or manipulated as desired.





	

