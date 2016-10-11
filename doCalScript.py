"""
An example script showing a basic non-polarised calibration if VAMPIRES data
using vampCalNP. Various example commmands are shown at the end.

Note: if returning to the console after execution, you may need to do plt.show() to
manipulate plot windows (at least on OSX)
"""


import vampCalNP as vc
import matplotlib.pyplot as plt

################## Input parameters ####################

# Location of 'templates' directory:
rootDir = '/Volumes/bnorris/code_svn/masking/'

# Path to science data, file prefix (before _1_A, etc.), and details:
#srcPath = '/Volumes/silo4/snert/VAMPIRESData_201603/Analysis/etaCrv_Multi_18hN_redoMF/'
srcPath = "/Users/bnorris/Don't Backup/etaCrv_Multi_18hN_redoMF/"
srcPrefix = 'bs_etaCrv_18hNCombined_20160320_750-50_18holeNudged_'
srcStartNum = 0
numSrc = 200
numSrc = 12 # A small number for testing
srcExtn = '.idlvar'
srcCubeInfoFilename = 'cubeinfoAug2016.idlvar'

# Path to calibrator data, file prefix (before _1_A, etc.), and details:
#calPath = '/Volumes/silo4/snert/VAMPIRESData_201603/Analysis/delCrv-CAL_Multi_18hN_redoMF/'
calPath = "/Users/bnorris/Don't Backup/delCrv-CAL_Multi_18hN_redoMF/"
calPrefix = 'bs_delCrv_18hNCombined_20160320_750-50_18holeNudged_'
calStartNum = 0
numCal = 96
numCal = 12 # A small number for testing
calExtn = '.idlvar'
calCubeInfoFilename = 'cubeinfoAug2016.idlvar'

makeSomeExamplePlots = False

###################################################################################################
# Example procedure:

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


# calibratedDataSepStatesAllDerot now contains the calibrated data set.
# It is a list object with one entry for each dataset (so 4x the number of input files due to 4 polz states per file)
# Each dataset has the following attributes:
#   .u_coords, .v_coords    - arrays of uv coordinates for that baseline
#   .vis2                   - array of squared visibilities
#   .vis2Err                - array of 1-sigma visibilty errors (NB covariance is not yet taken into account!)
#   .cp                     - array of closure phases
#   .cpErr                  - array of 1-sigma closure phase errors (NB covariance is not yet taken into account!)
#   .pa                     - *relative* position-angle for this data (i.e. true pa - instrument offset)
#   .mfFilename             - the name of the matched-filter file used in data reduction. Thie file contains
#                             information on u,v sampling, correspondence of closure phase triangles to baselines, etc.
#   .hwpMode                - Integer specifying the half-wave plate state for this data, as per the following:
#                             0 = 0 degrees; 1 = 22.5 degrees; 2 = 45 degrees; 3 = 67.5 degrees
#   .chanMode               - Integer specifying the polarising beamsplitter channel for this data (0 or 1)
#   .lcvrMode               - Integer specifying the lcvr state for this data (0 or 1)
#
# It also contains an object called mfFileData, which contains metadata regarding the physical geometry of the mask,
# used, for example, to work out which closure phases correspond to which baselines, etc. It has these attributes:
#   .mfFileData.u
#   .mfFileData.u           - The u,v coordinates for each baseline
#   .mfFileData.BL2H_IX     - The 'baseline-to-hole indices'. A (numBLs x 2) array which specifies the pair of hole
#                             numbers corresponding to each baseline.
#                             I.e., given a baseline, bl2h_ix gives the 2 holes that go to make it up
#   .mfFileData.H2BL_IX     - The 'hole-to-baseline indices'. An (numHoles x numHoles) array which baselines correspond
#                             to any pair of holes.
#                             I.e., given a pair of holes i,j H2BL_IX[i,j] gives the number of the baseline
#   .mfFileData.BS2BL_IX    - The 'bispectrum-to-baseline' indices. Relates bispectrum/closure phase triangle indices
#                             to baselines. I.e., given a point in the bispectrum, BS2BL_IX gives the 3 baselines
#                             which make the triangle.
#   .mfFileData.BL2BS_IX    - The 'baseline-to-bispectrum indexes'. BL2BS_IX gives the index of all points in the
#                             bispectrum containing a given baseline
#   .mfFileData.FILTER      - 2-element array giving the centre wavelength and bandpass of the filter


if makeSomeExamplePlots:
    # Let's make a plot of the calibrated data (only do this if you have ~10 input files, or graph
    # gets too cluttered)
    vc.plotData.plotByPolzstate(calibratedDataSepStatesAllDerot, xlims=[0, 1.2e7], ylims=[0,2],
                                alpha = 0.1, figNum=1)

    # Or make a 2D plot of the visibilities from the first file:
    vc.plotData.plot2DSingleVis2(calibratedDataSepStatesAll, interpolate=True, figNum=2)

    # Plot histogram of closure phases of calibrated data,
    # keeping 16 polarisation states separate
    vc.plotData.plotCPByPolzstate(calibratedDataSepStatesAllDerot, alpha=0.1, figNum=3)

    plt.show()

# See below for many more plotting options...


######################################################################################
# More examples:
#
# # All calibrator data combined to single dataset
# calDataCombined = vc.combineData(calData)
#
# # All calibrator data combined but kept in its 16 polarisation states
# calDataCombinedSepStates = vc.combineDataSeparated(calData)
#
# # Set of all science data calibrated against a single combined cal
# calibratedDataSingle = vc.calibrateSingle(sciData, calDataCombined)
#
# # Set of all science data calibrated against the matching cal polz modes
# calibratedDataSepStatesAll = vc.calibrateSeparately(sciData, calDataCombinedSepStates)
#
# # Flatten all calibrated data for plotting.
# # Just for testing - this doesn't de-rotate the data!!!!!!
# calibratedDataSepStates = vc.combineDataSeparated(calibratedDataSepStatesAll)
#
# # Flatten all calibrated data for plotting and don't keep the 16 separate polz states.
# # Just for testing - this doesn't de-rotate the data!!!!!!
# calibratedDataCombined = vc.combineData(calibratedDataSepStatesAll)
#
# # De-rotate calibrated data
# calibratedDataSepStatesAllDerot = vc.derotate(calibratedDataSepStatesAll)


########## Some plotting commands #########
#
# # Plot the combined calibrator data
# vc.plotData.plotSingleVis2(calDataCombined)
#
# # Plot the combined calibrator data which has the 16 polarisation states separate
# # One plot window per HWP position, each with 4 colours for the sub-states
# vc.plotData.plotByHWP(calDataCombinedSepStates, ylims=[0,0.5], colorBySet=False)
#
# # Plot the combined calibrator data which has the 16 polarisation states separate
# # Plot each polz state separately
# vc.plotData.plotByPolzstate(calDataCombinedSepStates, ylims=[0,0.5])
#
# # Plot the calibrated data, keeping 16 polarisation states separate
# vc.plotData.plotByPolzstate(calibratedDataSepStates, ylims=[0,2])
#
# # Plot the calibrated data, NOT keeping 16 polarisation states separate
# vc.plotData.plotSingleVis2(calibratedDataCombined)
#
# # Make a 2D plot showing the power spectrum (inerpolated or scatter plot)
# vc.plotData.plot2DSingleVis2(calibratedDataCombined, interpolate=True)
# vc.plotData.plot2DSingleVis2(calibratedDataSepStatesAll, interpolate=True)
# vc.plotData.plot2DSingleVis2(calibratedDataSepStatesAllDerot, interpolate=True, flatten=True)
#
# Look at uv coverage:
# vc.plotData.viewUVcoverage(calibratedDataSepStatesAllDerot, alpha=0.01)
#
# # Plot histogram of closure phases of calibrated data,
# # keeping 16 polarisation states separate
# vc.plotData.plotCPByPolzstate(calibratedDataSepStates)
#
# # Plot histogram of closure phases of calibrator star data,
# # keeping 16 polarisation states separate
# vc.plotData.plotCPByPolzstate(calDataCombinedSepStates)

