
"""
vampCalNP - VAMPIRES Calibration No-Polarisation

This is a simple code to perform conventional (i.e. not polarised-differential) calibration of VAMPIRES data, for a
given source and cal pair. This is meant to be a bare-bones script, offering a simple Python alternative to the IDL
calibrate_v2_cp.pro. It does not do anything clever!

It takes as input the bs_*.idlvar files and the cubeinfo file from the IDL masking pipeline,
and the correpsonding cubeinfo*.idlvar file.

IMPORTANT: Data must be in the standard VAMPIRES sequence, i.e. the first file is HWP=0, the second file is HWP=22.5,
then 45, then 67.5, then back to 0. Each file must have _1_A, _1_B, _2_A and _2_B variants.

A lot of useful metadata is in the IDL cubeinfo file, which is not currently used. Look at the olog
and plog structures in that file; variable names are fairly self-explanatory.

If returning to the console after execution, you may need to do plt.show() to manipulate plot
windows (plt.ion() doesn't work, at least not on osx)
"""

# ################## Example Input parameters ####################
#
# # Location of 'templates' directory:
# rootDir = '/Volumes/bnorris/code_svn/masking/'
#
# # Path to science data, file prefix (before _1_A, etc.), and details:
# #srcPath = '/Volumes/silo4/snert/VAMPIRESData_201603/Analysis/etaCrv_Multi_18hN_redoMF/'
# srcPath = "/Users/bnorris/Don't Backup/etaCrv_Multi_18hN_redoMF/"
# srcPrefix = 'bs_etaCrv_18hNCombined_20160320_750-50_18holeNudged_'
# srcStartNum = 0
# numSrc = 200
# srcExtn = '.idlvar'
# srcCubeInfoFilename = 'cubeinfoAug2016.idlvar'
#
# # Path to calibrator data, file prefix (before _1_A, etc.), and details:
# #calPath = '/Volumes/silo4/snert/VAMPIRESData_201603/Analysis/delCrv-CAL_Multi_18hN_redoMF/'
# calPath = "/Users/bnorris/Don't Backup/delCrv-CAL_Multi_18hN_redoMF/"
# calPrefix = 'bs_delCrv_18hNCombined_20160320_750-50_18holeNudged_'
# calStartNum = 0
# numCal = 96
# calExtn = '.idlvar'
# calCubeInfoFilename = 'cubeinfoAug2016.idlvar'


import numpy as np
from scipy import io
from time import sleep
from copy import deepcopy
from copy import copy
import matplotlib as mpl
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D


plt.ion()


def readDataSet(srcPath, srcPrefix, srcStartNum, numSrc, srcExtn, rootDir, srcCubeInfoFilename):
    # Covariances are not yet handled. cov matrices do exist in the .bs file...
    chanSuffixes = ['_1', '_2']
    lcvrSuffixes = ['_A', '_B']
    data = []
    cubeInfo = readCubeInfo(srcPath+srcCubeInfoFilename)
    mfFileData = readMFFile(rootDir+'templates/'+cubeInfo.mffile)
    count = 0 # To index arrays in cubeinfo
    hwpMode = 0
    for curNum in range(srcStartNum, srcStartNum+numSrc):
        filenumString = '%05d' % curNum
        for chanSuffix in chanSuffixes:
            for lcvrSuffix in lcvrSuffixes:
                curFilename = srcPath + srcPrefix + filenumString + chanSuffix + lcvrSuffix + srcExtn
                curData = readBSFile(curFilename)
                curData.pa = cubeInfo.pa[count]
                curData.hwpMode = hwpMode
                curData.mfFileData = mfFileData

                if chanSuffix == '_1':
                    curData.chanMode = 0
                else:
                    curData.chanMode = 1

                if lcvrSuffix == '_A':
                    curData.lcvrMode = 0
                else:
                    curData.lcvrMode = 1

                data.append(curData)
                count += 1

        if hwpMode == 3:
            hwpMode = 0
        else:
            hwpMode +=1
    return data


class readMFFile:
    def __init__(self, fileString):
        print fileString
        curFileObj = io.readsav(fileString[0], python_dict=False, verbose=False)
        self.BL2H_IX = curFileObj.BL2H_IX
        self.H2BL_IX = curFileObj.H2BL_IX
        self.BL2BS_IX = curFileObj.BL2BS_IX
        self.BS2BL_IX = curFileObj.BS2BL_IX
        self.HOLE_DIAM = curFileObj.HOLE_DIAM
        self.RAD_PIXEL = curFileObj.RAD_PIXEL
        self.FILTER = curFileObj.FILTER
        self.u = curFileObj.u
        self.v = curFileObj.v
        del (curFileObj)


class readBSFile:
    def __init__(self, fileString, empty=False):
        if not empty:
            print 'Reading file '+fileString
            curFileObj = io.readsav(fileString, python_dict=False, verbose=False)
            self.u_coords = curFileObj.u
            self.v_coords = curFileObj.v
            self.vis2 = curFileObj.v2
            self.vis2Err = curFileObj.v2_sig
            self.bispect = curFileObj.bs
            self.cp = curFileObj.cp
            self.cpErr = curFileObj.cp_sig
            self.mfFilename = curFileObj.mf_file1
            del (curFileObj)
        else:
            # Return an empty instance
            pass


class readCubeInfo:
    def __init__(self, cubeInfoFileString):
        # Get useful metadata from cubeinfo file
        cubeinfoObj = io.readsav(cubeInfoFileString, python_dict=False, verbose=False)
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
        self.pa = cubeinfoObj.olog.pa[0]
        del (cubeinfoObj)


def combineData(inData, chosenModes=[], llim=-1, ulim=-1):
    # Make weighted arithmetic mean of V2s and CPs
    # These arrays are of shape (nFiles, nV2s) (or CPs)
    # Note - SEM estimate does not yet take into account covariance!
    print 'Combining data with chosenMode:'
    print chosenModes

    curData = readBSFile('', empty=True)
    curData.u_coords = inData[0].u_coords
    curData.v_coords = inData[0].v_coords
    curData.mfFilename = inData[0].mfFilename
    curData.mfFileData = inData[0].mfFileData

    if len(chosenModes) == 0:
        chosenInds = range(len(inData))
    else:
        allHwpModes = np.asarray([s.hwpMode for s in inData])
        allChanModes = np.asarray([s.chanMode for s in inData])
        allLcvrModes = np.asarray([s.lcvrMode for s in inData])
        chosenHwpModes = chosenModes[0]
        chosenChanModes = chosenModes[1]
        chosenLcvrModes = chosenModes[2]
        chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                              & np.in1d(allLcvrModes, chosenLcvrModes))[0]

    chosenInds = np.asarray(chosenInds)
    if llim >= 0:
        chosenInds = chosenInds[np.where(chosenInds >= llim)]
    if ulim >= 0:
        chosenInds = chosenInds[np.where(chosenInds < ulim)]

    inDataChosen = [inData[i] for i in chosenInds]
    rawVis2 = np.asarray([s.vis2 for s in inDataChosen])
    rawVis2Err = np.asarray([s.vis2Err for s in inDataChosen])
    rawcp = np.asarray([s.cp for s in inDataChosen])
    rawcpErr = np.asarray([s.cpErr for s in inDataChosen])
    nFiles = len(inDataChosen)

    # If combining multiple modes, set these values to -1...
    if (sum(len(x) for x in chosenModes) > 3) or (len(chosenModes) == 0):
        curData.hwpMode = -1
        curData.chanMode = -1
        curData.lcvrMode = -1
    else:
        curData.hwpMode = chosenModes[0][0]
        curData.chanMode = chosenModes[1][0]
        curData.lcvrMode = chosenModes[2][0]

    vis2Weights = 1./(rawVis2Err**2)
    for i in range(nFiles):
        vis2Weights[i,:] = vis2Weights[i,:] / np.sum(vis2Weights[i,:])
    curData.vis2 = np.average(rawVis2, axis=0, weights=vis2Weights)
    curData.vis2NotWeighted = np.average(rawVis2, axis=0)
    # combinedVis2Err = np.zeros_like(combinedVis2)
    # for i in range(rawVis2.shape[1]):
    #     combinedVis2Err[i] = np.sqrt(np.average( (rawVis2[:,i]-combinedVis2[i])**2, weights=vis2Weights[:,i]))
    curData.vis2Err = np.sqrt(np.average((rawVis2 - curData.vis2) ** 2, axis=0, weights=vis2Weights)) \
                      / np.sqrt(nFiles)

    cpWeights = 1. / (rawcpErr ** 2)
    for i in range(nFiles):
        cpWeights[i, :] = cpWeights[i, :] / np.sum(cpWeights[i, :])
        curData.cp = np.average(rawcp, axis=0, weights=cpWeights)
        curData.cpErr = np.sqrt(np.average((rawcp - curData.cp) ** 2, axis=0, weights=cpWeights)) / np.sqrt(nFiles)

    return [curData]


def combineDataSeparated(inData):
    data = []
    for h in range(4):
        for c in range(2):
            for l in range(2):
                curChosenMode = [[h], [c], [l]]
                data.append(combineData(inData, curChosenMode)[0])
    return data


class plotData:
    def __init__(self):
        pass

    @staticmethod
    def plotSingleVis2(data, chosenModes = [], figNum = 1, subplot = 111, xlims=[0, 1.2e7], ylims=[0, 2],
                       clearPlot = True, fmt='bx', alpha=1, figsize = (6,4), flatten=False, errBars=True,
                       xlabel='Baseline length (sfu)', ylabel='V^2', title=''):
        if len(chosenModes) == 0:
            chosenInds = range(len(data))
        else:
            allHwpModes = np.asarray([s.hwpMode for s in data])
            allChanModes = np.asarray([s.chanMode for s in data])
            allLcvrModes = np.asarray([s.lcvrMode for s in data])
            chosenHwpModes = chosenModes[0]
            chosenChanModes = chosenModes[1]
            chosenLcvrModes = chosenModes[2]
            chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                                  & np.in1d(allLcvrModes, chosenLcvrModes))[0]

        blengths = np.sqrt(data[0].u_coords ** 2 + data[0].v_coords ** 2)
        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()
        plt.subplot(subplot)
        if flatten:
            chosenData = [data[i] for i in chosenInds]
            allU_coords = np.asarray([d.u_coords for d in chosenData]).flatten()
            allV_coords = np.asarray([d.v_coords for d in chosenData]).flatten()
            allBlengths = np.sqrt(allU_coords ** 2 + allV_coords ** 2)
            allVis2 = np.asarray([d.vis2 for d in chosenData]).flatten()
            allVis2Err = np.asarray([d.vis2Err for d in chosenData]).flatten()
            if errBars:
                plt.errorbar(allBlengths, allVis2, yerr=allVis2Err, fmt=fmt, alpha=alpha)
            else:
                plt.scatter(allBlengths, allVis2, marker='.', edgecolors='none', alpha=alpha)
        else:
            for k in chosenInds:
                if errBars:
                    plt.errorbar(blengths, data[k].vis2, yerr=data[k].vis2Err, fmt=fmt, alpha=alpha)
                else:
                    plt.scatter(blengths, data[k].vis2, marker='.', edgecolors='none', alpha=alpha)
                print 'HWP mode: %d, Chan mode: %d, LCVR mode: %d' % (
                        data[k].hwpMode, data[k].chanMode, data[k].lcvrMode)
                #plt.pause(0.001)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        #plt.plot(blengths, data.vis2NotWeighted, 'ro')
        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.pause(0.001)

    @staticmethod
    def plot2DSingleVis2(data, chosenModes=[], figNum=1, subplot=111, lim=1.2e7, cmap='viridis',
                       clearPlot=True, figsize=(6, 6), interpolate=True, step=1e5, alpha=0.5, sz=100,
                       xlabel='Baseline length (sfu)', ylabel='Baseline length (sfu)', title='', flatten=False,
                       vmin=-1, vmax=-1):
        if len(chosenModes) == 0:
            chosenInds = range(len(data))
        else:
            allHwpModes = np.asarray([s.hwpMode for s in data])
            allChanModes = np.asarray([s.chanMode for s in data])
            allLcvrModes = np.asarray([s.lcvrMode for s in data])
            chosenHwpModes = chosenModes[0]
            chosenChanModes = chosenModes[1]
            chosenLcvrModes = chosenModes[2]
            chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                                  & np.in1d(allLcvrModes, chosenLcvrModes))[0]

        if len(chosenInds) > 1 and flatten == False:
            print "WARNING: More than one dataset - only showing the first one. Try flatten=True"

        #blengths = np.sqrt(data[0].u_coords ** 2 + data[0].v_coords ** 2)
        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()
        plt.subplot(subplot)

        if flatten:
            chosenData = [data[i] for i in chosenInds]
            u_coords = np.asarray([d.u_coords for d in chosenData]).flatten()
            v_coords = np.asarray([d.v_coords for d in chosenData]).flatten()
            blengths = np.sqrt(u_coords ** 2 + v_coords ** 2)
            vis2 = np.asarray([d.vis2 for d in chosenData]).flatten()
            vis2Err = np.asarray([d.vis2Err for d in chosenData]).flatten()
        else:
            u_coords = data[chosenInds[0]].u_coords
            v_coords = data[chosenInds[0]].v_coords
            blengths = np.sqrt(u_coords ** 2 + v_coords ** 2)
            vis2 = data[chosenInds[0]].vis2

        if interpolate:
            grid_x, grid_y = np.mgrid[-lim:lim:step, -lim:lim:step]
            uAll = np.concatenate((u_coords, -u_coords))
            vAll = np.concatenate((v_coords, -v_coords))
            vis2All = np.concatenate((vis2, vis2))
            points = (uAll, vAll)
            griddedVis2 = griddata(points, vis2All, (grid_x, grid_y), method='cubic')
            if vmin > -1:
                griddedVis2[np.where(griddedVis2 <= vmin)] = vmin
            if vmax > -1:
                griddedVis2[np.where(griddedVis2 >= vmax)] = vmax
            plt.imshow(griddedVis2)
            plt.axis('equal')
        else:
            clippedVis2 = copy(vis2)
            if vmin > -1:
                clippedVis2[np.where(clippedVis2 <= vmin)] = vmin
            if vmax > -1:
                clippedVis2[np.where(clippedVis2 >= vmax)] = vmax
            cols = vis2 / np.max(vis2)
            uAll = np.concatenate((u_coords, -u_coords))
            vAll = np.concatenate((v_coords, -v_coords))
            colsAll = np.concatenate((cols, cols))
            plt.scatter(uAll, vAll, c=colsAll, marker='o', s=sz, edgecolors='none', alpha=alpha,
                        cmap=cmap)
            plt.axis('equal')
        plt.colorbar()


    @staticmethod
    def plotSurfSingleVis2(data, chosenModes=[], figNum=1, subplot=111, lim=1.2e7, cmap='viridis',
                       clearPlot=True, figsize=(6, 6), step=1e5, alpha=0.5, sz=100,
                       xlabel='Baseline length (sfu)', ylabel='Baseline length (sfu)', title='',
                        flatten=False, vmin=-1, vmax=-1):
        if len(chosenModes) == 0:
            chosenInds = range(len(data))
        else:
            allHwpModes = np.asarray([s.hwpMode for s in data])
            allChanModes = np.asarray([s.chanMode for s in data])
            allLcvrModes = np.asarray([s.lcvrMode for s in data])
            chosenHwpModes = chosenModes[0]
            chosenChanModes = chosenModes[1]
            chosenLcvrModes = chosenModes[2]
            chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                                  & np.in1d(allLcvrModes, chosenLcvrModes))[0]

        if len(chosenInds) > 1 and flatten == False:
            print "WARNING: More than one dataset - only showing the first one. Try flatten=True"

        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()
        ax=plt.subplot(subplot, projection='3d')

        if flatten:
            chosenData = [data[i] for i in chosenInds]
            u_coords = np.asarray([d.u_coords for d in chosenData]).flatten()
            v_coords = np.asarray([d.v_coords for d in chosenData]).flatten()
            blengths = np.sqrt(u_coords ** 2 + v_coords ** 2)
            vis2 = np.asarray([d.vis2 for d in chosenData]).flatten()
            vis2Err = np.asarray([d.vis2Err for d in chosenData]).flatten()
        else:
            u_coords = data[chosenInds[0]].u_coords
            v_coords = data[chosenInds[0]].v_coords
            blengths = np.sqrt(u_coords ** 2 + v_coords ** 2)
            vis2 = data[chosenInds[0]].vis2

        grid_x, grid_y = np.mgrid[-lim:lim:step, -lim:lim:step]
        uAll = np.concatenate((u_coords, -u_coords))
        vAll = np.concatenate((v_coords, -v_coords))
        vis2All = np.concatenate((vis2, vis2))
        points = (uAll, vAll)
        griddedVis2 = griddata(points, vis2All, (grid_x, grid_y), method='cubic')
        if vmin > 0:
            griddedVis2[np.where(griddedVis2 <= vmin)] = vmin
        if vmax > 0:
            griddedVis2[np.where(griddedVis2 >= vmax)] = vmax
        ax.plot_surface(grid_x, grid_y, griddedVis2)
        #ax.set_zlim(vmin, vmax)




    @staticmethod
    def plotByHWP(data, chosenModes = [[],[0,1],[0,1]], xlims=[0, 1.2e7], ylims=[0, 2], reClear = False,
                  colorBySet = False, waitTime = 0.1, errBars=True):
        # Note chosenModes will ignore HWP specified
        colors = ['b', 'g', 'r', 'c', 'm', 'k']
        colorCount = -1

        for curData in data:
            if curData.chanMode in chosenModes[1] and curData.lcvrMode in chosenModes[2]:
                if reClear:
                    if curData.hwpMode + curData.chanMode+ curData.lcvrMode == 0:
                        plt.clf()
                subplotNum = 221 + curData.hwpMode

                if not colorBySet:
                    if curData.chanMode == 0:
                        if curData.lcvrMode == 0:
                            f1 = 'r'
                        else:
                            f1 = 'g'
                    else:
                        if curData.lcvrMode == 0:
                            f1 = 'b'
                        else:
                            f1 = 'm'
                    f2 = 'x'
                else:
                    if curData.hwpMode + curData.chanMode + curData.lcvrMode == 0:
                        colorCount += 1
                        if colorCount == 6:
                            colorCount = 0
                    f1 = colors[colorCount]
                    f2 = 'x'

                titString = 'HWP position %d' % curData.hwpMode
                plotData.plotSingleVis2([curData], subplot = subplotNum, xlims = xlims, ylims = ylims,
                                clearPlot = False, fmt=f1+f2, alpha=0.5, figsize=(12, 8), title=titString,
                                errBars=errBars)
                sleep(waitTime)
        plt.tight_layout()


    @staticmethod
    def plotByPolzstate(data, xlims=[0, 1.2e7], ylims=[0, 2], alpha=1,
                  figNum = 1, clearPlot = True, fmt='bx', figsize = (15,10), errBars=True):
        # Makes a separate plot for each polz state (so 16 plots)
        count = 0
        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()

        for curData in data:
            subplotVal = curData.hwpMode*4 + curData.chanMode*2 + curData.lcvrMode
            plt.subplot(4, 4, subplotVal+1)
            # plt.tight_layout()
            # try:
            #     for subData in curData:
            #         blengths = np.sqrt(data[0][0].u_coords ** 2 + data[0][0].v_coords ** 2)
            #         plt.errorbar(blengths, subData.vis2, yerr=subData.vis2Err, fmt=fmt)
            #         titString = 'HWP: %d, Chan: %d, LCVR: %d' % (curData[0].hwpMode, curData[0].chanMode,
            #                                                      curData[0].lcvrMode)
            #         # print 'HWP mode: %d, Chan mode: %d, LCVR mode: %d' % (
            #         #     subData[k].hwpMode, subData[k].chanMode, subData[k].lcvrMode)
            # except:
            #     # Handle case where data is not a list
            blengths = np.sqrt(data[0].u_coords ** 2 + data[0].v_coords ** 2)
            if errBars:
                plt.errorbar(blengths, curData.vis2, yerr=curData.vis2Err, fmt=fmt, alpha=alpha)
            else:
                plt.scatter(blengths, curData.vis2, marker='.', s=10, edgecolors='none', alpha=alpha)
            titString = 'HWP: %d, Chan: %d, LCVR: %d' % (curData.hwpMode, curData.chanMode,
                                                         curData.lcvrMode)
            plt.xlim(xlims)
            plt.ylim(ylims)
            plt.title(titString)
            plt.xlabel('Baseline length (sfu)')
            plt.ylabel('V^2')
            plt.pause(0.001)
            count += 1
        plt.tight_layout()


    @staticmethod
    def plotSingleCP(data, chosenModes = [], figNum = 1, subplot = 111, xlims=(), ylims=(),
                       clearPlot = True, alpha=1, figsize = (6,4), bins=25, color='b',
                       xlabel='Closure phase (deg)', ylabel='Frequency', title=''):
        if len(chosenModes) == 0:
            chosenInds = range(len(data))
        else:
            allHwpModes = np.asarray([s.hwpMode for s in data])
            allChanModes = np.asarray([s.chanMode for s in data])
            allLcvrModes = np.asarray([s.lcvrMode for s in data])
            chosenHwpModes = chosenModes[0]
            chosenChanModes = chosenModes[1]
            chosenLcvrModes = chosenModes[2]
            chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                                  & np.in1d(allLcvrModes, chosenLcvrModes))[0]

        blengths = np.sqrt(data[0].u_coords ** 2 + data[0].v_coords ** 2)
        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()
        plt.subplot(subplot)
        for k in chosenInds:
            plt.hist(data[k].cp/np.pi*180, bins=bins, alpha=alpha, color=color)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if len(xlims) != 0:
                plt.xlim(xlims)
            if len(ylims) != 0:
                plt.ylim(ylims)
            curTit = title + 'sd = %.2f deg' % np.std(data[k].cp/np.pi*180)
            plt.title(curTit)
            plt.tight_layout()
            print 'HWP mode: %d, Chan mode: %d, LCVR mode: %d' % (data[k].hwpMode, data[k].chanMode, data[k].lcvrMode)
        plt.pause(0.001)


    @staticmethod
    def plotCPByHWP(data, chosenModes = [[],[0,1],[0,1]], reClear = False,
                  colorBySet = False, waitTime = 0.1):
        # Note chosenModes will ignore HWP specified
        colors = ['b', 'g', 'r', 'c', 'm', 'k']
        colorCount = -1

        for curData in data:
            if curData.chanMode in chosenModes[1] and curData.lcvrMode in chosenModes[2]:
                if reClear:
                    if curData.hwpMode + curData.chanMode+ curData.lcvrMode == 0:
                        plt.clf()
                subplotNum = 221 + curData.hwpMode

                if not colorBySet:
                    if curData.chanMode == 0:
                        if curData.lcvrMode == 0:
                            f1 = 'r'
                        else:
                            f1 = 'g'
                    else:
                        if curData.lcvrMode == 0:
                            f1 = 'b'
                        else:
                            f1 = 'm'
                else:
                    if curData.hwpMode + curData.chanMode + curData.lcvrMode == 0:
                        colorCount += 1
                        if colorCount == 6:
                            colorCount = 0
                    f1 = colors[colorCount]

                titString = 'HWP position %d --- ' % curData.hwpMode
                plotData.plotSingleCP([curData], subplot = subplotNum,
                               clearPlot = False, color=f1, alpha=0.2, figsize=(12, 8), title=titString)
                sleep(waitTime)
        plt.tight_layout()


    @staticmethod
    def plotCPByPolzstate(data, xlims=(), ylims=(), color='b', alpha=1,
                  bins=25, figNum = 1, clearPlot = True, figsize = (15,10)):
        # Makes a separate plot for each polz state (so 16 plots)
        count = 0
        plt.figure(figNum, figsize=figsize)
        if clearPlot:
            plt.clf()
        for curData in data:
            subplotVal = curData.hwpMode*4 + curData.chanMode*2 + curData.lcvrMode
            plt.subplot(4, 4, subplotVal+1)
            plt.hist(curData.cp / np.pi * 180, bins=bins, alpha=alpha, color=color)
            titString = 'HWP: %d, Chan: %d, LCVR: %d - sd = %.2f deg' % (curData.hwpMode,
                curData.chanMode, curData.lcvrMode, np.std(curData.cp/np.pi*180))
            if len(xlims) != 0:
                plt.xlim(xlims)
            if len(ylims) != 0:
                plt.ylim(ylims)
            plt.title(titString, fontsize=10)
            plt.xlabel('Closure phase (deg)')
            plt.ylabel('Frequency')
            plt.pause(0.001)
            count += 1
        plt.tight_layout()


    @staticmethod
    def viewUVcoverage(data, alpha=0.5):
        u = np.asarray([d.u_coords for d in data])
        v = np.asarray([d.v_coords for d in data])
        u = np.concatenate((u, -1*u))
        v = np.concatenate((v, -1*v))
        plt.scatter(u, v, marker='.', alpha=alpha)


def calibrateSingle(sciData, calData, chosenModes=[], llim=-1, ulim=-1):
    # sciData can be a set of data, but calData should be a single (i.e. combined) set
    # the polz modes of returned data will be that of the source data
    # NOTE - doesn't yet handle covariance!
    if len(calData) > 1:
        print "WARNING: calData has more than one entry - only using the first one"
    print 'Combining data with chosenMode:'
    print chosenModes

    data = []
    if len(chosenModes) == 0:
        chosenInds = range(len(sciData))
    else:
        allHwpModes = np.asarray([s.hwpMode for s in sciData])
        allChanModes = np.asarray([s.chanMode for s in sciData])
        allLcvrModes = np.asarray([s.lcvrMode for s in sciData])
        chosenHwpModes = chosenModes[0]
        chosenChanModes = chosenModes[1]
        chosenLcvrModes = chosenModes[2]
        chosenInds = np.where(np.in1d(allHwpModes, chosenHwpModes) & np.in1d(allChanModes, chosenChanModes)
                              & np.in1d(allLcvrModes, chosenLcvrModes))[0]
    chosenInds = np.asarray(chosenInds)
    if llim >= 0:
        chosenInds = chosenInds[np.where(chosenInds >= llim)]
    if ulim >= 0:
        chosenInds = chosenInds[np.where(chosenInds < ulim)]
    sciDataChosen = [sciData[i] for i in chosenInds]

    for curData in sciDataChosen:
        calibdData = readBSFile('', empty=True)
        calibdData.u_coords = curData.u_coords
        calibdData.v_coords = curData.v_coords
        calibdData.mfFilename = curData.mfFilename
        calibdData.mfFileData = curData.mfFileData
        calibdData.vis2 = curData.vis2 / calData[0].vis2
        calibdData.vis2Err = np.abs(calibdData.vis2) * np.sqrt( (curData.vis2Err / curData.vis2)**2 +
                                                                (calData[0].vis2Err / calData[0].vis2)**2 )
        calibdData.cp = curData.cp - calData[0].cp
        calibdData.cpErr = np.sqrt( curData.cpErr**2 + calData[0].cpErr**2 )
        calibdData.hwpMode = curData.hwpMode
        calibdData.chanMode = curData.chanMode
        calibdData.lcvrMode = curData.lcvrMode
        calibdData.pa = curData.pa
        data.append(calibdData)

    return data


def calibrateSeparately(sciData, calData, llim=-1, ulim=-1):
    # This will calibrate like-with-like. calData should be the list returned from combineRawDataSeparated,
    # with each polarisation mode only included once.
    # This version puts all data into one long list
    if not len(calData) == 16:
        print "WARNING: calData doesn't seem to have the right number of entries..."
    data = []
    for curCalData in calData:
        curChosenMode = [[curCalData.hwpMode], [curCalData.chanMode], [curCalData.lcvrMode]]
        curCalibdData = calibrateSingle(sciData, [curCalData], curChosenMode, llim=llim, ulim=ulim)
        for c in curCalibdData:
            data.append(c)
    return data


def derotate(data):
    derotatedData=[]
    data2 = deepcopy(data)
    for curData in data2:
        ang = -curData.pa /180*np.pi
        rmat = np.array([ [np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)] ])
        newCoords = np.dot(rmat, np.asarray([curData.u_coords, curData.v_coords]))
        curData.u_coords = newCoords[0,:]
        curData.v_coords = newCoords[1,:]
        derotatedData.append(curData)
    return derotatedData





###################################################################################################
# # Example procedure:
#
# # Read in data
# sciData = vc.readDataSet(srcPath, srcPrefix, srcStartNum, numSrc, srcExtn, rootDir, srcCubeInfoFilename)
# calData = vc.readDataSet(calPath, calPrefix, calStartNum, numCal, calExtn, rootDir, calCubeInfoFilename)
#
# # Do basic calibration:
# # Combine all calibrator data but keep it in its 16 polarisation states
# calDataCombinedSepStates = vc.combineDataSeparated(calData)
#
# # Make set of all science data calibrated against the matching calibrator polz states
# calibratedDataSepStatesAll = vc.calibrateSeparately(sciData, calDataCombinedSepStates)
#
# # De-rotate the calibrated data as per its PA
# calibratedDataSepStatesAllDerot = vc.derotate(calibratedDataSepStatesAll)
#
#
# # calibratedDataSepStatesAllDerot now contains the calibrated data set.
# # It is a list object with one entry for each dataset (so 4x the number of input files due to 4 polz states per file)
# # Each dataset has the following attributes:
# #   .u_coords, .v_coords    - arrays of uv coordinates for that baseline
# #   .vis2                   - array of squared visibilities
# #   .vis2Err                - array of 1-sigma visibilty errors (NB covariance is not yet taken into account!)
# #   .cp                     - array of closure phases
# #   .cpErr                  - array of 1-sigma closure phase errors (NB covariance is not yet taken into account!)
# #   .pa                     - *relative* position-angle for this data (i.e. true pa - instrument offset)
# #   .mfFilename             - the name of the matched-filter file used in data reduction. Thie file contains
# #                             information on u,v sampling, correspondence of closure phase triangles to baselines, etc.
# #   .hwpMode                - Integer specifying the half-wave plate state for this data, as per the following:
# #                             0 = 0 degrees; 1 = 22.5 degrees; 2 = 45 degrees; 3 = 67.5 degrees
# #   .chanMode               - Integer specifying the polarising beamsplitter channel for this data (0 or 1)
# #   .lcvrMode               - Integer specifying the lcvr state for this data (0 or 1)
#
#
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
#
#
# # Let's make a plot of the calibrated data (only do this if you have ~10 input files, or graph
# # gets too cluttered)
# vc.plotData.plotByPolzstate(calibratedDataSepStatesAllDerot, xlims=[0, 1.2e7], ylims=[0,2],
#                             alpha = 0.1, figNum=1)
#
# # Or make a 2D plot of the visibilities from the first file:
# vc.plotData.plot2DSingleVis2(calibratedDataSepStatesAll, interpolate=True, figNum=2)
#
# # Plot histogram of closure phases of calibrated data,
# # keeping 16 polarisation states separate
# vc.plotData.plotCPByPolzstate(calibratedDataSepStatesAllDerot, alpha=0.1, figNum=3)
#
# plt.show()
#
# # See below for many more plotting options...
#
#
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

