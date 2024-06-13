from vmtk import pypes
import os
import numpy as np


#Getting Target and Source Points From vmtknetworkextraction For Centerline Calculations
def getSource(identifier: str, savePath: str):
    file = open(savePath + 'centerPre/' + identifier + '.dat')
    lst = []
    lstTotal = None
    
    for line in file:
        lst += [line.split()]
    if len(lst) > 2:
        if lst[2][3] == 0:
            lstTotal = len(lst)-2
            xTarg = lst[3][0]
            yTarg = lst[3][1]
            zTarg = lst[3][2]
        else:
            lstTotal = len(lst)-2
            xTarg = lst[2][0]
            yTarg = lst[2][1]
            zTarg = lst[2][2]

        xSource = lst[lstTotal][0]
        ySource = lst[lstTotal][1]
        zSource = lst[lstTotal][2]
    else:
        xTarg = '1'
        yTarg = '1'
        zTarg = '1'
        
        xSource = '1'
        ySource = '1'
        zSource = '1'

    targs = [xTarg,yTarg,zTarg]
    source = [xSource,ySource,zSource]

    return targs,source

    
#VMTK Scripts for Getting Aorta Geometry
def getDims(pathToCT: str, tempPath: str, identifier: str, savePath: str):
    
    pypes.PypeRun('vmtkimagereader -ifile ' + pathToCT + ' -origin 0.0 0.0 0.0 --pipe vmtkimagewriter -ofile ' + tempPath + 'aorta-label.vti')
    #For multilabel add -upperthreshold and make both that and -lowerthreshold 52.0
    pypes.PypeRun('vmtkimageinitialization -ifile ' + tempPath + 'aorta-label.vti -interactive 0 -method threshold -lowerthreshold 1.0 --pipe vmtklevelsetsegmentation -ifile ' + tempPath + 'aorta-label.vti -iterations 300 -ofile ' + tempPath + 'testing.vti')
    pypes.PypeRun('vmtkmarchingcubes -ifile ' + tempPath + 'testing.vti -ofile ' + tempPath + 'mcSurf.vtp --pipe')
    pypes.PypeRun('vmtksurfacesmoothing -ifile ' + tempPath + 'mcSurf.vtp -passband 0.005 -iterations 30 -ofile ' + savePath + 'smoothedModel/' + identifier + '.vtp --pipe')
    pypes.PypeRun('vmtknetworkextraction -ifile ' + savePath + 'smoothedModel/' + identifier + '.vtp -advancementratio 1.10 -ofile ' + tempPath + 'centerPre.vtp')
    pypes.PypeRun('vmtksurfacewriter -ifile ' + tempPath + 'centerPre.vtp -ofile ' + savePath + 'centerPre/' + identifier + '.dat')

    targs,sources = getSource(identifier, savePath)

    source = sources[0] + ' ' + sources[1] + ' ' + sources[2]
    target = targs[0] + ' ' + targs[1] + ' ' + targs[2]


    pypes.PypeRun('vmtksurfacereader -ifile ' + savePath + 'smoothedModel/' + identifier + '.vtp --pipe vmtkcenterlines -seedselector pointlist -sourcepoints ' + source + ' -targetpoints ' + target + ' -ofile ' + tempPath + 'center.vtp')
    pypes.PypeRun('vmtksurfacereader -ifile ' + savePath + 'smoothedModel/' + identifier + '.vtp --pipe vmtkdistancetocenterlines -centerlinesfile ' + tempPath + 'center.vtp -ofile ' + tempPath + 'centerDist.vtp')   
    pypes.PypeRun('vmtkcenterlinegeometry -ifile ' + tempPath + 'center.vtp -ofile ' + savePath + 'geomRaw/' + identifier + '.vtp')
    pypes.PypeRun('vmtksurfacewriter -ifile ' + savePath + 'geomRaw/' + identifier + '.vtp -ofile ' + savePath + 'geomData/' + identifier + '.dat')
    pypes.PypeRun('vmtksurfacewriter -ifile ' + tempPath + 'centerDist.vtp -ofile ' + savePath + 'centerDist/' + identifier + '.dat')
    pypes.PypeRun('vmtksurfacewriter -ifile ' + tempPath + 'center.vtp -ofile ' + savePath + 'centerLines/' + identifier + '.vtp')

################################################################

#File Path to Where Labels Are Stores
base = '/mnt/e/aorta/chest/without_contrast/'
os.chdir(base)

files = os.listdir()

#Temporary File Storage
tempFilePath = '/mnt/e/tempFiles/'

#Path to save outputs
savePath = '/mnt/e/CenterlinesAndRadiusTest/'


#Generating Necessary File Paths
if os.path.isdir(tempFilePath) == False:
    os.makedirs(tempFilePath)
    
if os.path.isdir(savePath) == False:
    os.makedirs(savePath)
    os.makedirs(savePath + 'centerDist')
    os.makedirs(savePath + 'centerLines')
    os.makedirs(savePath + 'centerPre')
    os.makedirs(savePath + 'geomData')
    os.makedirs(savePath + 'geomRaw')
    os.makedirs(savePath + 'smoothedModel')

#Running by iterating through all files
for i in range(0,3):
    os.chdir(base+files[i])
    tempBase = base + files[i] + '/' + os.listdir()[0] + '/aorta.seg.nrrd' 
    getDims(tempBase,tempFilePath,os.listdir()[0],savePath)
    print(os.listdir()[0])

