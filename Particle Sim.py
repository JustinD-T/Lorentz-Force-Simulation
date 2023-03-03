# Libraries
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import random
random.seed()
from alive_progress import alive_bar
testMode = {'quickStart' : True, 'multiTest' : False, 'ranAngle' : False}

# Changes:
    # Change meshgrid to work for differntiating theta values

mainDir = 'C:\Coding Projects1\Particle Sim\Testing Varuibles'

def programInitialize(): #Initializes varuibles from json

    # Prints directory contents out for selection
    if testMode['quickStart'] == True:
        fNum = 0
    else:
        for i in range(0, len(os.listdir(mainDir))): 
            print('['+str(i)+'] '+os.listdir(mainDir)[i])        
        fNum = input('Choose a File: ')

    # Initializes json into a dictionary
    directory = os.path.join(mainDir, os.listdir(mainDir)[int(fNum)])
    with open(directory, 'r') as j:
        varDict = json.load(j, parse_float = lambda x: float(x))
    
    # Particle array creation [[x values], [y values], [z values]]
    
    return varDict

def magLineVis(varDict): #generates a figure with magentic lines
    # PROBLEMS:
        # For unit vectors negitive/0 or larger than ~4, lengths breakdown

    def magLineParamFunc(i, j, k, meshParam, color='m'): # Parameterized function for graphing

    # Meshgrid returns 2d array of points with 'mdis' distance between range 'mlow'-'mlim
        x, y, z = np.meshgrid(np.arange(meshParam['mlow'], meshParam['mlim'], meshParam['mdis']), np.arange(meshParam['ylow'], meshParam['ylim'], meshParam['ydis']), np.arange(meshParam['mlow'], meshParam['mlim'], meshParam['mdis'])) 
    
    # quiver.(origin x, y, z : ending i, j, k - determines vector)
        ax.quiver(x, y, z, i, j, k, color=color, alpha=0.75)


# Initializes figure and axis, with 3d projection
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    gBounds = [-5, 5] #controls bounds of graph    
    ax.set_xlim(gBounds)
    ax.set_ylim(gBounds)
    ax.set_zlim(gBounds)

    # Parameters for brute-testing
    ranLow = 0
    ranLim = 2.5

    # Generates and Creates meshgrid
    meshParam = {'mlim' : 5, 'mlow' : -5, 'mdis' : 2.5, 'ylim' : 5, 'ylow' : -5, 'ydis' : 3.5}
    if testMode['multiTest'] == True:
        colors = ['m', 'b', 'r']
        ranLow = 0
        ranLim = 2.5
        for c in colors:
            magLineParamFunc(round(random.uniform(ranLow, ranLim), 2), round(random.uniform(ranLow, ranLim), 2), round(random.uniform(ranLow, ranLim), 2), meshParam, c)
    else:
        if testMode['ranAngle'] == True:
            i = round(random.uniform(ranLow, ranLim), 2)
            j = round(random.uniform(ranLow, ranLim), 2)
            k = round(random.uniform(ranLow, ranLim), 2)
            print(i, j, k)
        else:
            i = varDict['chargeField']['unitVectors']['i']
            j = varDict['chargeField']['unitVectors']['j']
            k = varDict['chargeField']['unitVectors']['k']
        magLineParamFunc(i, j, k, meshParam)

# Plots initial particle vector
    # Arrow magnitude and origin hard-coded to normal magline quiver func, make dependant on meshParam
    ax.quiver(-5, -5, -5, varDict['particleVar']['unitVectors']['i']+5, varDict['particleVar']['unitVectors']['j']+5, varDict['particleVar']['unitVectors']['k']+5)

    plt.show()


def particlePath(varDict):

    # Loads common varuibles
    # Hard coded in mass, change to interpret from json
    mass = varDict['particleVar']['mass'] * 10 ** -27
    time = varDict['obsTime'] / varDict['nVal']
    q = varDict['particleVar']['chargeElem'] / 6.24*10**-18
    Bmag = varDict['chargeField']['bTesla']

    def vectorDecomp(i, j, k, m): #returns a decomposed vector array from values [i, j, k, m]
        try:
            # Decomposition of each vector into x, y, z
            x = np.cos( np.arctan(j / i) ) * np.cos( np.arctan(k / i)) * m
            z = np.cos( np.arctan(j / k) ) * np.cos( np.arctan(i / k)) * m
            y = np.sqrt(m**2 - x**2 - z**2)
        except:
            print('invalid arc-tan on itter', i, 'arc-tan value of', j / i, i / k)

        return x, y, z

    def forceCalc(varDict, i, j, k, v): # Calculates dimensional forces given particle unit vectors and velocity
        


        # Decomposes both particle varuibles and B varuibles into consituent dimensions
        Bx, By, Bz = vectorDecomp(varDict['chargeField']['unitVectors']['i'], varDict['chargeField']['unitVectors']['j'], varDict['chargeField']['unitVectors']['k'], varDict['chargeField']['bTesla'])
        Vx, Vy, Vz = vectorDecomp(i, j, k, v)

        # Performs vector crossproduct using dimensional consituentcies
            # Uses cosine law and distance formula to find angle between particle and mag field
        vectorDist = np.sqrt( (Bx - Vx)**2 + (By - Vy)**2 + (Bz - Vz)**2 )
        try:
            # Fails as arccos value is over one radian
                # on 14113: 0.8230597367305824
                # on 14114: 1.0046366127008255
            theta = np.arccos( (v**2 + Bmag**2 - vectorDist**2) / (2 * v * Bmag) )
        except RuntimeWarning:
            print('invalid arc-cos input:',  (v**2 + Bmag**2 - vectorDist**2) / (2 * v * Bmag) , 'itternum:', i)
        Fx, Fy, Fz = ( np.cross([Vx, Vy, Vz], [Bx, By, Bz]) * q ) * np.sin(theta)

        Vfx = np.sqrt( ( Fx / mass * time )**2 + Vx**2)
        Vfy = np.sqrt( ( Fy / mass * time )**2 + Vy**2)
        Vfz = np.sqrt( ( Fz / mass * time )**2 + Vz**2)

        Vf = np.sqrt( Vfx**2 + Vfy**2 + Vfz**2)      

        return Vfx, Vfy, Vfz, Vf


    velocityArray = []
    with alive_bar(varDict['nVal']) as bar:
        for i in range(varDict['nVal']):
            if i == 0:
                velocityArray.append(forceCalc(varDict, varDict['particleVar']['unitVectors']['i'], varDict['particleVar']['unitVectors']['j'], varDict['particleVar']['unitVectors']['k'], varDict['particleVar']['iVelocity']))
                bar()
            else:
                velocityArray.append(forceCalc(varDict, *velocityArray[i-1]))
                bar()

    print('final velocity', velocityArray[-1])



varDict = programInitialize()
magLineVis(varDict)
particlePath(varDict)
