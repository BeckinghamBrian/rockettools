# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:19:20 2020

@author: M27091
"""

from Engine_Class  import Engine

# basic engine starting info
OF = 1.6
Pc = 400 # psi
thrust = 3372.1342 # lbs
burnTime = 40 # seconds
height_of_optimization = 1000 # m above MSL

C = Engine(thrust, Pc, OF, height_of_optimization, burnTime)

C.engine_name = 'Hephaestus'

C.Lstar = 0.9 # m
C.alpha = 40 # degrees
C.beta = 30 # degrees
C.newOx = [['N2O4(L)'   ,'N 2 O 4'      ,    '96.5'],
            ['SiO2'      ,'Si 1 O 2'     ,    '3.5']]
C.newFuel = [['C2H5OH' ,'C 2 H 6 O 1'  ,    '75',          'Ethanol'],
              ['H2O'    ,'H 2 O 1'      ,    '25',          'Water']]
C.newFuelBool = 1

# Propellant/heat transfer stuffs
C.coolantTempStart = 310 # K temp of coolant from tank (37 C)
# C.Kwall = 385 # W/m*K thermal conductivity of chamber wall
# C.Kwall = 167
C.Kwall=398
C.Kc = 0.171 # W/m*K thermal conductivity of fuel
C.coefThermEx = 13*10**-6
C.youngMod = 140
# C.coefThermEx = 16.7*10**-6 # fractional expansion per degree C
# C.youngMod = 121 # GPa
C.chamberWallThickness = 1 # mm chamber wall thickness
C.channelHeight = 1 # mm channel height
C.channelWallThickness = 1 # mm thickness of wall between channels
C.numChannels = 40 # number of cooling channels (important for Re calcs)
C.temp_step = 0.1

# Isp vs OF plot criteria
C.OF_low = 0.1
C.OF_high = 2
C.temp_step = 1

# C.doItAll(300,'B:/brian/Documents/Rockets/Hephaestus/Engine/V2 Plots - Geo')
# C.getCeaProperties()
# C.engineDimensions()
# C.propFlowRates()
C.printChamberInfo()
# C.engineContour()
# C.plotNozzleContour()
# C.plotEngineContour()
# C.savePlots(300)
# C.allPlots()
# C.saveNozGeo('B:/brian/Documents/Rockets/Long Shot/Engine/V2 Plots - Geo')
# C.saveChamberGeo('B:/brian/Documents/Rockets/Hephaestus/Code')
# C.heatTransfer()
# C.heatTransferPlots()
# C.plotIspVOF()
