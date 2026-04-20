# -*- coding: utf-8 -*-
"""
Created on Sat Jan  8 13:41:48 2022

@author: brian
"""

from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import math
import MOC_nozzle as MOC
import pandas as pd
from os import path
import numpy as np
from molmass import Formula

class Engine():

    def __init__(self, thrust=0, Pc=0, OF=0, height_of_optimization=0, burnTime=0):
        
        self.engine_name = ''
        self.thrust = thrust * 4.44822162 # convert to N
        self.burnTime = burnTime
        self.Pc_psi = Pc 
        self.Pc = self.Pc_psi * 6894.75729 # convert to Pa
        self.OF = OF
        self.num_characteristics = 300
        self.Lstar = 1
        self.alpha = 30
        self.beta = 15
        self.height_of_optimization = height_of_optimization
        self.Pa = 101325*(1-(2.25577)*(10**(-5))*self.height_of_optimization)**5.25588
        self.Pe = self.Pa # assumed perfect expansion
        self.PcOvPe = self.Pc/self.Pe
        
        self.coolantTempStart = 300
        self.Kwall = 0
        self.Kc = 0
        self.chamberWallThickness = 0
        self.channelHeight = 0
        self.channelWallThickness = 0
        self.numChannels = 0
        self.temp_step = 1
        self.coefThermEx = 0
        self.youngMod = 0
        
        self.currFuel='Ethanol'
        self.currOx='O2'
        self.currMono='H2O2'
        self.oxCooled=False
        
        ################# ADD NEW MONOPROP   ########################################

        # Use code in monoprop mode
        # monoprop mode WILL override bi prop mode
        self.monoMode=0

        # 1 = new monoprop WILL be used for analysis, 0 = new monoprop will NOT be used
        # monoprop defaults to peroxide
        self.newMonoBool=0

        # name for new fuel
        self.newMonoName = ''

        # mono_component_name, elemental_composition, percent_weight
        # add as many rows as needed
        self.newMono = [['H2O2(L)'   ,'H 2 O 2'      ,    '100']]
        
        #################     ADD NEW FUEL    #######################################

        # 1 = new fuel WILL be used for analysis, 0 = new fuel will NOT be used
        # fuel defaults to pure ethanol
        self.newFuelBool=0

        # name for new fuel
        self.newFuelName = '75_ETH'

        # fuel_component_name, elemental_composition, percent_weight, chemical name
        # add as many rows as needed
        self.newFuel = [['C2H5OH(L)' ,'C 2 H 6 O 1'  ,    '75',          'Ethanol'],
                        ['H2O(L)'    ,'H 2 O 1'      ,    '25',          'Water']]
        
        #################     ADD NEW OXIDIZER    #######################################

        # 1 = new fuel WILL be used for analysis, 0 = new fuel will NOT be used
        # fuel defaults to liquid oxygen
        self.newOxBool=0

        # name for new ox
        self.newOxName = 'GelN2O4'

        # ox_component_name, elemental_composition, percent_weight
        # add as many rows as needed
        self.newOx = [['N2O4(L)'   ,'N 2 O 4'      ,    '96.5'],
                      ['SiO2'      ,'Si 1 O 2'     ,    '3.5']]


        
    def getCeaProperties(self):
        
        if self.newMonoBool:
            
            # create new mono card for rocketCEA and add together the prop weight percentages
            mono_weight_percent=0
            self.newMonoCard=''
            for i in range(len(self.newMono)):
                mono_weight_percent = mono_weight_percent + float(self.newMono[i][2])
                self.newMonoCard=self.newMonoCard+'name '+self.newMono[i][0]+' '+self.newMono[i][1]+' wt%='+self.newMono[i][2]+' '
    
            # throw error if components do not add to 100% weight
            if mono_weight_percent != 100:
                raise ValueError('NEW PROPELLANT COMPONENTS DO NOT ADD TO 100%')


        if self.newFuelBool:
            
            # create new fuel card for rocketCEA and add together the fuel weight percentages
            fuel_weight_percent=0
            self.newFuelCard=''
            for i in range(len(self.newFuel)):
                fuel_weight_percent = fuel_weight_percent + float(self.newFuel[i][2])
                self.newFuelCard=self.newFuelCard+'fuel '+self.newFuel[i][0]+' '+self.newFuel[i][1]+' wt%='+self.newFuel[i][2]+' '
    
            # throw error if components do not add to 100% weight
            if fuel_weight_percent != 100:
                raise ValueError('NEW FUEL COMPONENTS DO NOT ADD TO 100%')
            

        if self.newOxBool:
            
            # create new fuel card for rocketCEA and add together the fuel weight percentages
            ox_weight_percent=0
            self.newOxCard=''
            for i in range(len(self.newOx)):
                ox_weight_percent = ox_weight_percent + float(self.newOx[i][2])
                self.newOxCard=self.newOxCard+'oxid '+self.newOx[i][0]+' '+self.newOx[i][1]+' wt%='+self.newOx[i][2]+' '
    
            # throw error if components do not add to 100% weight
            if ox_weight_percent != 100:
                raise ValueError('NEW OXIDIZER COMPONENTS DO NOT ADD TO 100%')
        
        # monoprop mode overrides bi prop mode
        if self.monoMode:
            # if a new mono prop was created, use that for CEA
            if self.newMonoBool:
                add_new_propellant(self.newMonoName, self.newMonoCard)
                self.currMono=self.newMonoName
            self.rocketcea_eng_obj = CEA_Obj(propName=self.currMono)
            self.eps = self.rocketcea_eng_obj.get_eps_at_PcOvPe(self.Pc_psi, self.OF, self.PcOvPe)
            self.Tc = self.rocketcea_eng_obj.get_Tcomb(self.Pc_psi, self.OF) * 0.555556
            self.MM = self.rocketcea_eng_obj.get_Chamber_MolWt_gamma(self.Pc_psi, self.OF, self.eps)[0] / 1000
            self.gamma = self.rocketcea_eng_obj.get_Chamber_MolWt_gamma(self.Pc_psi, self.OF, self.eps)[1]
            self.cea_Me = self.rocketcea_eng_obj.get_MachNumber(self.Pc_psi, self.OF, self.eps)
            self.cea_Isp = self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, self.OF, self.eps)[0]
            self.Cstar = self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, self.OF, self.eps)[1] * 0.3048
            return
        
        # if new oxidizer was created, use that for CEA
        if self.newOxBool:
            add_new_oxidizer(self.newOxName, self.newOxCard)
            self.currOx=self.newOxName
            
        # if new fuel was created, use that for CEA
        if self.newFuelBool:
            add_new_fuel(self.newFuelName, self.newFuelCard)
            self.currFuel=self.newFuelName

        self.rocketcea_eng_obj = CEA_Obj(oxName=self.currOx, fuelName=self.currFuel)
        self.eps = self.rocketcea_eng_obj.get_eps_at_PcOvPe(self.Pc_psi, self.OF, self.PcOvPe)
        self.Tc = self.rocketcea_eng_obj.get_Tcomb(self.Pc_psi, self.OF) * 0.555556
        self.MM, self.gamma = self.rocketcea_eng_obj.get_Chamber_MolWt_gamma(self.Pc_psi, self.OF, self.eps)
        self.cea_Me = self.rocketcea_eng_obj.get_MachNumber(self.Pc_psi, self.OF, self.eps)
        self.cea_Isp = self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, self.OF, self.eps)[0]
        self.Cstar = self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, self.OF, self.eps)[1] * 0.3048
        return
    
    def plotIspVOF(self):
        
        self.getCeaProperties()
        
        eps = [self.eps-5, self.eps, self.eps+5]
        for e in eps:
            ispArr = []

            MR = self.OF_low
            mrArr = []
            while MR < self.OF_high:
                ispArr.append(self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, MR, e)[0])
                mrArr.append(MR)
                MR += 0.05
            plt.plot(mrArr, ispArr, label='AreaRatio %d'%e)
            
        plt.title('Isp v. MR')
        plt.xlabel('Mixture Ratio')
        plt.ylabel('Isp')
        plt.legend(loc='best')
        plt.gcf().set_dpi(300)
        plt.show()
    
    def engineDimensions(self):
        
        self.getCeaProperties()
        
        R=8.314
        
        self.MM = self.MM/1000
        self.Rbar = R/self.MM
        pi = math.pi

        # Calculations
        self.To = self.Tc # stagnation temp K assumed to be adiabatic flame temp in chamber
        self.Pt = self.Pc*(1+((self.gamma-1)/2))**(-self.gamma/(self.gamma-1)) # pressure at throat Pa
        self.Tt = self.Tc*(1/(1+((self.gamma-1)/2))) # Throat temp K
        self.Me = math.sqrt((2/(self.gamma-1))*(((self.Pc/self.Pa)**((self.gamma-1)/self.gamma))-1)) # Mach at exit
        self.Po = self.Pt*(1/(1+(((self.gamma-1)/2)))**(-self.gamma/(self.gamma-1))) # stagnation pressure Pa using Pt
        self.Ueq = math.sqrt(((2*self.gamma)/(self.gamma-1))*self.Rbar*self.To*(1-((self.Pe/self.Po)**((self.gamma-1)/self.gamma)))) # exit velocity m/s
        self.mdot = self.thrust/self.Ueq # kg

        self.mdotFuel = self.mdot / (self.OF + 1)
        self.mdotOx = self.mdot - self.mdotFuel

        bigGamma = math.sqrt(self.gamma*((2/(self.gamma+1))**((self.gamma+1)/(self.gamma-1))))
        self.At = (self.mdot*math.sqrt(self.Rbar*self.To))/(self.Po*bigGamma) # area of the throat m2
        self.Ae = (self.At/self.Me)*((1+((self.gamma-1)/2)*(self.Me**2))/((self.gamma+1)/2))**((self.gamma+1)/(2*(self.gamma-1))) # area of the exit m2

        self.Isp = self.Ueq/9.81

        self.cp = (self.gamma/(self.gamma-1))*self.Rbar

        self.Dt = math.sqrt((4/pi)*self.At)
        self.De = math.sqrt((4/pi)*self.Ae)
        self.Dc = self.Dt*3


        self.Ac = (math.pi * self.Dc**2)/4
        self.Vc = self.Lstar*self.At
        self.Lc = self.Vc / (1.1 * self.Ac)

        # solving for chamber diameter iteratively
        self.theta = math.tan(math.radians(self.beta))
        i = 0
        while i < 10:
            self.Dc = math.sqrt((self.Dt**3 + (24/math.pi)*self.theta*self.Vc) / (self.Dc + 6*self.theta*self.Lc))
            i = i + 1

        # nozzle geometry curves
        self.Re = self.De / 2
        self.Rt = self.Dt / 2
        self.Rc3 = 0.382 * self.Rt
        self.Rc2 = 1.5 * self.Rt
        self.Rc1 = self.Rc2 * 2
        self.Rc = self.Dc / 2

        # distances from throat (negative)
        tan = math.tan(math.radians(self.alpha))
        # sin = math.sin(math.radians(self.alpha))
        cos = math.cos(math.radians(self.alpha))
        self.throat = 0
        self.endStraight = -1*((self.Rc2 - (self.Rc2*cos)) / math.tan(math.radians(self.alpha/2)))*1000
        conv_vert_portion = self.Rc - self.Rt - (self.Rc2 - (self.Rc2*cos)) - (self.Rc1 - (self.Rc1*cos))
        self.endRc1 = self.endStraight - (conv_vert_portion / tan)*1000
        self.endChamber = self.endRc1 - ((self.Rc1 - (self.Rc1*cos)) / math.tan(math.radians(self.alpha/2)))*1000
        self.startChamber = self.endChamber - self.Lc*1000
        self.endRc3 = 1000*(self.Rc3 * math.cos(math.radians(90-self.beta))) / (math.tan(math.radians(self.beta)))
           
    def propFlowRates(self):
        
        self.getCeaProperties()
        self.engineDimensions()
        
        if self.monoMode == 1:
            raise Exception('PROPELLANT RATES NOT AVAILABLE IN MONOPROP MODE\nSET ''monoMode = 0'' AND USE BI PROPELLANT MODE')
        
        if len(self.newFuel) == 2:
        
            # still needed for CoolProp
            # fuel_comp_1 = float(self.newFuel[0][2]) / 100 # fraction of 1
            fuel_comp_2 = float(self.newFuel[1][2]) / 100
            comp_1 = self.newFuel[0][1]
            comp_1.replace(' ', '')
            comp_1 = Formula(comp_1)
            fuel_comp_1_MM = comp_1.mass
            comp_2 = self.newFuel[1][1]
            comp_2.replace(' ', '')
            comp_2 = Formula(comp_2)
            fuel_comp_2_MM = comp_2.mass
            # waterMM = 18.02 # g/mol
            # ethanolMM = 46.07 # g/mol
    
            # find values for CoolProp in heat transfer calcs
            fuel_comp_1_frac = round(((1-fuel_comp_2) * fuel_comp_1_MM) / (((1-fuel_comp_2) * fuel_comp_1_MM) + (fuel_comp_2 * fuel_comp_2_MM)),4)
            fuel_comp_2_frac = round(1 - fuel_comp_1_frac, 4)
            # remove the 0 before the decimal for CoolProp
            fuel_comp_1_frac = str(fuel_comp_1_frac)[1:]
            fuel_comp_2_frac = str(fuel_comp_2_frac)[1:]
    
            self.fuelComposition = self.newFuel[0][3]+'['+fuel_comp_1_frac+']&'+self.newFuel[1][3]+'['+fuel_comp_2_frac+']'
            
        else:
            
            self.fuelComposition = self.newFuel[0][3]
        
        if (self.currOx != 'O2'):
            print("Ox cooling only available for liquid oxygen")
            
        else:
            self.oxComposition = 'Oxygen'
            
        # density and volume/sec calcs from tanks
        self.rhoLOX = 1141 # kg/m3
        self.rhoFuel = PropsSI('D', 'T', 300, 'P', 101325, self.fuelComposition) # kg/m3

        self.VdotLOX = (self.mdotOx / self.rhoLOX) * 1000 # L/s
        self.VdotFuel = (self.mdotFuel / self.rhoFuel) * 1000 # L/s

        self.volumeFuel = self.VdotFuel * self.burnTime # L
        self.volumeLOX = self.VdotLOX * self.burnTime # L
        self.massFuel = self.volumeFuel * self.rhoFuel / 1000 # kg
        self.massLOX = self.volumeLOX * self.rhoLOX / 1000 # kg
        self.massProp = self.massFuel + self.massLOX
        
    def printChamberInfo(self):
        
        self.getCeaProperties()
        self.engineDimensions()
        self.propFlowRates()
        
        print("Ae/At: ",round(self.Ae/self.At,3))
        print("Chamber pressure (kPa): ",round(self.Pc/1000,2))
        print("Chamber temp (K): ", round(self.Tc,2))
        print("Thrust (N): ",round(self.thrust,2))
        print("Pressure at throat (kPa): ",round(self.Pt/1000,2))
        print("Temp at throat (K): ",round(self.Tt,2))
        print("Mach at exit: ",round(self.Me,4))
        print("CEA Mach at exit: ",round(self.cea_Me,4))
        print("Stag pressure (kPa): ",round(self.Po/1000,2))
        print("Stag temp (K): ",round(self.To,2))
        print("Exit velocity (m/s): ",round(self.Ueq,2))
        # print("Exit temp (K): ",round(Te,2))
        print("mdot (kg/s): ",round(self.mdot,3))
        print("fuel mdot (kg/s): ",round(self.mdotFuel,3))
        print("ox mdot (kg/s): ",round(self.mdotOx,3))
        print("fuel Vdot (L/s): ",round(self.VdotFuel,4))
        print("ox Vdot (L/s): ",round(self.VdotLOX,4))
        print("fuel volume (L): ",round(self.volumeFuel,2))
        print("ox volume (L): ",round(self.volumeLOX,2))
        print("fuel mass (kg): ",round(self.massFuel,2))
        print("ox mass (kg): ",round(self.massLOX,2))
        print("propellant mass (kg): ",round(self.massProp,2))
        print("Chamber volume (m3): ",round(self.Vc,6))
        print("Chamber length (m): ",round(self.Lc,6))
        print("Diameter of chamber (m): ",round(self.Dc,6))
        print("Area of throat (m2): ",round(self.At,6))
        print("Diameter of throat (m): ",round(self.Dt,6))
        print("Area of exit (m2): ",round(self.Ae,6))
        print("Diameter of exit (m): ",round(self.De,6))
        print("c* (m/s): ",round(self.Cstar,2))
        print("Cp: ",round(self.cp,2))
        print("\nEngine Profile:")
        print("Rc: ",round(self.Rc*1000,6))
        print("Rc1: ",round(self.Rc1*1000,6))
        print("Rc2: ",round(self.Rc2*1000,6))
        print("Rt: ",round(self.Rt*1000,6))
        print("Rc3: ",round(self.Rc3*1000,6))
        print("Re: ",round(self.Re*1000,6))
        print("Isp: ",round(self.Isp,2))
        print("cea_Isp: ",round(self.cea_Isp,2))
        
    def engineContour(self):
    # given the geometry of the engine this will return the profile of the engine
    # as well as plot the profile. Nozzle geometry is found using Method of Characteristics
    # startChamber is negative dist from throat that the chamber starts in mm
    # endChamber is neg dist from throat that the chamber ends in mm
    # endRc1 is neg dist from throat that first conv curve ends in mm
    # endStraight is neg dist from throat that straight conv section ends in mm
    # throat is location of the throat, assumed to be at origin in mm
    # endRc3 is pos dist from throat that first div curve ends in mm
    # Rc is radius of chamber in m
    # Rt is radius of throat in m
    # Re is radius of exit in m
    # Me is design mach number at exit
    # n is number of characteristics you want to calculate for nozzle
        
        self.getCeaProperties()
        self.engineDimensions()
        
        # convert to m
        Rc = self.Rc * 1000
        Rt = self.Rt * 1000
        Re = self.Re * 1000
        
        # these no longer get used
        # thetaN = math.radians(40) # degrees to rad
        # thetaE = math.radians(20) # degrees to rad
        
        # sets the start point of the engine and initializes lists
        self.engineX = [self.startChamber]
        self.engineY = []
        
        x = self.engineX[0]
        
        i=0
        # solves for the curves that make up the engine up to the throat in 0.01 mm intervals
        while True:

            # chamber
            # constant radius (cylindrical)
            if self.engineX[i] <= self.endChamber:
                self.engineY.append(Rc)
                
            # first curve in conv section
            # this has double the radius of the curve leading into the nozzle
            if self.endChamber < self.engineX[i] <= self.endRc1:
                self.engineY.append(math.sqrt(abs(((3*Rt)**2) - (self.engineX[i] - self.endChamber)**2)) + (Rc - 3 * Rt))
            
            # straight conv section
            # straight line connecting the two converging section curves
            if self.endRc1 < self.engineX[i] <= self.endStraight:
                
                # finding the slope of the line between the end and start points
                # (y - y1) = m(x - x1) + b
                x1 = self.endRc1
                y1 = math.sqrt(abs(((3*Rt)**2) - (self.endRc1 - self.endChamber)**2)) + (Rc - 3 * Rt)
                
                x2 = self.endStraight
                y2 = (Rt + 1.5 * Rt) - math.sqrt(abs(((1.5*Rt)**2) - (self.endStraight)**2))
                
                # slope of the line
                m = (y1 - y2) / (x1 - x2)
                
                b = y1 - m * x1
                
                # y = mx + b
                self.engineY.append(m * self.engineX[i] + b)
                
            # curve into throat
            # curve based on 1.5 times Rt
            if self.endStraight < self.engineX[i] <= (self.throat - self.throat):
                self.engineY.append((Rt + 1.5 * Rt) - math.sqrt(abs(((1.5*Rt)**2) - (self.engineX[i])**2)))
                
            # once we reach the throat, break the loop and get nozzle geometry from MOC
            if self.engineX[i] == 0:
                break
            
        #          This part is replaced with the Method of Characteristics   
        #
        ###############################################################################
        ###############################################################################
                
            # circular curve after throat
        #     # curve based on 0.382 time Rt
        #     if (throat - throat) < X[i] <= endRc3:
        #         Y.append((Rt + 0.382 * Rt) - math.sqrt(abs(((0.382*Rt)**2) - (X[i])**2)))
                
        #     # parabolic nozzle
        #     # takes thetaN, Rn, thetaE, and Re and returns the parabolic curve to fit
        #     if endRc3 < X[i]:
                
        #         # finding Rn (the last point in the curve after the throat)
        #         Rn = (Rt + 0.382 * Rt) - math.sqrt(abs(((0.382*Rt)**2) - (endRc3)**2))
        
                
        # #           Bit off, take 1        
        # #############################################################################        
        #         # # linear sys of eqns used to solve for constants for the parabola
        #         # A = np.array([[2*Rn, 1, 0],[2*Re, 1, 0],[Rn**2, Rn, 1]])
        #         # B = np.array([[(1/(math.tan(thetaN)))],[(1/(math.tan(thetaE)))],[endRc3]])
                
        #         # # constants seen in x = ay^2 + by + c
        #         # a, b, c = np.linalg.inv(A).dot(B)
                
        #         # # quadratic formula to solve for y (the radius)
        #         # F = (b**2) - 4 * a * (c - X[i])
        #         # Y.append(((-1*b) + math.sqrt(F)) / (2 * a))
        # #############################################################################
        
        # #           Take , better but needs points from CAD
        # #############################################################################
        #         # # linear sys of eqns used to solve for constants for the parabola
        #         # # uses vertex from CAD in eqn1, Rn in eqn2, and exit in eqn3
        #         # A = np.array([[8.6257**2, 8.6257, 1],[Rn**2, Rn, 1],[Re**2, Re, 1]])
        #         # B = np.array([[-5.585],[endRc3],[52.356]])
                
        #         # # constants seen in x = ay^2 + by + c
        #         # a, b, c = np.linalg.inv(A).dot(B)
                
        #         # # quadratic formula to solve for y (the radius)
        #         # F = (b**2) - 4 * a * (c - X[i])
        #         # Y.append(((-1*b) + math.sqrt(F)) / (2 * a))
        # #############################################################################
                
                
                
        #         # if the radius of the exit is reached, exit the loop and give the 
        #         # values for the nozzle curve to X and Y to be plot
        #         if Y[i] >= Re:                
        #             break
        
        ###############################################################################
        ###############################################################################
                
            x = round(x + 0.01,2) 
            self.engineX.append(round(x,2))
            i = i + 1
        
        # get nozzle geometry from MOC_nozzle and 
        # append the points onto X and Y       
        self.nozX, self.nozY = MOC.nozzle(self.Me, self.num_characteristics, Rt, Re)
        self.nozXContour = []
        self.nozYContour = []
        self.nozZContour = []
        for i in range(len(self.nozX)):
            self.engineX.append(self.nozX[i])
            self.engineY.append(self.nozY[i])
            self.nozXContour.append(self.nozX[i])
            self.nozYContour.append(self.nozY[i])
            
        self.nozZContour = [0] * len(self.nozXContour)
        
    def saveNozGeo(self, filePath):
    # saves nozzle geometry in a txt file that can be imported into CAD 
    # software to model the engine
        
        self.engineContour()

        j = True
        k=1
        while j == True:
            filename = filePath + '/nozzle_geo' + str(k) + '.txt'
            j = path.exists(filename)
            k = k + 1
            
        df = pd.DataFrame(data={"xVals":self.nozXContour, "yVals":self.nozYContour, "zVals":self.nozZContour})
        df.to_csv(filename, sep=' ', index = False, header=False)
        
    def saveChamberGeo(self, filePath):
    # saves nozzle geometry in a txt file that can be imported into CAD 
    # software to model the engine
        
        self.engineContour()

        j = True
        k=1
        while j == True:
            filename = filePath + '/chamber_geo' + str(k) + '.txt'
            j = path.exists(filename)
            k = k + 1
            
        self.engineZ = [x*0 for x in self.engineY]
        df = pd.DataFrame(data={"xVals":self.engineX, "yVals":self.engineY, "zVals":self.engineZ})
        df.to_csv(filename, sep=' ', index = False, header=False)
        
    def heatTransfer(self):
        
        if self.coolantTempStart == 0:
            raise ValueError('INPUT VALUE FOR coolandTempStart THIS IS THE TEMPERATURE (K) OF THE FUEL IN THE TANK')
        if self.Kwall == 0:
            raise ValueError('INPUT VALUE FOR Kwall THIS IS THE THERMAL CONDUCTIVITY (W/m*K) OF THE CHAMBER WALL')
        if self.Kc == 0:
            raise ValueError('INPUT VALUE FOR Kc THIS IS THE THERMAL CONDUCTIVITY (W/m*K) OF THE FUEL')
        if self.chamberWallThickness == 0:
            raise ValueError('INPUT VALUE FOR chamberWallThickness THIS IS THE THICKNESS (mm) OF THE CHAMBER WALL')
        if self.channelHeight == 0:
            raise ValueError('INPUT VALUE FOR channelHeight THIS IS THE DEPTH (mm) OF THE COOLING CHANNELS')
        if self.channelWallThickness == 0:
            raise ValueError('INPUT VALUE FOR channelWallThickness THIS IS THE WALL THICKNESS (mm) BETWEEN COOLING CHANNELS')
        if self.numChannels == 0:
            raise ValueError('INPUT VALUE FOR numChannels THIS IS THE NUMBER OF COOLING CHANNELS AROUND THE CHAMBER')

        self.engineContour()
        
        if self.monoMode == 1:
            raise Exception('REGEN COOLING NOT AVAILABLE IN MONOPROP MODE\nSET ''monoMode = 0'' AND USE BI PROPELLANT MODE')
        
        self.propFlowRates()
        
        X = [a/1000 for a in self.engineX]
        Y = [b/1000 for b in self.engineY]
        Y = Y[::-1]
        X = X[::-1]
        
        # heat transfer stuff
        muE = 6.5974511*(10**-5) # mu at throat kg/ms cp of air at Te of 2124 K
        w = 0.6 # some constant for the HX eqn I really don't know
        r = 0.9 # some other constant fr dunno
        
        self.Te = self.Tc / (1 + ((self.gamma-1)/2)*self.Me**2) # exit temp (works for freestream temp too)
        
        self.oTo = [0] * len(X) # K location temp averaged from local temp T and temp of wall Tw
        self.hg = [1] * len(X) # convection coef for combustion gases
        self.qw =[0] * len(X) # W/m2 heat flux thru wall
        self.qw_safer = [0] * len(X) # heat flux with added safety factor
        self.M =[0] * len(X) # Mach number along flow thru nozzle
        # area=[0] * len(X)
        self.T = [0] * len(X) # K local temp along flow
        self.Tr = [0] * len(X) # K local recovery temp 
        self.Tw1 = [350] * len(X)
        self.Tw2 = [0] * len(X)
        self.thermalStress = [0] * len(X)
        self.TCoolant = [self.coolantTempStart] * (len(X)+1)
        
        self.hc = [1] * len(X)
        self.coolV = [0] * len(X)
        
        self.channelWidth = [0] * len(X)
        
        xStep = X[1] - X[2] # m x step along engine
        pastThroat = False
        circum = 0

        numChannelsChamber = self.numChannels / 1
        
        # k = ((self.gamma+1)/(2*(self.gamma-1)))
        const1 = (0.026/(self.Dt**0.2))*((self.Pc/self.Cstar)**0.8)
        
        ####################################################
        ####################################################
        
        # this shit broke (no longer a jank fix)

        m1 = (1 - self.Me) / (1 - (self.Ae/self.At))
        b1 = 1-abs(m1)
        
        m2 = (0.05 - 1) / ((self.Ac/self.At) - 1)
        b2 = 1+abs(m2)
        
        ####################################################
        ####################################################
        
        self.CpDict = {}
        self.rhoDict = {}
        self.muDict = {}
        
        tempCp = 0
        tempMu = 0
        tempRho = 0
        tempArray = np.arange(self.coolantTempStart, 500, self.temp_step)
        
        for q in range(len(tempArray)):
            tempArray[q] = round(tempArray[q],2)
        
        if self.currOx == 'O2' and self.oxCooled:
            # ethanol and water by mole fractions, change this to not be hard coded
            rho = PropsSI('D', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
            mu = PropsSI('V', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
            Cp = PropsSI('C', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
            print(rho)
            print(mu)
            print(Cp)
        elif self.currOx != 'O2' and self.oxCooled:
            print("ONLY ABLE TO COOL WITH LIQUID OXYGEN")
            return
        elif not self.oxCooled:
            # ethanol and water by mole fractions, change this to not be hard coded
            rho = PropsSI('D', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
            mu = PropsSI('V', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
            Cp = PropsSI('C', 'T', tempArray, 'P', self.Pc, self.fuelComposition)
        
        for t in range(len(tempArray)):
            self.CpDict[tempArray[t]] = Cp[t]
            self.muDict[tempArray[t]] = mu[t]
            self.rhoDict[tempArray[t]] = rho[t]
             
        for i in range(len(X)):
            
            # set flag for when the calculations pass nozzle to allow for
            # number of channels to change along chamber
            if X[i] == 0:
                pastThroat = True
                
        

            while (abs(self.hg[i]*(self.Tr[i]-self.Tw1[i])) - (self.hc[i]*(self.Tw2[i] - self.TCoolant[i])) > 1):
                A = np.pi*(Y[i]**2)
                
                # # this shit is broke and fucks everything up
                # #############################################
                # if X[i] >= 0:
                #     M[i] = ((((((A/Astar)**(1/k)) * (1 + ((gamma-1)/2))) - 1)/((gamma-1)/2))**(1/(1+k)))
                # if X[i] < 0:
                #     M[i] = 0.5
                # #############################################
                
                # # fix testing
                # #############################################

                # M[i] = (((((A/Astar) * ((1 + ((gamma-1)/2)))**k) - 1)/(((gamma-1)/2)**k))**(1/(k - 1)))
                # # M[i] = (((((A/Astar)/((2/(gamma+1))**k)) - 1) / (((gamma-1)/2)**k))**(1/(2*k)))
                # # if X[i] == 0:
                # #     print(M[i])

                #############################################
                
                # fix testing x2
                # this assumes area ratio vs Ma is linear
                #############################################

                if X[i] >= 0:
                    self.M[i] = (m1 * (A/self.At)) + b1
                    
                if X[i] < 0:
                    self.M[i] = (m2 * (A/self.At)) + b2
                    
                #############################################
                
                self.T[i] = (self.Tc / (1+ ((self.gamma-1)/2) * self.M[i]**2))
                self.Tr[i] = (self.T[i]*(1+((self.gamma-1)/2)*r*(self.M[i]**2)))
                self.oTo[i] = ((self.T[i]+self.Tw1[i])/2)
                
                
                # the arith. mean temp has negligible impact from error
                # when Ma is higher than low supersonic
                # if M[i] < 1.1:
                #     oTo[i] = ((T[i]+Tw1[i])/2)
                    
                # # this limits the mean arth. temp to lower mach numbers
                # # not known yet if this is a reasonable fix or not
                # else:
                #     Tss = (Tc / (1 + ((gamma-1)/2) * 1.1**2))
                #     oTo[i] = ((Tss+Tw1[i])/2)                
                
                self.hg[i] = (const1*((self.Dt/(2*Y[i]))**1.8)*self.cp*(muE**0.2)*((self.Te/self.oTo[i])**(0.8-0.2*w)))
                
                self.qw[i] = (self.hg[i]*(self.Tr[i]-self.Tw1[i]))
                
                # use 0.45 if the oTo is not compensated for high mach numbers
                self.qw_safer[i] = (self.qw[i]*0.66)
                
                # temp of wall on the coolant side in order 
                # for the combustion side to stay at temp we want
                self.Tw2[i] = (self.Tw1[i] - (self.qw_safer[i]/(self.Kwall/(self.chamberWallThickness/1000))))
                
                tempCp = self.CpDict[round(self.TCoolant[i],0)]
                tempMu = self.muDict[round(self.TCoolant[i],0)]
                tempRho = self.rhoDict[round(self.TCoolant[i],0)]

                circum = 2*(Y[i] + (self.chamberWallThickness/1000))*np.pi
                totalChannelSpace = circum - (self.numChannels*(self.channelWallThickness/1000))
                self.channelWidth[i] = (totalChannelSpace/self.numChannels)
            
                # Prandtl number 
                Pr = (tempMu * tempCp) / self.Kc
                
                # the area of the coolant passage
                Ap = (self.channelHeight/1000)*self.channelWidth[i]
                
                # the wetted perimeter of the coolant passage
                Pw = 2*(self.channelHeight/1000) + 2*self.channelWidth[i]
                
                # the 'diameter' of the passage
                D = (4*Ap) / Pw
                
                # velocity in cooling passages
                u = (self.mdot / self.numChannels) / (tempRho * Ap)
                self.coolV[i] = (u)
                
                # Renolds num in cooling passages
                Re = (D * u * tempRho) / tempMu
                
                # calculating hc based on channel geometry
                self.hc[i] = (0.023 * tempCp * ((self.mdot/self.numChannels)/Ap) * (Re**(-0.2)) * (Pr**(2/3)))
                
                self.Tw1[i] = self.Tw1[i] + self.temp_step
                
            self.TCoolant[i+1] = ((self.qw_safer[i] * (self.channelWidth[i]*xStep))/((Ap*tempRho*u)*tempCp)) + self.TCoolant[i]
            self.thermalStress[i] = ((self.youngMod * 1000) * self.coefThermEx * (abs(self.Tw1[i] - self.Tw2[i])))

            # once the calculations pass the nozzle and reach the radius
            # of the chamber, lower the number of channels
            if pastThroat == True and Y[i] == self.Rc:
                self.numChannels = numChannelsChamber


        self.qw = self.qw[::-1]
        self.qw_safer = self.qw_safer[::-1]
        self.hc = self.hc[::-1]
        self.Tw1 = self.Tw1[::-1]
        self.Tw2 = self.Tw2[::-1]
        self.thermalStress = self.thermalStress[::-1]
        self.channelWidth = self.channelWidth[::-1]
        self.channelWidth = [x*1000 for x in self.channelWidth]
        self.TCoolant = self.TCoolant[::-1]
        self.coolV = self.coolV[::-1]
        self.M = self.M[::-1]
        self.TCoolant.pop(-1)
        
        print('Max Thermal Stress in Nozzle (MPa): ',round(max(self.thermalStress),0))
        
    def plotEngineContour(self):
        
        self.engineContour()
        
        x=np.arange(-200,150,10, dtype=(int))
        y=np.arange(0,100,5,dtype=(int))
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.engineY)
        ax.set_aspect(1.0)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.yticks(y)
        plt.title('Engine Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
        
    def plotNozzleContour(self):
        
        self.engineContour()
        
        x=np.arange(-200,150,10, dtype=(int))
        y=np.arange(0,100,5,dtype=(int))
        
        fig, ax = plt.subplots(figsize=(25,10))
        plt.plot(self.nozX, self.nozY,"b")
        ax.set_aspect(1.0)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.yticks(y)
        plt.title('Nozzle Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
        
    def heatTransferPlots(self):
        
        self.engineContour()
        
        if self.monoMode == 1:
            raise Exception('HEAT TRANSFER PLOTS NOT AVAILABLE IN MONOPROP MODE\nSET ''monoMode = 0'' AND USE BI PROPELLANT MODE')
        
        self.propFlowRates()
        self.heatTransfer()
        
        x=np.arange(-200,150,10, dtype=(int))
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.M)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Mach Number Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Mach Num')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.qw,label='qw')
        plt.plot(self.engineX,self.qw_safer,label='qw - safer')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Flux Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Flux (W/m^2)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.Tw1,label='Hot Side')
        plt.plot(self.engineX,self.Tw2,label='Fuel Side')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Wall Temp')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.thermalStress,label='Thermal Stress')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Thermal Stress')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Thermal Stress (kPa)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.coolV)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Coolant Velocity')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Velocity (m/s)')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.hc,label='hc')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Transfer Coefficient Thru Cooling Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Transfer Coefficient (W/(m^2K))')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.TCoolant)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Temp of Coolant Along Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.channelWidth)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Width of Coolant Passages Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
        
    def allPlots(self):
        
        self.engineContour()
        
        x=np.arange(-200,150,10, dtype=(int))
        y=np.arange(0,100,5,dtype=(int))
            
        fig, ax = plt.subplots(figsize=(25,10))
        plt.plot(self.nozX, self.nozY,"b")
        ax.set_aspect(1.0)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.yticks(y)
        plt.title('Nozzle Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
            
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.engineY)
        ax.set_aspect(1.0)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.yticks(y)
        plt.title('Engine Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
        
        if self.monoMode == 1:
            return
        
        self.propFlowRates()
        self.heatTransfer()
        self.plotIspVOF()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.M)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Mach Number Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Mach Num')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.qw,label='qw')
        plt.plot(self.engineX,self.qw_safer,label='qw - safer')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Flux Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Flux (W/m^2)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.Tw1,label='Hot Side')
        plt.plot(self.engineX,self.Tw2,label='Fuel Side')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Wall Temp')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.thermalStress,label='Thermal Stress')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Thermal Stress')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Thermal Stress (kPa)')
        plt.legend(loc='upper right')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.coolV)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Coolant Velocity')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Velocity (m/s)')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.hc,label='hc')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Transfer Coefficient Thru Cooling Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Transfer Coefficient (W/(m^2K))')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.TCoolant)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Temp of Coolant Along Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.show()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.channelWidth)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Width of Coolant Passages Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.show()
        
    def savePlots(self, dpi):
        
        self.engineContour()
        
        x=np.arange(-160,100,10, dtype=(int))
        y=np.arange(0,100,5,dtype=(int))
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.engineY)
        ax.set_aspect(1.0)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.yticks(y)
        plt.title('Engine Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.savefig('engine_contour.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.nozX, self.nozY,"b")
        ax.set_aspect(1.0)
        plt.title('Nozzle Geometry')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.savefig('nozzle_contour.png', dpi=dpi)
        
        if self.monoMode == 1:
            return
        
        self.propFlowRates()
        self.heatTransfer()
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.M)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Mach Number Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Mach Num')
        plt.savefig('mach_num_thru_engine.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.qw,label='qw')
        plt.plot(self.engineX,self.qw_safer,label='qw - safer')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Flux Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Flux (W/m^2)')
        plt.legend(loc='upper right')
        plt.savefig('heat_flux.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.Tw1,label='Hot Side')
        plt.plot(self.engineX,self.Tw2,label='Fuel Side')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Wall Temp')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.legend(loc='upper right')
        plt.savefig('wall_temp.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.thermalStress,label='Thermal Stress')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Thermal Stress')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Thermal Stress (kPa)')
        plt.legend(loc='upper right')
        plt.savefig('thermal_stress.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.coolV)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Coolant Velocity')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Velocity (m/s)')
        plt.savefig('coolant_velocity.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.hc,label='hc')
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Heat Transfer Coefficient Thru Cooling Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Heat Transfer Coefficient (W/(m^2K))')
        plt.savefig('hx_coef.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.TCoolant)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Temp of Coolant Along Passages')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Temperature (K)')
        plt.savefig('coolant_temp.png', dpi=dpi)
        
        fig, ax = plt.subplots(figsize=(25,10), dpi=300)
        plt.plot(self.engineX,self.channelWidth)
        plt.grid(color='k',alpha=0.25)
        plt.xticks(x)
        plt.title('Width of Coolant Passages Along Engine')
        plt.xlabel('Distance Along Engine Axially (mm)')
        plt.ylabel('Distance (mm)')
        plt.savefig('passage_width.png', dpi=dpi)
        
    def doItAll(self,dpi,path):
        
        self.printChamberInfo()
        self.savePlots(dpi)
        self.saveNozGeo(path)
        self.plotIspVOF()