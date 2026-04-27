# -*- coding: utf-8 -*-
"""
Engine_Class.py

This module takes user input and returns heat transfer simulation plots,
performance plots, and geometric sizing for chemical rocket engines.

run program with input file or import class and use within other functions
to run directly:   python Engine_Class.py my_engine.toml [dpi] [output_path]

takes arguments
# argv[1]  : path to TOML input file            (required)
# argv[2]  : DPI for saved plots                (optional, default 150)
# argv[3]  : output directory for nozzle geo    (optional, default '.')

Raises:
    ImportError: checks for tomllib to be installed
    ValueError: checks if mono propellant combination adds to 100%
    ValueError: checks if fuel propellant combination adds to 100%
    ValueError: checks if oxidizer propellant combination adds to 100%
    Exception: propellant flow rates are unavailable in monopropellant mode
    ValueError: checks for starting temperature of propellant in tank
    ValueError: checks for the thermal conductivity of the chamber wall
    ValueError: checks for the thermal conductivity of the regen fluid
    ValueError: checks for the chamber wall thickness
    ValueError: checks for the regen channel height
    ValueError: checks for the regen channel wall thickness
    ValueError: checks for the number of regen channels
    Exception: regen cooling is unavailable in monopropellant mode
    Exception: heat transfer plots are unavailable in monopropellant mode

Output:
    generates a javasccript report of file type docx that captures all input data
    and the resulting outputs.
"""

from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import MOC_nozzle as MOC
import pandas as pd
import shutil
from os import path
import numpy as np
from molmass import Formula

# fluid_names provides canonical ↔ CEA/CoolProp/RocketProps name translation
# and a unified get_props() dispatcher (fluid_dict → rocketprops → CoolProp).
try:
    import importlib.util as _ilu, os as _fno
    _fn_path = _fno.path.join(_fno.path.dirname(_fno.path.abspath(__file__)), 'fluid_names.py')
    if _fn_path:
        _fn_spec = _ilu.spec_from_file_location('fluid_names', _fn_path)
        _fn_mod  = _ilu.module_from_spec(_fn_spec)
        _fn_spec.loader.exec_module(_fn_mod)
        import sys as _fns; _fns.modules['fluid_names'] = _fn_mod
        from fluid_names import translate, canonical, get_props as _get_fluid_props
        _HAS_FLUID_NAMES = True
except Exception as _fn_err:
    _HAS_FLUID_NAMES = False
    def translate(name, target): return name
    def canonical(name): return name

class Engine():

    def __init__(self, thrust=0.0, Pc=0.0, OF=0.0, height_of_optimization=0.0, burnTime=0.0,
                 input_file=None):
        """
        Args:
            thrust (float, optional): 
                Designed thrust of. Defaults to 0.
            Pc (float, optional): 
                Designed chamber pressure. Defaults to 0.
            OF (float, optional): 
                Designed mixture ratio (oxidizer / fuel) by mass. Defaults to 0.
            height_of_optimization (float, optional): 
                Altitude that the engine is perfectly expanded. Defaults to 0.
            burnTime (float, optional): 
                Burntime in seconds. Defaults to 0.
            input_file (_type_, optional): 
                Path to a .toml input file.  Any value in the file overrides the
                corresponding argument or default set below.  The five required
                parameters (OF, Pc, thrust, height_of_optimization, burnTime) may
                be omitted from the constructor call entirely when input_file is
                provided. Defaults to None.
        """

        # Required engine parameters (may be overridden by input_file)
        self.OF                   = OF
        self.Pc_psi               = Pc
        self.thrust_lbf           = thrust        # lbf value before conversion
        self.height_of_optimization = height_of_optimization
        self.burnTime             = burnTime

        self.engine_name          = ''
        self.num_characteristics  = 300
        self.Lstar                = 1
        self.alpha                = 30
        self.beta                 = 15

        self.Kwall                = 0
        self.Kc                   = 0
        self.chamberWallThickness = 0
        self.channelHeight        = 0
        self.channelWallThickness = 0
        self.numChannels          = 0
        self.temp_step            = 1
        self.coefThermEx          = 0
        self.youngMod             = 0

        self.currFuel             = 'Ethanol'
        self.currOx               = 'O2'
        self.currMono             = 'H2O2'
        self.oxCooled             = False
        self.oxTankTemp           = 90 # K
        self.fuTankTemp           = 300 # K
        
        if self.oxCooled:
            self.coolantTempStart = self.oxTankTemp
        else:
            self.coolantTempStart = self.fuTankTemp

        # OF sweep bounds for plotIspVOF (overrideable)
        self.OF_low               = 1.0
        self.OF_high              = 4.0

        # plot section, controls whether that plot is generated in generate_report()
        # overrideable via input file
        self.plot_engine_contour    = True
        self.plot_nozzle_contour    = True
        self.plot_isp_vs_of         = True
        self.plot_mach              = True
        self.plot_heat_flux         = True
        self.plot_wall_temp         = True
        self.plot_thermal_stress    = True
        self.plot_coolant_velocity  = True
        self.plot_hx_coefficient    = True
        self.plot_coolant_temp      = True
        self.plot_channel_width     = True

        # report output settings
        self.report_output_dir      = '.'   # directory where report files are saved
        self.report_dpi             = 150   # DPI for plot images embedded in report

        ################# ADD NEW MONOPROP   ########################################
        self.monoMode             = 0
        self.newMonoBool          = 0
        self.newMonoName          = ''
        self.newMono              = [['H2O2(L)', 'H 2 O 2', '100']]

        #################     ADD NEW FUEL    #######################################
        self.newFuelBool          = 0
        self.newFuelName          = '75_ETH'
        self.newFuel              = [['C2H5OH(L)', 'C 2 H 6 O 1', '75', 'Ethanol'],
                                     ['H2O(L)',    'H 2 O 1',     '25', 'Water']]

        #################     ADD NEW OXIDIZER    ###################################
        self.newOxBool            = 0
        self.newOxName            = 'GelN2O4'
        self.newOx                = [['N2O4(L)', 'N 2 O 4',  '96.5'],
                                     ['SiO2',    'Si 1 O 2', '3.5']]

        # checck if input file was provided
        if input_file is not None:
            self._load_input_file(input_file)

        # take class input variables
        self.thrust    = self.thrust_lbf * 4.44822162          # convert lbf → N
        self.Pc        = self.Pc_psi * 6894.75729              # convert psi → Pa
        self.Pa        = 101325 * (1 - 2.25577e-5 * self.height_of_optimization) ** 5.25588
        self.Pe        = self.Pa                               # assumed perfect expansion
        self.PcOvPe    = self.Pc / self.Pe

    # parse input file and accept variable values
    def _load_input_file(self, filepath):
        """
        Parse a TOML input file and apply its values to this Engine instance.

        TOML is Python-native via the built-in `tomllib` (Python 3.11+) or the
        `tomli` back-port for earlier versions.  It is a simple key=value format
        that maps directly to Python types (int, float, bool, str, list-of-lists)
        with zero custom parsing logic required.

        Supported keys
        --------------
        Required (at least these should appear in every input file):
            OF, Pc, thrust, height_of_optimization, burnTime

        Optional overrides (all Engine defaults apply when absent):
            engine_name, num_characteristics, Lstar, alpha, beta,
            coolantTempStart, Kwall, Kc, chamberWallThickness,
            channelHeight, channelWallThickness, numChannels, temp_step,
            coefThermEx, youngMod, currFuel, currOx, oxCooled,
            monoMode, newMonoBool, newMonoName, newMono,
            newFuelBool, newFuelName, newFuel,
            newOxBool, newOxName, newOx,
            OF_low, OF_high
        """
        try:
            import tomllib          # Python 3.11+
        except ImportError:
            try:
                import tomli as tomllib   # pip install tomli  (Python < 3.11)
            except ImportError:
                raise ImportError(
                    "TOML support requires Python 3.11+ (built-in tomllib) "
                    "or the 'tomli' package (pip install tomli)."
                )

        with open(filepath, 'rb') as f:
            data = tomllib.load(f)

        # whitelist attributes the input file is allowed to set
        # Pc, Pa, thrust are recomputed after this
        allowed = {
            # required mission parameters
            'OF', 'Pc_psi', 'thrust_lbf', 'height_of_optimization', 'burnTime',
            # geometry / design
            'engine_name', 'num_characteristics', 'Lstar', 'alpha', 'beta',
            # heat transfer / cooling
            'coolantTempStart', 'Kwall', 'Kc', 'chamberWallThickness',
            'channelHeight', 'channelWallThickness', 'numChannels',
            'temp_step', 'coefThermEx', 'youngMod',
            # propellant selection
            'currFuel', 'currOx', 'oxCooled',
            'newFuelBool', 'newFuelName', 'newFuel',
            'newOxBool', 'newOxName', 'newOx',
            'monoMode', 'newMonoBool', 'newMonoName', 'newMono',
            'oxTankTemp', 'fuTankTemp',
            # OF sweep
            'OF_low', 'OF_high',
            # plot selection flags
            'plot_engine_contour', 'plot_nozzle_contour', 'plot_isp_vs_of',
            'plot_mach', 'plot_heat_flux', 'plot_wall_temp', 'plot_thermal_stress',
            'plot_coolant_velocity', 'plot_hx_coefficient', 'plot_coolant_temp',
            'plot_channel_width',
            # # wall material for mat_dict lookup
            # 'wall_material',
            # report settings
            'report_output_dir', 'report_dpi',
        }

        for key, value in data.items():
            if key in allowed:
                setattr(self, key, value)
            else:
                print(f"[Engine] WARNING: unknown input key '{key}' ignored.")
                
        if self.oxCooled:
            self.coolantTempStart = self.oxTankTemp
        else:
            self.coolantTempStart = self.fuTankTemp

    def getCeaProperties(self):
        
        if self.newMonoBool:
            
            # pull the weight percentages out of the new propellant data and sum to check that it meets 100%        
            # throw error if components do not add to 100% weight
            if np.sum(np.array([row[2] for row in self.newMono], dtype=float)) != 100:
                raise ValueError('NEW PROPELLANT COMPONENTS DO NOT ADD TO 100%')

            # build new propellant card string
            self.newMonoCard = ' '.join(
                'name ' + row[0] + ' ' + row[1] + ' wt%=' + row[2]
                for row in self.newMono
            )


        if self.newFuelBool:
            
            # pull the weight percentages out of the new propellant data and sum to check that it meets 100%
            # throw error if components do not add to 100% weight
            if np.sum(np.array([row[2] for row in self.newFuel], dtype=float)) != 100:
                raise ValueError('NEW FUEL COMPONENTS DO NOT ADD TO 100%')

            # build new propellant card string
            self.newFuelCard = ' '.join(
                'fuel ' + row[0] + ' ' + row[1] + ' wt%=' + row[2]
                for row in self.newFuel
            )
            

        if self.newOxBool:
    
            # pull the weight percentages out of the new propellant data and sum to check that it meets 100%
            # throw error if components do not add to 100% weight
            if np.sum(np.array([row[2] for row in self.newOx], dtype=float)) != 100:
                raise ValueError('NEW OXIDIZER COMPONENTS DO NOT ADD TO 100%')

            # build new propellant card string
            self.newOxCard = ' '.join(
                'oxid ' + row[0] + ' ' + row[1] + ' wt%=' + row[2]
                for row in self.newOx
            )
        
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
        else:
            # translate currOx from any naming convention to the CEA name
            # e.g. "Oxygen" → "O2", "LOX" → "O2", "N2O4" → "N2O4(L)"
            self.currOx = translate(self.currOx, 'cea')

        # if new fuel was created, use that for CEA
        if self.newFuelBool:
            add_new_fuel(self.newFuelName, self.newFuelCard)
            self.currFuel=self.newFuelName
        else:
            # translate currFuel from any naming convention to the CEA name
            # e.g. "Ethanol" → "C2H5OH(L)", "Methane" → "CH4(L)"
            self.currFuel = translate(self.currFuel, 'cea')

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
        
        # create an array of several expansion ratios to analytically see how changes effect performance
        eps = [self.eps-5, self.eps, self.eps+5]
        for e in eps:
            
            # create an array of mixture ratios between levels specified in user input
            # these will be combined with the expansion ratio with the isps plotted for various points
            mrArr = np.arange(self.OF_low, self.OF_high, 0.05)  
            ispArr = [self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, MR, e)[0] for MR in mrArr]

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
        pi = np.pi

        # Calculations
        self.To = self.Tc # stagnation temp K assumed to be adiabatic flame temp in chamber
        self.Pt = self.Pc*(1+((self.gamma-1)/2))**(-self.gamma/(self.gamma-1)) # pressure at throat Pa
        self.Tt = self.Tc*(1/(1+((self.gamma-1)/2))) # Throat temp K
        self.Me = np.sqrt((2/(self.gamma-1))*(((self.Pc/self.Pa)**((self.gamma-1)/self.gamma))-1)) # Mach at exit
        self.Po = self.Pt*(1/(1+(((self.gamma-1)/2)))**(-self.gamma/(self.gamma-1))) # stagnation pressure Pa using Pt
        self.Ueq = np.sqrt(((2*self.gamma)/(self.gamma-1))*self.Rbar*self.To*(1-((self.Pe/self.Po)**((self.gamma-1)/self.gamma)))) # exit velocity m/s
        self.mdot = self.thrust/self.Ueq # kg

        self.mdotFuel = self.mdot / (self.OF + 1)
        self.mdotOx = self.mdot - self.mdotFuel

        bigGamma = np.sqrt(self.gamma*((2/(self.gamma+1))**((self.gamma+1)/(self.gamma-1))))
        self.At = (self.mdot*np.sqrt(self.Rbar*self.To))/(self.Po*bigGamma) # area of the throat m2
        self.Ae = (self.At/self.Me)*((1+((self.gamma-1)/2)*(self.Me**2))/((self.gamma+1)/2))**((self.gamma+1)/(2*(self.gamma-1))) # area of the exit m2

        self.Isp = self.Ueq/9.81

        self.cp = (self.gamma/(self.gamma-1))*self.Rbar

        self.Dt = np.sqrt((4/pi)*self.At)
        self.De = np.sqrt((4/pi)*self.Ae)
        self.Dc = self.Dt*3


        self.Ac = (np.pi * self.Dc**2)/4
        self.Vc = self.Lstar*self.At
        self.Lc = self.Vc / (1.1 * self.Ac)

        # solving for the chamber diameter iteratively
        self.theta = np.tan(np.radians(self.beta))
        for _ in np.arange(10):
            self.Dc = np.sqrt((self.Dt**3 + (24/np.pi)*self.theta*self.Vc) / (self.Dc + 6*self.theta*self.Lc))

        # nozzle geometry curves
        self.Re = self.De / 2
        self.Rt = self.Dt / 2
        self.Rc3 = 0.382 * self.Rt
        self.Rc2 = 1.5 * self.Rt
        self.Rc1 = self.Rc2 * 2
        self.Rc = self.Dc / 2

        # distances from throat (negative)
        tan = np.tan(np.radians(self.alpha))
        # sin = np.sin(np.radians(self.alpha))
        cos = np.cos(np.radians(self.alpha))
        self.throat = 0
        self.endStraight = -1*((self.Rc2 - (self.Rc2*cos)) / np.tan(np.radians(self.alpha/2)))*1000
        conv_vert_portion = self.Rc - self.Rt - (self.Rc2 - (self.Rc2*cos)) - (self.Rc1 - (self.Rc1*cos))
        self.endRc1 = self.endStraight - (conv_vert_portion / tan)*1000
        self.endChamber = self.endRc1 - ((self.Rc1 - (self.Rc1*cos)) / np.tan(np.radians(self.alpha/2)))*1000
        self.startChamber = self.endChamber - self.Lc*1000
        self.endRc3 = 1000*(self.Rc3 * np.cos(np.radians(90-self.beta))) / (np.tan(np.radians(self.beta)))
           
    def propFlowRates(self):
        
        self.getCeaProperties()
        self.engineDimensions()
        
        if self.monoMode == 1:
            raise Exception('PROPELLANT RATES NOT AVAILABLE IN MONOPROP MODE\nSET ''monoMode = 0'' AND USE BI PROPELLANT MODE')
        
        if self.newFuelBool:
            if len(self.newFuel) == 2:
            
                # still needed for CoolProp
                fuel_comp_2 = float(self.newFuel[1][2]) / 100
                fuel_card_comp_1 = self.newFuel[0][1]
                fuel_card_comp_1.replace(' ', '')
                fuel_card_comp_1 = Formula(fuel_card_comp_1)
                fuel_comp_1_MM = fuel_card_comp_1.mass
                fuel_card_comp_2 = self.newFuel[1][1]
                fuel_card_comp_2.replace(' ', '')
                fuel_card_comp_2 = Formula(fuel_card_comp_2)
                fuel_comp_2_MM = fuel_card_comp_2.mass
        
                # find values for CoolProp in heat transfer calcs
                fuel_comp_1_frac = round(((1-fuel_comp_2) * fuel_comp_1_MM) / 
                                        (((1-fuel_comp_2) * fuel_comp_1_MM) + (fuel_comp_2 * fuel_comp_2_MM)),4)
                fuel_comp_2_frac = round(1 - fuel_comp_1_frac, 4)
                # remove the 0 before the decimal for CoolProp
                fuel_comp_1_frac = str(fuel_comp_1_frac)[1:]
                fuel_comp_2_frac = str(fuel_comp_2_frac)[1:]
        
                self.fuelComposition = self.newFuel[0][3]+'['+fuel_comp_1_frac+']&'+self.newFuel[1][3]+'['+fuel_comp_2_frac+']'
                
            else:
                
                self.fuelComposition = self.newFuel[0][3]
                
        else:
            self.fuelComposition = self.currFuel
        
        if self.newOxBool:
            if len(self.newOx) == 2:
            
                # still needed for CoolProp
                ox_comp_2 = float(self.newOx[1][2]) / 100
                ox_card_comp_1 = self.newOx[0][1]
                ox_card_comp_1.replace(' ', '')
                ox_card_comp_1 = Formula(ox_card_comp_1)
                ox_comp_1_MM = ox_card_comp_1.mass
                ox_card_comp_2 = self.newOx[1][1]
                ox_card_comp_2.replace(' ', '')
                ox_card_comp_2 = Formula(ox_card_comp_2)
                ox_comp_2_MM = ox_card_comp_2.mass
        
                # find values for CoolProp in heat transfer calcs
                ox_comp_1_frac = round(((1-ox_comp_2) * ox_comp_1_MM) / 
                                        (((1-ox_comp_2) * ox_comp_1_MM) + (ox_comp_2 * ox_comp_2_MM)),4)
                ox_comp_2_frac = round(1 - ox_comp_1_frac, 4)
                # remove the 0 before the decimal for CoolProp
                ox_comp_1_frac = str(ox_comp_1_frac)[1:]
                ox_comp_2_frac = str(ox_comp_2_frac)[1:]
        
                self.oxComposition = self.newOx[0][3]+'['+ox_comp_1_frac+']&'+self.newOx[1][3]+'['+ox_comp_2_frac+']'
                
            else:
                
                self.oxComposition = self.newOx[0][3]
                
        else:
            # NOTE This is rather poorly designed
            # translate currOx to a fluid name map
            # e.g. "O2" (CEA) → "Oxygen" (CoolProp), "LOX" → "Oxygen"
            try:
                self.oxComposition = translate(self.currOx, 'coolprop')
            except ValueError as e:
                print(e)
            try:
                self.oxComposition = translate(self.currOx, 'cea')
            except ValueError as e:
                print(e)
            try:
                self.oxComposition = translate(self.currOx, 'rocketprops')
            except ValueError as e:
                print(e)
            try:
                self.oxComposition = translate(self.currOx, 'fluid_dict')
            except ValueError as e:
                print(e)
                
                

        # density and volume/sec calcs from tanks
        # use get_props() which tries fluid_dict → rocketprops → CoolProp in order
        if _HAS_FLUID_NAMES:
            _ox_rho_props  = _get_fluid_props(self.oxComposition,  self.oxTankTemp, 101325,
                                              ['density_kg_m3'])
            _fu_rho_props  = _get_fluid_props(self.fuelComposition, self.fuTankTemp, 101325,
                                              ['density_kg_m3'])
            self.rhoOX   = _ox_rho_props['density_kg_m3']   # kg/m3
            self.rhoFuel = _fu_rho_props['density_kg_m3']   # kg/m3
            
        else:
            
            self.rhoOX   = PropsSI('D', 'T', self.oxTankTemp,  'P', 101325, self.oxComposition) # kg/m3
            self.rhoFuel = PropsSI('D', 'T', self.fuTankTemp, 'P', 101325, self.fuelComposition) # kg/m3

        self.VdotOX = (self.mdotOx / self.rhoOX) * 1000 # L/s
        self.VdotFuel = (self.mdotFuel / self.rhoFuel) * 1000 # L/s

        self.volumeFuel = self.VdotFuel * self.burnTime # L
        self.volumeOX = self.VdotOX * self.burnTime # L
        self.massFuel = self.volumeFuel * self.rhoFuel / 1000 # kg
        self.massOX = self.volumeOX * self.rhoOX / 1000 # kg
        self.massProp = self.massFuel + self.massOX
        
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
        print("ox Vdot (L/s): ",round(self.VdotOX,4))
        print("fuel volume (L): ",round(self.volumeFuel,2))
        print("ox volume (L): ",round(self.volumeOX,2))
        print("fuel mass (kg): ",round(self.massFuel,2))
        print("ox mass (kg): ",round(self.massOX,2))
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
        
        """Given the geometry of the engine this will return the profile of the engine
        as well as plot the profile. Nozzle geometry is found using Method of Characteristics

        Returns:
            list[float]: returns three lists that ccan be imported into CAD
            to create splines of the engine contour
        """

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
        # thetaN = np.radians(40) # degrees to rad
        # thetaE = np.radians(20) # degrees to rad

        # full x array up to the throat
        engineX_conv = np.round(
            np.arange(self.startChamber, 0.0 + 0.005, 0.01), 2
        )

        # Pre-compute straight-section slope & intercept
        x1 = self.endRc1
        y1 = np.sqrt(abs(((3*Rt)**2) - (self.endRc1 - self.endChamber)**2)) + (Rc - 3 * Rt)
        x2 = self.endStraight
        y2 = (Rt + 1.5 * Rt) - np.sqrt(abs(((1.5*Rt)**2) - (self.endStraight)**2))
        m  = (y1 - y2) / (x1 - x2)
        b  = y1 - m * x1

        # CHANGE: define a scalar function for the piecewise radius, then apply via
        # np.vectorize so each x-value is evaluated without a Python-level for-loop.
        # Math inside is unchanged — identical to the original if/elif branches.
        def _radius_at(x):
            # chamber — constant radius (cylindrical)
            if x <= self.endChamber:
                return Rc
            # first curve in conv section — double-radius curve
            if x <= self.endRc1:
                return np.sqrt(abs(((3*Rt)**2) - (x - self.endChamber)**2)) + (Rc - 3 * Rt)
            # straight conv section — line connecting the two curves
            if x <= self.endStraight:
                return m * x + b
            # curve into throat — 1.5*Rt radius circle
            return (Rt + 1.5 * Rt) - np.sqrt(abs(((1.5*Rt)**2) - x**2))

        # apply all contour math across the length of the chamber and compute y points
        engineY_conv = np.vectorize(_radius_at)(engineX_conv)

        # get nozzle geometry from MOC_nozzle
        nozX_arr, nozY_arr = MOC.nozzle(self.Me, self.num_characteristics, Rt, Re)

        # merge all cells into final lists of the engine profile
        self.engineX = np.concatenate([engineX_conv, nozX_arr]).tolist()
        self.engineY = np.concatenate([engineY_conv, nozY_arr]).tolist()

        self.nozX        = nozX_arr
        self.nozY        = nozY_arr
        self.nozXContour = nozX_arr.tolist()
        self.nozYContour = nozY_arr.tolist()

        # create z values of 0 for 3d spline importing
        self.nozZContour = np.zeros(len(self.nozXContour)).tolist()
        
    def saveNozGeo(self, filePath):
        
        """Given the geometry of the nozzle this will return the profile of the engine
        as well as plot the profile. Nozzle geometry is found using Method of Characteristics

        Returns:
            list[float]: returns three lists that ccan be imported into CAD
            to create splines of the nozzle contour
        """
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
        """Takes a file path and saves the chamber geometry to that path.

        Args:
            filePath (str): path to where the geometry will be saved
        """
    # saves nozzle geometry in a txt file that can be imported into CAD 
    # software to model the engine
        
        self.engineContour()

        j = True
        k=1
        while j == True:
            filename = filePath + '/chamber_geo' + str(k) + '.txt'
            j = path.exists(filename)
            k = k + 1
            
        # create z values of 0 for 3d spline importing
        self.engineZ = np.zeros_like(self.engineY).tolist()
        df = pd.DataFrame(data={"xVals":self.engineX, "yVals":self.engineY, "zVals":self.engineZ})
        df.to_csv(filename, sep=' ', index = False, header=False)
        
    def heatTransfer(self):
        """Iteratively solves for wall and coolant temperatures at each step along
        the chamber.

        Raises:
            ValueError: INPUT VALUE FOR coolandTempStart THIS IS THE TEMPERATURE (K) OF THE FUEL IN THE TANK
            ValueError: INPUT VALUE FOR Kwall THIS IS THE THERMAL CONDUCTIVITY (W/m*K) OF THE CHAMBER WALL
            ValueError: INPUT VALUE FOR Kc THIS IS THE THERMAL CONDUCTIVITY (W/m*K) OF THE FUEL
            ValueError: INPUT VALUE FOR chamberWallThickness THIS IS THE THICKNESS (mm) OF THE CHAMBER WALL
            ValueError: INPUT VALUE FOR channelHeight THIS IS THE DEPTH (mm) OF THE COOLING CHANNELS
            ValueError: INPUT VALUE FOR channelWallThickness THIS IS THE WALL THICKNESS (mm) BETWEEN COOLING CHANNELS
            ValueError: INPUT VALUE FOR numChannels THIS IS THE NUMBER OF COOLING CHANNELS AROUND THE CHAMBER
            Exception: REGEN COOLING NOT AVAILABLE IN MONOPROP MODE, SET ''monoMode = 0'' AND USE BI PROPELLANT MODE
        """
        
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
        # if self.oxCooled and self.oxComposition != "O2":
        #     raise ValueError('OX REGEN COOLING ONLY SUPPORTED WITH OXYGEN ("O2" MUST BE USED FOR CEA)')

        self.engineContour()
        
        if self.monoMode == 1:
            raise Exception('REGEN COOLING NOT AVAILABLE IN MONOPROP MODE\nSET ''monoMode = 0'' AND USE BI PROPELLANT MODE')
        
        self.propFlowRates()
        
        # #
        # #
        # #
        # # FIX
        # #
        # #
        # #

        # # ── mat_dict integration ──────────────────────────────────────────────
        # # If wall_material names a key in mat_dict, Kwall / coefThermEx /
        # # youngMod are looked up per-station from the dictionary instead of
        # # being held at the scalar TOML values.
        # #
        # # Loading strategy (lazy / low-memory):
        # #   mat_dict is imported once as a Python module.  Python's module cache
        # #   (sys.modules) means the file is read from disk only on the very first
        # #   call; every subsequent call in the same session reuses the already-
        # #   loaded dict at zero I/O cost.  Only the three properties actually
        # #   consumed by heatTransfer() are extracted from each entry — the rest
        # #   of the dict is never touched.
        # #
        # # Temperature-dependence strategy:
        # #   Wall material properties (especially Kwall and youngMod) change
        # #   significantly with temperature.  Rather than picking one value for
        # #   the whole engine, we call get_mat_props(material, Tw1[i]) at each
        # #   station so that the hot-gas side uses the property at the local wall
        # #   temperature.  get_mat_props() linearly interpolates between the two
        # #   nearest registered temperature points, so it is fast and requires no
        # #   additional CoolProp calls.
        # #
        # # Fallback:
        # #   If wall_material is '' or mat_dict cannot be found, the scalar
        # #   values from self.Kwall / self.coefThermEx / self.youngMod are used
        # #   unchanged — identical to the original behaviour.

        # _use_mat_dict = False
        # _mat_get = None   # will be get_mat_props function if available

        # if self.wall_material:
        #     try:
        #         import os as _os, sys as _sys, importlib.util as _ilu, importlib as _il
        #         _mat_path = _os.path.join(
        #             _os.path.dirname(_os.path.abspath(__file__)), 'mat_dict.py')
        #         if 'mat_dict' not in _sys.modules and _os.path.exists(_mat_path):
        #             _spec = _ilu.spec_from_file_location('mat_dict', _mat_path)
        #             _md_mod = _il.util.module_from_spec(_spec)
        #             _spec.loader.exec_module(_md_mod)
        #             _sys.modules['mat_dict'] = _md_mod
        #         if 'mat_dict' in _sys.modules:
        #             _mat_get = _sys.modules['mat_dict'].get_mat_props
        #             # Verify the material key actually exists before committing
        #             if _mat_get(self.wall_material, 298.15) is not None:
        #                 _use_mat_dict = True
        #                 print(f"[Engine] mat_dict: using '{self.wall_material}' "
        #                       f"for Kwall, coefThermEx, youngMod")
        #             else:
        #                 print(f"[Engine] WARNING: wall_material '{self.wall_material}' "
        #                       f"not found in mat_dict — using TOML scalar values.")
        #     except Exception as _e:
        #         print(f"[Engine] WARNING: mat_dict load failed ({_e}). "
        #               f"Using TOML scalar values for wall properties.")
        
        # change from mm to m for heat transfer math
        X = np.array(self.engineX) / 1000  
        Y = np.array(self.engineY) / 1000  
        # flip arrays around to be able to iterate in a single counter flow direction
        Y = Y[::-1]
        X = X[::-1]
        
        # heat transfer stuff
        muE = 6.5974511*(10**-5) # mu at throat kg/ms cp of air at Te of 2124 K
        w = 0.6 # some constant for the HX eqn I really don't know
        r = 0.9 # some other constant fr dunno
        
        self.Te = self.Tc / (1 + ((self.gamma-1)/2)*self.Me**2) # exit temp (works for freestream temp too)

        # create arrays of 0s for the math to be completed on
        n = len(X)

        self.oTo          = np.zeros(n)           
        self.hg           = np.ones(n)            
        self.qw           = np.zeros(n)           
        self.qw_safer     = np.zeros(n)           
        self.M            = np.zeros(n)           
        self.T            = np.zeros(n)           
        self.Tr           = np.zeros(n)           
        self.Tw1          = np.full(n, 350.0)     
        self.Tw2          = np.zeros(n)           
        self.thermalStress= np.zeros(n)           
        # TCoolant has n+1 elements (one per station plus the outlet)
        self.TCoolant     = np.full(n + 1, float(self.coolantTempStart))  
        self.hc           = np.ones(n)            
        self.coolV        = np.zeros(n)           
        self.channelWidth = np.zeros(n)           
        
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
        
        tempCp  = 0
        tempMu  = 0
        tempRho = 0
        tempArray = np.round(np.arange(self.coolantTempStart, 500, self.temp_step),2)

        # select coolant for regen cooling
        if self.oxCooled:
            coolant_fluid = self.oxComposition
        else:
            coolant_fluid = self.fuelComposition

        # ── Unified fluid property lookup ─────────────────────────────────────
        # get_props() from fluid_names.py tries sources in priority order:
        #   1. fluid_dict  — in-memory cache, zero I/O after first call
        #   2. rocketprops — fast, covers LOX, RP-1, MMH, N2H4, N2O4, etc.
        #   3. CoolProp    — broadest coverage including mixtures; slowest
        #
        # Name translation is handled inside get_props() so coolant_fluid can
        # be any format (CoolProp, CEA, RocketProps, canonical) — it is
        # resolved to the correct string for whichever source answers first.
        #
        # Results are written back into fluid_dict on a cache miss so
        # subsequent calls within the same session are fully served from memory.
        _P = float(self.Pc)
        _props_needed = ['density_kg_m3', 'dynamic_viscosity_Pa_s', 'specific_heat_J_kgK']
        rho_list, mu_list, Cp_list = [], [], []

        # Load fluid_dict for write-back on cache misses
        try:
            import os as _os2, sys as _sys2, importlib.util as _ilu2, importlib as _il2
            _fd_path2 = _os2.path.join(
                _os2.path.dirname(_os2.path.abspath(__file__)), 'fluid_dict.py')
            if 'fluid_dict' not in _sys2.modules and _os2.path.exists(_fd_path2):
                _spec2 = _ilu2.spec_from_file_location('fluid_dict', _fd_path2)
                _fd_mod2 = _il2.util.module_from_spec(_spec2)
                _spec2.loader.exec_module(_fd_mod2)
                _sys2.modules['fluid_dict'] = _fd_mod2
            _add_state = _sys2.modules['fluid_dict'].add_fluid_state
        except Exception:
            _add_state = None

        for _T in tempArray:
            _T = float(_T)
            if _HAS_FLUID_NAMES:
                # try:
                _p = _get_fluid_props(coolant_fluid, _T, _P, _props_needed)
                rho_list.append(_p['density_kg_m3'])
                mu_list.append(_p['dynamic_viscosity_Pa_s'])
                Cp_list.append(_p['specific_heat_J_kgK'])
                # write back to fluid_dict if the answer came from rocketprops or coolprop
                if _p.get('_source') != 'fluid_dict' and _add_state is not None:
                    _fd_k = fluid_dict_key(coolant_fluid) if _HAS_FLUID_NAMES else coolant_fluid
                    try:
                        from fluid_names import fluid_dict_key
                        _fd_k = fluid_dict_key(coolant_fluid)
                    except Exception:
                        _fd_k = coolant_fluid
                    try:
                        _k_val = float(PropsSI('L', 'T', _T, 'P', _P,
                                        translate(coolant_fluid, 'coolprop')))
                    except Exception:
                        _k_val = self.Kc
                    _add_state(_fd_k, round(_T, 4), round(_P, 2), {
                        'temperature_K':             _T,
                        'pressure_Pa':               _P,
                        'density_kg_m3':             _p['density_kg_m3'],
                        'dynamic_viscosity_Pa_s':    _p['dynamic_viscosity_Pa_s'],
                        'specific_heat_J_kgK':       _p['specific_heat_J_kgK'],
                        'thermal_conductivity_W_mK': _k_val,
                        'prandtl_number': (
                            _p['dynamic_viscosity_Pa_s'] * _p['specific_heat_J_kgK']
                            / _k_val if _k_val else 0),
                        'notes': f'Auto-populated via {_p.get("_source","?")}'
                                    f' @ T={_T}K P={_P}Pa',
                    })
                # except Exception:
                #     # last-resort direct CoolProp call
                #     cp_str = translate(coolant_fluid, 'coolprop')
                #     rho_list.append(float(PropsSI('D', 'T', _T, 'P', _P, cp_str)))
                #     mu_list.append(float(PropsSI('V', 'T', _T, 'P', _P, cp_str)))
                #     Cp_list.append(float(PropsSI('C', 'T', _T, 'P', _P, cp_str)))
            else:
                # fluid_names not available — direct CoolProp (original behaviour)
                rho_list.append(float(PropsSI('D', 'T', _T, 'P', _P, coolant_fluid)))
                mu_list.append(float(PropsSI('V', 'T', _T, 'P', _P, coolant_fluid)))
                Cp_list.append(float(PropsSI('C', 'T', _T, 'P', _P, coolant_fluid)))
            rho = np.array(rho_list)
            mu  = np.array(mu_list)
            Cp  = np.array(Cp_list)
        else:
            rho = PropsSI('D', 'T', tempArray, 'P', self.Pc, coolant_fluid)
            mu  = PropsSI('V', 'T', tempArray, 'P', self.Pc, coolant_fluid)
            Cp  = PropsSI('C', 'T', tempArray, 'P', self.Pc, coolant_fluid)

        # map the temperature arrays to the specific heat, density, and viscousity
        self.CpDict  = dict(zip(tempArray, Cp))   
        self.muDict  = dict(zip(tempArray, mu))   
        self.rhoDict = dict(zip(tempArray, rho))  

        for i in range(len(X)):
            
            # set flag for when the calculations pass nozzle to allow for
            # number of channels to change along chamber
            if X[i] == 0:
                pastThroat = True
                
            if i > 0:
                self.Tw1[i] = max(self.Tw1[i - 1]-2, float(self.TCoolant[i]))
                

            while (abs(self.hg[i]*(self.Tr[i]-self.Tw1[i])) - (self.hc[i]*(self.Tw2[i] - self.TCoolant[i])) > 1):

                # #
                # #
                # #
                # # FIX
                # #
                # #
                # #
                # # ── Per-station wall material properties ─────────────────────────
                # # If mat_dict is active, interpolate Kwall / coefThermEx / youngMod
                # # at the current wall temperature (Tw1[i]).  On the very first
                # # iteration Tw1[i] is the physics-based initial guess; as the inner
                # # while-loop converges it climbs toward the true wall temperature,
                # # but the material properties are only re-fetched at the *station*
                # # level (once per station), not inside the inner convergence loop.
                # # This keeps the property lookup cost O(stations), not
                # # O(stations × convergence_iterations).
                # #
                # # If mat_dict is not active, _Kwall / _CTE / _E stay equal to the
                # # TOML scalars set during __init__ — no change to the equations.
                # if _use_mat_dict:
                #     _mat_props = _mat_get(self.wall_material, float(self.Tw1[i]))
                #     if _mat_props is not None:
                #         _Kwall   = _mat_props.get('thermal_conductivity_W_mK', self.Kwall)
                #         _CTE     = _mat_props.get('coef_thermal_expansion_1_K', self.coefThermEx)
                #         _E_GPa   = _mat_props.get('youngs_modulus_GPa', self.youngMod)
                #     else:
                #         # Tw1[i] may be outside the registered range — use scalars
                #         _Kwall, _CTE, _E_GPa = self.Kwall, self.coefThermEx, self.youngMod
                # else:
                #     _Kwall, _CTE, _E_GPa = self.Kwall, self.coefThermEx, self.youngMod
                    
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
                
                self.T[i]   = (self.Tc / (1+ ((self.gamma-1)/2) * self.M[i]**2))
                self.Tr[i]  = (self.T[i]*(1+((self.gamma-1)/2)*r*(self.M[i]**2)))
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
                self.qw_safer[i] = (self.qw[i]*1.1)
                
                # temp of wall on the coolant side in order 
                # for the combustion side to stay at temp we want
                self.Tw2[i] = (self.Tw1[i] - (self.qw_safer[i]/(self.Kwall/(self.chamberWallThickness/1000))))
                
                tempCp  = self.CpDict[round(self.TCoolant[i],0)]
                tempMu  = self.muDict[round(self.TCoolant[i],0)]
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

        self.qw           = self.qw[::-1]
        self.qw_safer     = self.qw_safer[::-1]
        self.hc           = self.hc[::-1]
        self.Tw1          = self.Tw1[::-1]
        self.Tw2          = self.Tw2[::-1]
        self.thermalStress= self.thermalStress[::-1]
        self.channelWidth = self.channelWidth[::-1]
        # change bacck to mm
        self.channelWidth = self.channelWidth* 1000
        self.TCoolant     = self.TCoolant[::-1]
        self.coolV        = self.coolV[::-1]
        self.M            = self.M[::-1]
        # shorten the TCoolant array by 1
        self.TCoolant = self.TCoolant[:-1]
        
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

    def generate_report(self):
        """
        Run all analyses, save selected plots as PNG files, write a JSON data
        file, then call generate_report.js (Node.js / docx-js) to build the
        final Word document (.docx).

        Report sections
        ---------------
        1. INPUTS        – every user-settable attribute as valid TOML, copy-
                        pasteable to create a new input file.
        2. CHAMBER INFO  – all outputs equivalent to printChamberInfo().
        3. PLOTS         – only the plots whose flag is True, two per page,
                        sized to fit US Letter (9 × 4 inches each).

        Files written to self.report_output_dir:
            report_data.json              consumed by generate_report.js
            plot_<name>.png              one per selected plot
            <engine_name>_report.docx    final Word document
        """
        import os, json, subprocess
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        out_dir = self.report_output_dir
        os.makedirs(out_dir, exist_ok=True)
        dpi = self.report_dpi

        # FIG_H 4.0 → 3.125 (= 300 px ÷ 96 dpi) to match IMG_H_PX = 300 in JS.
        # The PNG physical size must equal the docx image slot so Word renders
        # at native resolution with no scaling.  3.125 in × 2 plots = 6.25 in,
        # leaving 0.25 in of headroom inside the 6.5 in landscape content area.
        FIG_W, FIG_H = 9.0, 3.125
        saved_plots = []   # [{title, path}, …] consumed by JS builder
        
        def _save(flag, title, plot_fn):
            if not flag:
                return
            fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=dpi)
            plot_fn(fig, ax)
            slug  = title.lower().replace(' ', '_').replace('/', '_')
            fname = os.path.join(out_dir, f'plot_{slug}.png')
            fig.tight_layout()
            fig.savefig(fname, dpi=dpi)
            plt.close(fig)
            saved_plots.append({'title': title, 'path': os.path.abspath(fname)})


        # run analyses
        self.engineContour()
        
        _save(self.plot_engine_contour, 'Engine Contour', lambda fig, ax: (
            ax.plot(self.engineX, self.engineY),
            ax.set_aspect(1.0),
            ax.grid(color='k', alpha=0.25),
            ax.set_title('Engine Geometry'),
            ax.set_xlabel('Distance Along Engine Axially (mm)'),
            ax.set_ylabel('Radius (mm)'),
        ))

        _save(self.plot_nozzle_contour, 'Nozzle Contour', lambda fig, ax: (
            ax.plot(self.nozX, self.nozY, 'b'),
            ax.set_aspect(1.0),
            ax.grid(color='k', alpha=0.25),
            ax.set_title('Nozzle Geometry (MOC)'),
            ax.set_xlabel('Distance Along Engine Axially (mm)'),
            ax.set_ylabel('Radius (mm)'),
        ))

        if not self.monoMode and not self.newOxBool:
            self.propFlowRates()
            self.heatTransfer()
            if self.plot_isp_vs_of:
                fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=dpi)
                for e in [self.eps - 5, self.eps, self.eps + 5]:
                    mr_arr  = np.arange(self.OF_low, self.OF_high, 0.05)
                    isp_arr = [self.rocketcea_eng_obj.get_IvacCstrTc(self.Pc_psi, MR, e)[0]
                                for MR in mr_arr]
                    ax.plot(mr_arr, isp_arr, label=f'AreaRatio {e:.0f}')
                ax.set_title('Isp vs Mixture Ratio')
                ax.set_xlabel('Mixture Ratio (O/F)')
                ax.set_ylabel('Isp (s)')
                ax.legend()
                ax.grid(alpha=0.25)
                fig.tight_layout()
                fname = os.path.join(out_dir, 'plot_isp_vs_of.png')
                fig.savefig(fname, dpi=dpi)
                plt.close(fig)
                saved_plots.append({'title': 'Isp vs Mixture Ratio',
                                    'path': os.path.abspath(fname)})

            _save(self.plot_mach, 'Mach Number', lambda fig, ax: (
                ax.plot(self.engineX, self.M),
                ax.grid(alpha=0.25),
                ax.set_title('Mach Number Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Mach Number'),
            ))
            _save(self.plot_heat_flux, 'Heat Flux', lambda fig, ax: (
                ax.plot(self.engineX, self.qw,       label='qw'),
                ax.plot(self.engineX, self.qw_safer, label='qw safer'),
                ax.legend(), ax.grid(alpha=0.25),
                ax.set_title('Heat Flux Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Heat Flux (W/m²)'),
            ))
            _save(self.plot_wall_temp, 'Wall Temperature', lambda fig, ax: (
                ax.plot(self.engineX, self.Tw1, label='Hot Side'),
                ax.plot(self.engineX, self.Tw2, label='Coolant Side'),
                ax.legend(), ax.grid(alpha=0.25),
                ax.set_title('Wall Temperature Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Temperature (K)'),
            ))
            _save(self.plot_thermal_stress, 'Thermal Stress', lambda fig, ax: (
                ax.plot(self.engineX, self.thermalStress),
                ax.grid(alpha=0.25),
                ax.set_title('Thermal Stress Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Thermal Stress (kPa)'),
            ))
            _save(self.plot_coolant_velocity, 'Coolant Velocity', lambda fig, ax: (
                ax.plot(self.engineX, self.coolV),
                ax.grid(alpha=0.25),
                ax.set_title('Coolant Velocity Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Velocity (m/s)'),
            ))
            _save(self.plot_hx_coefficient, 'HX Coefficient', lambda fig, ax: (
                ax.plot(self.engineX, self.hc),
                ax.grid(alpha=0.25),
                ax.set_title('Heat Transfer Coefficient (Cooling Passages)'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('h_c  W/(m²·K)'),
            ))
            _save(self.plot_coolant_temp, 'Coolant Temperature', lambda fig, ax: (
                ax.plot(self.engineX, self.TCoolant),
                ax.grid(alpha=0.25),
                ax.set_title('Coolant Temperature Along Passages'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Temperature (K)'),
            ))
            _save(self.plot_channel_width, 'Channel Width', lambda fig, ax: (
                ax.plot(self.engineX, self.channelWidth),
                ax.grid(alpha=0.25),
                ax.set_title('Cooling Channel Width Along Engine'),
                ax.set_xlabel('Distance Along Engine Axially (mm)'),
                ax.set_ylabel('Width (mm)'),
            ))

        # build toml compatible inputs
        inputs = {
            'OF':                    self.OF,
            'Pc_psi':                self.Pc_psi,
            'thrust_lbf':            self.thrust_lbf,
            'height_of_optimization': self.height_of_optimization,
            'burnTime':              self.burnTime,
            'engine_name':           self.engine_name,
            'num_characteristics':   self.num_characteristics,
            'Lstar':                 self.Lstar,
            'alpha':                 self.alpha,
            'beta':                  self.beta,
            'OF_low':                self.OF_low,
            'OF_high':               self.OF_high,
            'currFuel':              self.currFuel,
            'currOx':                self.currOx,
            'oxCooled':              bool(self.oxCooled),
            'coolantTempStart':      self.coolantTempStart,
            'Kwall':                 self.Kwall,
            'Kc':                    self.Kc,
            'chamberWallThickness':  self.chamberWallThickness,
            'channelHeight':         self.channelHeight,
            'channelWallThickness':  self.channelWallThickness,
            'numChannels':           self.numChannels,
            'temp_step':             self.temp_step,
            'coefThermEx':           self.coefThermEx,
            'youngMod':              self.youngMod,
            'monoMode':              self.monoMode,
            'newMonoBool':           self.newMonoBool,
            'newMonoName':           self.newMonoName,
            'newMono':               self.newMono,
            'newFuelBool':           self.newFuelBool,
            'newFuelName':           self.newFuelName,
            'newFuel':               self.newFuel,
            'newOxBool':             self.newOxBool,
            'newOxName':             self.newOxName,
            'newOx':                 self.newOx,
            'plot_engine_contour':   bool(self.plot_engine_contour),
            'plot_nozzle_contour':   bool(self.plot_nozzle_contour),
            'plot_isp_vs_of':        bool(self.plot_isp_vs_of),
            'plot_mach':             bool(self.plot_mach),
            'plot_heat_flux':        bool(self.plot_heat_flux),
            'plot_wall_temp':        bool(self.plot_wall_temp),
            'plot_thermal_stress':   bool(self.plot_thermal_stress),
            'plot_coolant_velocity': bool(self.plot_coolant_velocity),
            'plot_hx_coefficient':   bool(self.plot_hx_coefficient),
            'plot_coolant_temp':     bool(self.plot_coolant_temp),
            'plot_channel_width':    bool(self.plot_channel_width),
            'report_output_dir':     self.report_output_dir,
            'report_dpi':            self.report_dpi,
        }

        # build chamber info dict
        chamber = {
            'Ae_over_At':    round(self.Ae / self.At, 3),
            'Pc_kPa':        round(self.Pc / 1000, 2),
            'Tc_K':          round(self.Tc, 2),
            'thrust_N':      round(self.thrust, 2),
            'Pt_kPa':        round(self.Pt / 1000, 2),
            'Tt_K':          round(self.Tt, 2),
            'Me':            round(self.Me, 4),
            'cea_Me':        round(self.cea_Me, 4),
            'Po_kPa':        round(self.Po / 1000, 2),
            'To_K':          round(self.To, 2),
            'Ueq_m_s':       round(self.Ueq, 2),
            'mdot_kg_s':     round(self.mdot, 3),
            'Isp_s':         round(self.Isp, 2),
            'cea_Isp_s':     round(self.cea_Isp, 2),
            'Cstar_m_s':     round(self.Cstar, 2),
            'cp_J_kgK':      round(self.cp, 2),
            'Vc_m3':         round(self.Vc, 6),
            'Lc_m':          round(self.Lc, 6),
            'Dc_m':          round(self.Dc, 6),
            'At_m2':         round(self.At, 6),
            'Dt_m':          round(self.Dt, 6),
            'Ae_m2':         round(self.Ae, 6),
            'De_m':          round(self.De, 6),
            'Rc_mm':         round(self.Rc  * 1000, 4),
            'Rc1_mm':        round(self.Rc1 * 1000, 4),
            'Rc2_mm':        round(self.Rc2 * 1000, 4),
            'Rt_mm':         round(self.Rt  * 1000, 4),
            'Rc3_mm':        round(self.Rc3 * 1000, 4),
            'Re_mm':         round(self.Re  * 1000, 4),
            'mdotFuel_kg_s':  round(self.mdotFuel, 3),
            'mdotOx_kg_s':    round(self.mdotOx, 3),
        }
        if not self.monoMode and not self.newOxBool:
            chamber.update({
                'VdotFuel_L_s':   round(self.VdotFuel, 4),
                'VdotOx_L_s':     round(self.VdotOX, 4),
                'volFuel_L':      round(self.volumeFuel, 2),
                'volOx_L':        round(self.volumeOX, 2),
                'massFuel_kg':    round(self.massFuel, 2),
                'massOx_kg':      round(self.massOX, 2),
                'massProp_kg':    round(self.massProp, 2),
            })

        # json for JS builder
        report_name = (self.engine_name.strip() or 'engine') + '_report'
        data = {
            'engine_name': self.engine_name or 'Engine Report',
            'report_name': report_name,
            'out_dir':     os.path.abspath(out_dir),
            'inputs':      inputs,
            'chamber':     chamber,
            'plots':       saved_plots,
        }
        json_path = os.path.join(out_dir, 'report_data.json')
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
            

        # JS docx builder
        # script_dir = Path(__file__).parent
        script_dir = os.path.dirname(os.path.abspath(__file__))
        js_script  = os.path.join(script_dir, 'generate_report.js')
        # check that the path exists
        if not os.path.exists(js_script):
            print(f"[Engine] WARNING: generate_report.js not found at {js_script}.")
            print(f"[Engine] JSON data written to: {json_path}")
            return        
        
        # check for node.js installed
        node_exe = shutil.which('node')
        if node_exe is None:
            raise EnvironmentError(
                "Node.js not found. Install it from https://nodejs.org or add it to PATH."
            )
        result = subprocess.run(['node', js_script, json_path],
                                capture_output=True, text=True)
        if result.returncode != 0:
            print('[Engine] ERROR from generate_report.js:\n', result.stderr)
        else:
            print(f"[Engine] Report: {os.path.join(out_dir, report_name + '.docx')}")

    @classmethod
    def from_file(cls, filepath):
        """
        Create an Engine object entirely from a TOML file.

        Usage
        -----
            eng = Engine.from_file('my_engine.toml')

        The five required parameters (OF, Pc_psi, thrust_lbf,
        height_of_optimization, burnTime) must be present in the file.
        All other parameters fall back to class defaults if omitted.
        """
        return cls(input_file=filepath)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python Engine_Class.py <input_file.toml> [output_dir]")
        sys.exit(1)

    engine = Engine.from_file(sys.argv[1])
    if len(sys.argv) > 2:
        engine.report_output_dir = sys.argv[2]
    engine.generate_report()
