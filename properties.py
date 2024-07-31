# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:33:45 2023

@author: wangz
"""

import pandas as pd
import math
import numpy as np

from dash.exceptions import PreventUpdate
import streamlit as st

class ReservoirProperties:
    
    """Entry for Reservoir Properties

        Attributes:
           P0, k, temp, Thickness, Porosity, k
           Rock cr, Salinity

        Properties units:
            xxx
            
        Also includes validation check. 
        """
    def __init__(self, P0, temp, Thickness, salinity, Porosity, k, cr, XL, YL, BXL, BYL, FracP, BHPext):
        
        self.P0 = P0 
        self.temp = temp  # in degC
        self.Thickness = Thickness  
        self.salinity = salinity  
        self.Porosity = Porosity  # unitless
        self.k = k  # in mDarcy
        self.cr = cr
        self.XL= XL  
        self.YL = YL
        self.BXL = BXL
        self.BYL = BYL
        self.FracP = FracP
        self.BHPext = BHPext
        
        
        
class SimProperties:
    def __init__(self, nWXmin, nWXmax, nWE, rW, rE, SimTime, CP, nWX_npv_max):
        self.nWXmin = nWXmin
        self.nWXmax = nWXmax
        self.nWE = nWE
        self.rW = rW
        self.rE = rE
        #self.InjRate = InjRate
        #self.ExtRate = ExtRate
        self.SimTime = SimTime 
        self.CP = CP
        self.nWX_npv_max = nWX_npv_max
        
class RelakProperties:
    def __init__(self, kra0, krg0, m, n, Sar, Sgc):
        self.kra0 = kra0
        self.krg0 = krg0
        self.m = m
        self.n = n  
        self.Sar =  Sar 
        self.Sgc = Sgc
        

class NPVProperties:
    def __init__(self, initinv, opcostup, opcostdwn, txcr, disrate):
        self.initinv = initinv
        self.opcostup = opcostup
        self.opcostdwn = opcostdwn
        self.txcr = txcr
        self.disrate = disrate



def read_rsvr_input(units):

    rsvrprop = ReservoirProperties(    
       units.P0.values[0], units.temp.values[0], units.Thickness.values[0], units.salinity.values[0], units.Porosity.values[0], 
       units.k.values[0], units.cr.values[0], math.sqrt(units.Area_res.values[0]), math.sqrt(units.Area_res.values[0]), 
       math.sqrt(units.Area_Basin.values[0]), math.sqrt(units.Area_Basin.values[0]), units.FracP.values[0], 
       units.BHPext.values[0])
    
    simprop = SimProperties(
       units.nWXmin.values[0].astype(int), units.nWXmax.values[0].astype(int), units.nWE.values[0].astype(int), 
       units.rW.values[0], math.sqrt(units.Area_Basin.values[0]/math.pi),
       units.SimTime.values[0], units.CP.values[0], units.nWXmax.values[0].astype(int))

    
    relakprop = RelakProperties(
       units.kra0.values[0], units.krg0.values[0], units.m.values[0], units.n.values[0], units.Sar.values[0], units.Sgc.values[0])

    
    npvprop = NPVProperties(
       units.initinv.values[0], units.opcostup.values[0], units.opcostdwn.values[0], units.txcr.values[0], units.disrate.values[0])

    return rsvrprop, simprop, relakprop, npvprop

def Unit_conversion(input_df,UnitSys): 
    '''
    # The coefficient conversion for the following properties, converted from SI unit or field unit
    
    nWXmin,	nWXmax,	nWE, Area_res, Area_Basin, rW, P0, InjRate,
	ExtRate,cr	,Thickness,	k	,Porosity,	temp, salinity	
    CP , SimTime, kra0	, krg0, m ,n , Sar, Sgc, FracP, BHPext	
    drco, drcoExt, txcr, maco, macoExt, moco
                    
    The input xlsx locaiton is hard-coded. When the input format subject to change, please also
    update code here. 
    
    Salinity: assume Nacl as saline, 58.448 g/mol. water density 997g/L
    
    Update: Since S1 has been deactivated, several line are modified below, to remove InjRate and ExtRate inputs
    '''
    

   
    SIcoef = np.array([1, 1, 1, 1000**2, 1000**2, 1, 1, 1, 1, 9.87E-16,
              1, 1, 1, 1, 31536000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    
 
    Fieldcoef = np.array([1, 1, 1, 1609.344**2, 1609.344**2, 0.3048, 0.00689476,
                 0.000145038, 0.3048, 9.87E-16, 1, 	
                 1, 1.71607e-5, 0.00689476, 31536000, 1, 1, 1, 1, 1, 1, 0.00689476, 0.00689476, 
                 1, 1, 1, 1, 1, 1])


    output_df = input_df.copy()
    input_array = input_df.iloc[1,:].astype(float)

    SIcoef_array = np.array(SIcoef, dtype=float)
    Fieldcoef_array = np.array(Fieldcoef, dtype=float)

    old_cols_names = ['P0', 'temp', 'Thickness', 'k', 'salinity', 'm', 'n', 'cr', 'FracP', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc',
                      'nWXmin','nWXmax','nWE','Area_res','Area_Basin','rW','SimTime','BHPext',
                      'initinv','opcostup','opcostdwn','txcr','disrate','CP','simmodule']
  
    new_cols_names = ['Pressure', 'Temperature', 'Thickness', 'Permeability', 'Salinity', 'Water Corey Exponent (m)', 'Gas Corey Exponent (n)', 'Rock Compressibility', 'Maximum Allowable Injection Pressure', 'Endpoint Water Relative Permeability', 
                      'Endpoint Gas Relative Permeability', 'Porosity', 'Residual Water Saturation', 'Critical Gas Saturation',
                      '# of Injection wells (min)','# of Injection wells (max)','# of Extraction wells', 'Project Area','Reservoir Area','Injection Well Radius','Injection Duration', 'Min Extraction Pressure',
                      'Initial Investment','Annual Operational Cost (Upstream)','Annual Operational Cost (Downstream)','Tax Credit', 'Discount Rate','Critical Pressure Increase','Simulation Module']

    units_SI = [ 'MPa','C','m','mD', 'mol/Kg', '-','-','1/Pa','MPa','-','-','-','-','-','-',
                '-','-','km^2','km^2','m','year','MPa','$M','$M','$M/well','$/ton','-','MPa','-']
    units_field = [ 'psi','F','ft','mD', 'mg/L', '-','-','1/psi','psi','-','-','-','-','-','-',
                '-','-','mile^2','mile^2','ft','year','psi','$M','$M','$M/well','$/ton','-','psi','-']


    if UnitSys == 2: # if SI unit
        output_array = SIcoef_array * input_array  # multiply coefficient
        output_df.iloc[1,:] = output_array.copy() 
        output_df = output_df.rename(columns=output_df.iloc[0]).drop(output_df.index[0])
        # drop the unnecessary row, and rename the columns
       
        output_df['nWXmin'] = output_df['nWXmin'].astype(int)
        output_df['nWXmax'] = output_df['nWXmax'].astype(int)
        output_df['nWE'] = output_df['nWE'].astype(int)
        
        # then, to prepare the displayed input with units
        output_print = output_df.copy()
        output_print.iloc[0, :] = input_array
        output_print = output_print.reindex(columns=old_cols_names)
        output_print['nWXmin'] = output_df['nWXmin'] **2
        output_print['nWXmax'] = output_df['nWXmax'] **2
        for old_name, new_name in zip(old_cols_names, new_cols_names):
            output_print.rename(columns={old_name: new_name}, inplace=True)
        
        units_row = pd.DataFrame([units_SI], columns=output_print.columns)
        #column_dtypes = {col: output_print[col].dtype for col in units_row.columns}
        #units_row = units_row.astype(column_dtypes)
        output_print = pd.concat([output_print, units_row.rename(index=dict(zip(units_row.index, old_cols_names)))], axis=0, ignore_index=True)
        
        return output_df, output_print
    
    elif UnitSys == 1: # if Field Unit
       output_array = Fieldcoef_array * input_array # multiply coefficient
       output_df.iloc[1,:] = output_array.copy()
       output_df = output_df.rename(columns=output_df.iloc[0]).drop(output_df.index[0])
       # drop the unnecessary row, and rename the columns
       
       output_df['nWXmin'] = output_df['nWXmin'].astype(int)
       output_df['temp']  = (output_df['temp']-32)*5/9 # special treatment for temperature, from f to c
       output_df['nWXmax'] = output_df['nWXmax'].astype(int)
       output_df['nWE'] = output_df['nWE'].astype(int)
       
       # then, to prepare the displayed input with units
       output_print = output_df.copy()
       output_print.iloc[0, :] = input_array
       output_print = output_print.reindex(columns=old_cols_names)
       output_print['nWXmin'] = output_df['nWXmin'] **2
       output_print['nWXmax'] = output_df['nWXmax'] **2
#       output_print = output_df.reindex(old_cols_names, axis=1)
       for old_name, new_name in zip(old_cols_names, new_cols_names):
           output_print.rename(columns={old_name: new_name}, inplace=True)
       
       units_row = pd.DataFrame([units_field], columns=output_print.columns)
       output_print = pd.concat([output_print, units_row.rename(index=dict(zip(units_row.index, old_cols_names)))], axis=0, ignore_index=True)

       return output_df, output_print
    
    else:
        raise Exception("Wrong Unit System selection!")


def check_parameter_ranges(output_df, UnitSys):

    #Errors for input values out of range
    output_df = output_df.astype(float)
    error_tracker = 'False'
      
    if output_df.iloc[0, 0] < 1 or output_df.iloc[0, 0] > 10:
        st.warning('Number of injection wells should be between less than 10x10.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 1] < 1 or output_df.iloc[0, 1] > 10:
        st.warning('Number of injection wells should be between less than 10x10.', icon="⚠️")
        error_tracker = 'True'
    #if output_df.iloc[0, 2] < 0 or output_df.iloc[0, 2] > 16:
    #    st.warning('Number of extraction wells should be between less than 4x4.', icon="⚠️")
    #    error_tracker = 'True'
    if output_df.iloc[0, 1]**2 <= output_df.iloc[0, 2]:
        st.warning('Maximun number of injectors must be bigger than number of extractors.',  icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 0]**2 < output_df.iloc[0, 2]:
        st.warning('Minimum number of injectors should not be less than number of extractors. General geometry option gives you more choices to control number of injection and extraction wells.',  icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 3] > output_df.iloc[0, 4]:
        st.warning('Reservoir Area can not be smaller than Project Area.', icon="⚠️")
        error_tracker = 'True'

    valid_numbers = [0, 4, 8]
    if output_df.iloc[0, 2] not in valid_numbers:
        st.warning('Number of extractors must be one of the following numbers: 0, 4, 8.', icon="⚠️")
        error_tracker = 'True'
    
    if output_df.iloc[0, 5] < 0.05 or output_df.iloc[0, 5] > 0.5:
        if UnitSys == 2:
            st.warning('Wellbore radius should be >0.05 m and <0.5 m.', icon="⚠️")
        elif UnitSys == 1:
            st.warning('Wellbore radius should be >0.164 ft and <1.64 ft.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 6] < 7 or output_df.iloc[0, 6] > 55:
        if UnitSys == 2:
            st.warning('Initial reservoir pressure should be >7 MPa and <55 MPa.', icon="⚠️")
        elif UnitSys == 1:
            st.warning('Initial reservoir pressure should be >1015 psi and <7977 psi.', icon="⚠️")
        error_tracker = 'True'
    
    #if output_df.iloc[0, 7] < -0.01 or output_df.iloc[0, 7] > 1000:
    #    if UnitSys == 2:
    #        st.warning('Injection rate <1000 ton/day.', icon="⚠️")
    #    elif UnitSys == 1:
    #        st.warning('Injection rate <0.365 MMT/yr.', icon="⚠️")
    #    error_tracker = 'True'
    #if output_df.iloc[0, 8] < -0.01 or output_df.iloc[0, 8] > 1000:
    #    if UnitSys == 2:
    #        st.warning('Extraction rate should be <1000m^3/day.', icon="⚠️")
    #    elif UnitSys == 1:
    #        st.warning('Extraction rate should be <6290 bbl/day.', icon="⚠️")
    #    error_tracker = 'True'
    
    if output_df.iloc[0, 7] < -0.01 or output_df.iloc[0, 7] > 1e-3:
        if UnitSys == 2:
            st.warning('Compressibility should be <1E-3 1/Pa.', icon="⚠️")
        elif UnitSys == 1:
            st.warning('Compressibility should be <6.895 1/Psi.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 8] < -0.01 or output_df.iloc[0, 8] > 500:
        if UnitSys == 2:
            st.warning('Reservoir thickness should be <501 m. If you have thick reservoir, try to break it down '
                       'into multiple sections.', icon="⚠️")
        elif UnitSys == 1:
            st.warning('Reservoir thickness should be <1640 ft. If you have thick reservoir, try to break it down '
                       'into multiple sections.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 9] < -0.01 or output_df.iloc[0, 9] > 1e-11:
        st.warning('Permeability value is out of range.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 10] < 0 or output_df.iloc[0, 10] > 1:
        st.warning('Porosity value is out of range.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 11] < 31 or output_df.iloc[0, 11] > 107:
        st.warning('Temperature value is out of range.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 12] < -0.01 or output_df.iloc[0, 12] > 4:
        st.warning('Salinity value is out of range.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 13] < -0.01 or output_df.iloc[0, 13] > 50:
        if UnitSys == 2:
            st.warning('Critical pressure should be <50 MPa.', icon="⚠️")
        elif UnitSys == 1:
            st.warning('Critical pressure should be <7252 psi.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 14] < -0.01 or output_df.iloc[0, 14] > 3.2e9:
        st.warning('Injection Duration is too large.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 15] < 0.1 or output_df.iloc[0, 15] > 1:
        st.warning('Kra0 should be <1 and >0.1.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 16] < 0.1 or output_df.iloc[0, 16] > 1:
        st.warning('Krg0 should be <1 and >0.1.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 17] < 1 or output_df.iloc[0, 17] > 6:
        st.warning('m should be <6 and >1.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 18] < 1 or output_df.iloc[0, 18] > 6:
        st.warning('n should be <6 and >1.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 19] < 0.05 or output_df.iloc[0, 19] > 0.9:
        st.warning('Residual water saturation should be <0.9 and >0.05.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0, 20] < 0.05 or output_df.iloc[0, 20] > 0.9:
        st.warning('Residual gas saturation should be <0.9 and >0.05.', icon="⚠️")
        error_tracker = 'True'
    
    if output_df.iloc[0, 21] < output_df.iloc[0, 6] or output_df.iloc[0, 21] > 2.0 *output_df.iloc[0, 6]:
        st.warning('Maximum allowable injection pressure should be >initial reservoir pressure but <2.0 x initial reservoir pressure.', icon="⚠️")
        error_tracker = 'True'
    if output_df.iloc[0,28] == 2 and output_df.iloc[0, 2] >0:
        if output_df.iloc[0, 22] < 7 or output_df.iloc[0, 22] > output_df.iloc[0, 6] :
            if UnitSys == 2:
                st.warning('Extraction well BHP should be < initial reservoir pressure and > 7 MPa', icon="⚠️")
            elif UnitSys == 1:
                st.warning('Extraction well BHP should be < initial reservoir pressure and > 1015 psi', icon="⚠️")
            error_tracker = 'True'

    if error_tracker == 'True':
        st.stop()
        
def userexcel_conversion(userinput_df,userinput_inj,userinput_ext):
    input_df_raw = pd.read_excel('ET_input.xlsx',sheet_name='Parameters',header=None)
    input_df = input_df_raw.copy()
    input_df.iloc[1,0] = round(userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Number of Injection wells (min)'].index[0], 1]**0.5)
    input_df.iloc[1,1] = round(userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Number of Injection wells (max)'].index[0], 1]**0.5)
    input_df.iloc[1,2] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Number of Extraction wells'].index[0], 1]
    input_df.iloc[1,3] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Project Area'].index[0], 1]
    input_df.iloc[1,4] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Reservoir Area'].index[0], 1]
    input_df.iloc[1,5] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Injection Well Radius'].index[0], 1]
    input_df.iloc[1,6] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Initial Pressure'].index[0], 1]
    #input_df.iloc[1,7] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Injection Rate'].index[0], 1]
    #input_df.iloc[1,8] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Extraction Rate'].index[0], 1]
    input_df.iloc[1,7] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Rock Compressibility'].index[0], 1]
    input_df.iloc[1,8] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Thickness'].index[0], 1]
    input_df.iloc[1,9] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Permeability'].index[0], 1]
    input_df.iloc[1,10] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Porosity'].index[0], 1]
    input_df.iloc[1,11] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Temperature'].index[0], 1]
    input_df.iloc[1,12] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Salinity'].index[0], 1]
    input_df.iloc[1,13] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Critical Pressure Increase (AOR calculation)'].index[0], 1]
    input_df.iloc[1,14] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Injection Duration'].index[0], 1]
    input_df.iloc[1,15] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Endpoint Water Relative Permeability (kra0)'].index[0], 1]
    input_df.iloc[1,16] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Endpoint Gas Relative Permeability (krg0)'].index[0], 1]
    input_df.iloc[1,17] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Water Corey Exponent (m)'].index[0], 1]
    input_df.iloc[1,18] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Gas Corey Exponent (n)'].index[0], 1]
    input_df.iloc[1,19] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Residual Water Saturation (Sar)'].index[0], 1]
    input_df.iloc[1,20] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Critical Gas Saturation (Sgc)'].index[0], 1]
    input_df.iloc[1,21] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Max Allowable Injection Pressure'].index[0], 1]
    input_df.iloc[1,22] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Min Extraction Pressure'].index[0], 1]
    input_df.iloc[1,23] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Initial Investment'].index[0], 1]
    input_df.iloc[1,24] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Annual Operational Cost (Upstream)'].index[0], 1]
    input_df.iloc[1,25] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Annual Operational Cost (Downstream)'].index[0], 1]
    input_df.iloc[1,26] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Tax Credit'].index[0], 1]
    input_df.iloc[1,27] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Discount Rate'].index[0], 1]
    #input_df.iloc[1,30] = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Monitoring Cost'].index[0], 1]
    proj_name = userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Project Name'].index[0], 1]
    if userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Select a simulation module'].index[0], 1] == 'User Given Inputs':
        input_df.iloc[1, 28] = 1
        input_df.iloc[1,21] = userinput_inj.iloc[1:, 4].max()
        input_df.iloc[1,22] = userinput_ext.iloc[1:, 4].min()
    elif userinput_df.iloc[userinput_df[userinput_df.iloc[:, 0] == 'Select a simulation module'].index[0], 1] == 'Maximum Storage Capacity':
        input_df.iloc[1, 28] = 2
    else:
        # Set a default value if neither condition is met
        input_df.iloc[1, 28] = 0  # Replace 0 with any desired default value
    return input_df, proj_name


'''
from properties import Unit_conversion, read_rsvr_input
import pandas as pd
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_field.xlsx'
UnitSys = 1
BCValue = 2

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
userinput_inj = pd.read_excel(uploaded_file, sheet_name='Injectors', header=None)
userinput_ext = pd.read_excel(uploaded_file, sheet_name='Extractors', header=None)

input_df, proj_name = userexcel_conversion(userinput_df,userinput_inj,userinput_ext)
output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)
check_parameter_ranges(output_df, UnitSys)
'''