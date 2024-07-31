# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 11:17:55 2023
Scenario 1: Given Constant Rates
@author: wangz
"""
from Matcal import NumInj_build, Distance_Mat_Builder, A_ConstRate, B_ConstRate
import numpy as np
import math
from reservoir_funcs import Fluid_Prop 
from scipy.special import exp1



###############################################
#TESTING

def Scenario1(rsvrprop, simprop, relakprop, npvprop, BCValue):
    
    simprop.rE = math.sqrt(rsvrprop.BXL*rsvrprop.BYL/math.pi)
    
    if BCValue == 1:
        rED1 = simprop.rE / simprop.rW  # closed
        rED2 = simprop.rE / simprop.rW  # closed
    elif BCValue == 2:
        rED1 = 1e38  # open
        rED2 = 1e38  # open
    else:
        raise ValueError('Error: invalid BCValue.') # or use an appropriate exception type and error message


    Area = rsvrprop.XL * rsvrprop.YL
    
    NumInj = NumInj_build(simprop.nWXmin,simprop.nWXmax)
    
    TotalInjRate = [0] * simprop.nWXmax
    Capacity = [0] * simprop.nWXmax
    NPV= [0] * simprop.nWXmax
    
    
    ## Scenario 1: Constant Rate Test
    Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
    geometry = 0
    Qsave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
    Psave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
    rT1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
    rL1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
    
    for nWX in range(simprop.nWXmin,simprop.nWXmax + 1):    
        NumInj_t = int(NumInj[nWX-1])
        nWT = NumInj_t + simprop.nWE 
        Psave = np.zeros(nWT,order = 'F')
        chi_BL = np.zeros(NumInj_t,order = 'F')
        chi_dry = np.zeros(NumInj_t,order = 'F')
        Psave[:nWT] = np.full(nWT,rsvrprop.P0)
        P_ave = rsvrprop.P0
        err = 1
        it = 1
    
        
    
        while err > 1e-2:
    
            
            Fluid = Fluid_Prop(P_ave, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)
    
            Fluid.Spycher()
    
            mug, cg, rhog, rhoc = Fluid.CO2Prop()
    
            cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
    
            tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
            
            Qinj=simprop.InjRate/86400/rhoc*1000 # rm^3/s
            Qext=simprop.ExtRate/86400/rhob*rhobstd # rm^3/s
            epsilon=np.multiply(Qinj,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
            chi_BL[:NumInj_t]=np.full(NumInj_t,(1/4)*epsilon*dfgdSgL)
            chi_dry[:NumInj_t]=np.full(NumInj_t,(1/4)*epsilon*dfgdSgT)
            #NumInj_t = int(NumInj[nWX-1])
        
            A = A_ConstRate(nWX, simprop, rsvrprop, NumInj_t, simprop.rW, chi_BL, chi_dry, Lamda_g,
                                Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue, geometry)
            
            P, P_ave, Psave, err, it = B_ConstRate(A, simprop.nWE, NumInj_t, Qinj,Qext,mug,mub,rsvrprop.Thickness,rsvrprop.k,relakprop.kra0,relakprop.krg0,rsvrprop.P0,Psave,it)
#            print(err)
            if it>200:
                break
            
    
        for i in range(NumInj_t):
            rT1[i, nWX-1] = np.sqrt(zT * Qinj * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
            rL1[i, nWX-1] = np.sqrt(zL * Qinj * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
            Qsave2[i, nWX-1] = simprop.InjRate # ton per day
            Psave2[i, nWX-1] = np.copy(P[i])
        
        for i in range(NumInj_t, nWT):
#            Qsave2[i, nWX-1] = -1 * simprop.ExtRate * rhobstd / 1000 # ton per day
            Qsave2[i, nWX-1] = -1 * simprop.ExtRate  # m^3 per day
            Psave2[i, nWX-1] = np.copy(P[i])
            
        L = simprop.InjRate * NumInj_t
        TotalInjRate[nWX-1] = L # ton per day
        Capacity[nWX-1] = TotalInjRate[nWX-1] * simprop.SimTime / 86400 / 1000000 # Mton

        NPV[nWX - 1] = 0
        for years in range(1, int(simprop.SimTime / 365 / 24 / 60 / 60) + 1):
            cashflow = (Capacity[nWX - 1] / int(
                simprop.SimTime / 365 / 24 / 60 / 60) * 1000000 * npvprop.txcr - npvprop.opcostup * 1000000 - npvprop.opcostdwn * 1000000 * (
                                simprop.nWE + NumInj_t))
            annual_npv = cashflow / (1 + npvprop.disrate) ** years
            NPV[nWX - 1] = NPV[nWX - 1] + annual_npv
        NPV[nWX - 1] = (NPV[nWX - 1] - npvprop.initinv * 1000000) / 1000000  # Million $


        # NPV[nWX-1] = (Capacity[nWX-1] * 1000000 * npvprop.txcr - npvprop.drco * 1000000 * NumInj_t - npvprop.drcoExt * 1000000 * simprop.nWE
        #         - npvprop.maco * 1000 * NumInj_t * simprop.SimTime / (24 * 60**2 * 365) - npvprop.macoExt * 1000 * simprop.nWE * simprop.SimTime
        #         / (24 * 60**2 * 365) - npvprop.moco * 1000 * simprop.SimTime / (24 * 60**2 * 365) * Area / 1000000) / 1000000 # million $
   

    
    
    return Psave2, Qsave2, rL1, Capacity, NPV


'''
## TESTING
from properties import Unit_conversion, read_rsvr_input, userexcel_conversion
import pandas as pd
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_SI.xlsx'
UnitSys = 2
BCValue = 2

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
input_df, proj_name = userexcel_conversion(userinput_df)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)

Psave2, Qsave2, rL1, Capacity, NPV = Scenario1(rsvrprop, simprop, relakprop, npvprop, BCValue)

'''
