# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:38:33 2023
Scenario 2: Given Pressures
@author: wangz
"""
from Matcal import NumInj_build, Distance_Mat_Builder, A_ConstRate
import numpy as np
import math
from reservoir_funcs import Fluid_Prop 

###############################################
#TESTING

def Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi):
     
    if Sensi == 0:
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
        
        
        ## Scenario 2: When final pressures are , flow rates to be calcualted
        Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
        geometry = 0
        Fluid = Fluid_Prop(rsvrprop.FracP, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)
        Fluid.Spycher()
        mug, cg, rhog, rhoc = Fluid.CO2Prop()
        cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
        tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
        
        MaxP = (rsvrprop.FracP-rsvrprop.P0)*1e6
        MinP = (rsvrprop.BHPext-rsvrprop.P0)*1e6
        
        Qsave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
        Qsave_res = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F') # rm^3 for Inj and Ext
        Psave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
        rT1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
        rL1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
        for nWX in range(simprop.nWXmin,simprop.nWXmax + 1):  
            
            NumInj_t = int(NumInj[nWX-1])
            nWT = NumInj_t + simprop.nWE 
            
            Q = np.zeros(nWT,order = 'F')
            Qsave = np.zeros(nWT,order = 'F')
            B = np.zeros(nWT,order = 'F')
            chi_BL = np.zeros(nWT,order = 'F')
            chi_dry = np.zeros(nWT,order = 'F')
        
            Q[:nWT] = np.full(nWT,0.001)
            Qsave[:nWT] = np.full(nWT,1)
            Qinj = [0] * NumInj_t
            Qext = [0] * simprop.nWE
            
        
            err = 1
            it = 1
        
            NumInj_t = int(NumInj[nWX-1])
           
            while err > 1e-3:
        
                
                epsilon=np.multiply(Q,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
                chi_BL[:nWT]=np.full(nWT,(1/4)*epsilon*dfgdSgL)
                chi_dry[:nWT]=np.full(nWT,(1/4)*epsilon*dfgdSgT)
                for i in range(NumInj_t):
                    B[i] = 2 * math.pi * rsvrprop.Thickness * rsvrprop.k * relakprop.krg0 / mug * MaxP            
                for i in range(NumInj_t,nWT):
                    B[i] = 2 * math.pi * rsvrprop.Thickness * rsvrprop.k * relakprop.kra0 / mub * MinP
            
                A = A_ConstRate(nWX, simprop, rsvrprop, NumInj_t, simprop.rW, chi_BL, chi_dry, Lamda_g,
                                    Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue,geometry)
            
                Atrans = np.transpose(A)
                Qtrans = np.linalg.solve(Atrans, B)
                Q = np.transpose(Qtrans)
                err = np.linalg.norm(Qsave - Q, 2)
                Qsave = np.copy(Q)
                it += 1
               
                #print(err)
                if it>200:
                   break
               
                   
            for i in range(NumInj_t):
                rT1[i, nWX-1] = np.sqrt(zT * Qsave[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
                rL1[i, nWX-1] = np.sqrt(zL * Qsave[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
                Qsave_res[i, nWX-1] = Qsave[i]
                Qinj[i] = Qsave[i]*86400*rhoc/1000
                Qsave2[i, nWX-1] = np.copy(Qinj[i]) # ton per day
                Psave2[i, nWX-1] = rsvrprop.FracP
                
            for i in range(simprop.nWE):
                Qsave_res[i + NumInj_t,  nWX-1] = Qsave[i + NumInj_t]
                Qext[i] = Qsave[i + NumInj_t] * rhob / rhobstd * 86400 # m^3 per day
     #           Qsave2[i + NumInj_t, nWX-1] = -1 * Qext[i] *rhobstd/1000 # ton per day
                Qsave2[i + NumInj_t, nWX-1] = -1 * Qext[i]  # m^3 per day
                Psave2[i + NumInj_t, nWX-1] = rsvrprop.BHPext
                
                
                
                
            L = sum(Qinj)
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
        return Psave2, Qsave2, rL1, rT1, Capacity, NPV, Qsave_res
    
    if Sensi == 1 :
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
        
        nWX = simprop.nWX_npv_max
        
        NumInj = NumInj_build(simprop.nWXmin,simprop.nWXmax)
       
        
        ## Scenario 2: When final pressures are , flow rates to be calcualted
        Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
        geometry = 0
        Fluid = Fluid_Prop(rsvrprop.FracP, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)
        Fluid.Spycher()
        mug, cg, rhog, rhoc = Fluid.CO2Prop()
        cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
        tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
        
        MaxP = (rsvrprop.FracP-rsvrprop.P0)*1e6
        MinP = (rsvrprop.BHPext-rsvrprop.P0)*1e6
        
        Qsave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
        Psave2 = np.zeros((simprop.nWXmax**2 + simprop.nWE, simprop.nWXmax),order = 'F')
        rT1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
        rL1 = np.zeros((simprop.nWXmax**2+simprop.nWE, simprop.nWXmax),order = 'F')
       
       
        NumInj_t = int(NumInj[nWX-1])
        nWT = NumInj_t + simprop.nWE 
           
        Q = np.zeros(nWT,order = 'F')
        Qsave = np.zeros(nWT,order = 'F')
        B = np.zeros(nWT,order = 'F')
        chi_BL = np.zeros(nWT,order = 'F')
        chi_dry = np.zeros(nWT,order = 'F')
       
        Q[:nWT] = np.full(nWT,0.001)
        Qsave[:nWT] = np.full(nWT,1)
        Qinj = [0] * NumInj_t
        Qext = [0] * simprop.nWE
           
       
        err = 1
        it = 1
       
        NPV = 0   
        while err > 1e-3:
       
               
            epsilon=np.multiply(Q,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
            chi_BL[:nWT]=np.full(nWT,(1/4)*epsilon*dfgdSgL)
            chi_dry[:nWT]=np.full(nWT,(1/4)*epsilon*dfgdSgT)
            for i in range(NumInj_t):
                B[i] = 2 * math.pi * rsvrprop.Thickness * rsvrprop.k * relakprop.krg0 / mug * MaxP            
            for i in range(NumInj_t,nWT):
                B[i] = 2 * math.pi * rsvrprop.Thickness * rsvrprop.k * relakprop.kra0 / mub * MinP
           
            A = A_ConstRate(nWX, simprop, rsvrprop, NumInj_t, simprop.rW, chi_BL, chi_dry, Lamda_g,
                                   Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue,geometry)
           
            Atrans = np.transpose(A)
            Qtrans = np.linalg.solve(Atrans, B)
            Q = np.transpose(Qtrans)
            err = np.linalg.norm(Qsave - Q, 2)
            Qsave = np.copy(Q)
            it += 1
           
            #print(err)
            if it>200:
               break
              
                  
        for i in range(NumInj_t):
            rT1[i, nWX-1] = np.sqrt(zT * Qsave[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
            rL1[i, nWX-1] = np.sqrt(zL * Qsave[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
            Qinj[i] = Qsave[i]*86400*rhoc/1000
            Qsave2[i, nWX-1] = np.copy(Qinj[i]) # ton per day
            Psave2[i, nWX-1] = rsvrprop.FracP
               
        for i in range(simprop.nWE):
            Qext[i] = Qsave[i + NumInj_t] * rhob / rhobstd * 86400 # m^3 per day
            #Qsave2[i + NumInj_t, nWX-1] = -1 * Qext[i] *rhobstd/1000 # ton per day
            Qsave2[i + NumInj_t, nWX-1] = -1 * Qext[i]  # m^3 per day
            Psave2[i + NumInj_t, nWX-1] = rsvrprop.BHPext
               
                    
                    
        L = sum(Qinj)
        TotalInjRate = L # ton per day
        Capacity = TotalInjRate * simprop.SimTime / 86400 / 1000000 # Mton
             
        for years in range(1, int(simprop.SimTime / 365 / 24 / 60 / 60) + 1):
            cashflow = (Capacity / int(
                simprop.SimTime / 365 / 24 / 60 / 60) * 1000000 * npvprop.txcr - npvprop.opcostup * 1000000 - npvprop.opcostdwn * 1000000 * (
                                simprop.nWE + NumInj_t))
            annual_npv = cashflow / (1 + npvprop.disrate) ** years
            NPV = NPV + annual_npv
        NPV = (NPV - npvprop.initinv * 1000000) / 1000000  # Million $
        
        #NPV = (Capacity * 1000000 * npvprop.txcr - npvprop.drco * 1000000 * NumInj_t - npvprop.drcoExt * 1000000 * simprop.nWE 
        #            - npvprop.maco * 1000 * NumInj_t * simprop.SimTime / (24 * 60**2 * 365) - npvprop.macoExt * 1000 * simprop.nWE * simprop.SimTime 
        #            / (24 * 60**2 * 365) - npvprop.moco * 1000 * simprop.SimTime / (24 * 60**2 * 365) * Area / 1000000) / 1000000 # million $
        return Capacity, NPV
    else: 
        raise Exception('Wrong Sensitivity Analysis Selection!')
        




'''
########
## TESTING
from properties import Unit_conversion, read_rsvr_input, userexcel_conversion
import pandas as pd
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_SI.xlsx'
UnitSys = 2
BCValue = 2
Sensi = 0

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
input_df, proj_name = userexcel_conversion(userinput_df)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)

Psave2, Qsave2, rL1, rT1, Capacity, NPV, Qsave_res = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
print(Capacity)
'''