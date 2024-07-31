# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:38:14 2023

@author: wangz
"""
from typing import Any

from Matcal import Distance_Mat_Builder, A_Geometry, B_Geometry
import pandas as pd
import numpy as np
import math
from reservoir_funcs import Fluid_Prop 
import matplotlib.pyplot as plt
from scipy.special import exp1
from scipy.interpolate import griddata
import warnings
import statistics as sts
import streamlit as st
import sys
# Scenario 3: When given geometry and fixed flow rates

###############################################



def Scenario3(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, output_print, proj_name):    
    
    if UnitSys == 2:
        lencoef = 1000
    elif UnitSys == 1:
        lencoef = 1609.344
    else: 
        raise Exception ("Wrong Unit System selection!")
    
#    NumInj = NumInj_build(simprop.nWXmin,simprop.nWXmax)
    TotalInjRate = [0] * simprop.nWXmax
    Capacity = [0] * simprop.nWXmax
    NPV: Any= [0] * simprop.nWXmax
    
    

    df_inj = pd.read_excel(uploaded_file,sheet_name='Injectors')
    df_ext = pd.read_excel(uploaded_file,sheet_name='Extractors')
    
    WellNumber_inj = df_inj["Well Number"]
    WellNumber_ext = df_ext["Well Number"]
    
    #Xwell_UTM = df_inj["UTM_Well_X"]
    #Ywell_UTM = df_inj["UTM_Well_Y"]
    
    Xwell = df_inj["UTM_Well_X"]
    Xwell = pd.concat([df_inj["UTM_Well_X"], df_ext["UTM_Well_X"]], ignore_index=True)
    Ywell = df_inj["UTM_Well_Y"]
    Ywell = pd.concat([df_inj["UTM_Well_Y"], df_ext["UTM_Well_Y"]], ignore_index=True)      
    
    Xwell_input = Xwell.copy()
    Ywell_input = Ywell.copy()
    
    AllReservoirs_orig = pd.read_excel(uploaded_file, sheet_name='Geometry')
    # Check if the last row is not equal to the first row
    
    
    # Convert AllReservoirs values to km
    AllReservoirs_nonan = AllReservoirs_orig.dropna(axis=1, how='all')
    AllReservoirs_m = AllReservoirs_nonan.values  #in meter
    AllReservoirs = AllReservoirs_m / lencoef   # in km or mile
    # Get the number of reservoirs
    nReservoirs = AllReservoirs.shape[1] // 2
    
    # To store all reservoirs' geometry locations
    Xres = []
    Yres = []
    tolerance = 0.0001  # 0.01% tolerance
    for i in range(nReservoirs):
        x = (AllReservoirs_m[:, 2*i]) 
        y = (AllReservoirs_m[:, 2*i+1])
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        
        # Check if last coordinates are not equal to the first ones
        if abs(x[-1] - x[0]) >= tolerance * abs(x[0]) or abs(y[-1] - y[0]) >= tolerance * abs(y[0]):
            # Append first coordinates to arrays
            x = np.append(x, x[0])
            y = np.append(y, y[0])
        Xres.extend(x)
        Yres.extend(y)
    Xres = np.array(Xres)
    Yres = np.array(Yres)


                       
    try: 
        if UnitSys == 2:
            InjRate = df_inj["Injection Rate (Ton/day)"]
            ExtRate = df_ext["Extraction Rate (m^3/day)"]
            MaxAllowedP = df_inj["Max Injection Pressure (Mpa)"]
            MinAllowedP = df_ext["Min Extraction Pressure (Mpa)"]
        elif UnitSys == 1:
            InjRate = df_inj["Injection Rate (MMT/yr)"] * 2739.726027
            ExtRate = df_ext["Extraction Rate (bbl/day)"] * 0.158987
            MaxAllowedP = df_inj["Max Injection Pressure (psi)"] * 0.00689476
            MinAllowedP = df_ext["Min Extraction Pressure (psi)"] * 0.00689476
    except (KeyError, ValueError, ZeroDivisionError) as e:
        st.warning('Error in column names or data. Please check your input units. :exclamation:')
        sys.exit()
              
    simprop.rE = math.sqrt(rsvrprop.BXL*rsvrprop.BYL/math.pi)

    Area = rsvrprop.XL * rsvrprop.YL # Project Area
    
    # Shift well and reservoir locations based on the new center, from previous rectangle's center to average well loc coordiates in UTM. 
    
    combined_X = np.concatenate([Xwell.values, Xres])
    combined_Y = np.concatenate([Ywell.values, Yres])
    
    X_center = sts.mean(combined_X)
    Y_center = sts.mean(combined_Y)

    # To get the Basin boundary based on the well locations
    try:
        GeoRatio = (max(combined_X) - min(combined_X)) / (max(combined_Y) - min(combined_Y))
    except ZeroDivisionError:
        GeoRatio = 1
        st.warning('Please check your geometry tabs. Missing or incorrect well inputs :exclamation:')
        
    # Check if Xwell and Ywell have only 1 element
    if len(combined_X) == 1 or len(combined_Y) == 1:
        GeoRatio = 1
        st.warning(' Please check your geometry tabs. Missing or incorrect well/geometry inputs :exclamation:. The reservoir is set to be a square shape')
    # Make sure the shape of reservoir is not too strange. 
    if GeoRatio < 0.5: 
        GeoRatio = 0.5
    if GeoRatio > 2:
        GeoRatio = 2
    
        
    Area_Basin = rsvrprop.BXL*rsvrprop.BYL # Based on inputs, both reservoir and basin are set to be square first, then adjusted into rectangles
    
    rsvrprop.BXL = math.sqrt(GeoRatio*Area_Basin)
    rsvrprop.BYL = rsvrprop.BXL/GeoRatio
    
    # In case the basin cannot cover all the project area. 
    long_x_half = max(X_center-min(combined_X),max(combined_X)-X_center)
    long_y_half = max(Y_center-min(combined_Y),max(combined_Y)-Y_center) 
    if long_x_half>rsvrprop.BXL/2 or long_y_half>rsvrprop.BYL/2:
        rsvrprop.BXL = max(combined_X)-min(combined_X)
        rsvrprop.BYL = max(combined_Y)-min(combined_Y)
        X_center = (max(combined_X)+min(combined_X))/2
        Y_center = (max(combined_Y)+min(combined_Y))/2
        #Area_Basin_new = rsvrprop.BXL*rsvrprop.BYL/lencoef/lencoef
        st.warning('Please check your geometry data. The reservoir area input is too small, or there are outliers in your well/geometry input.\
                   The reservoir area has been adjusted for simulation.')
        st.warning('The Sensitivity Analysis is based on the adjusted reservoir area shown in Output Figures.')
    
    shift_x = rsvrprop.BXL/2 - X_center
    shift_y = rsvrprop.BYL/2 - Y_center


    Xwell_shift = Xwell + shift_x
    Ywell_shift = Ywell + shift_y
    
    Xwell = Xwell_shift.copy()
    Ywell = Ywell_shift.copy()
    
   
    nWT = len(WellNumber_inj)+len(WellNumber_ext)
    nWI = len(WellNumber_inj)
    
    
    
    simprop.nWE = len(WellNumber_ext)
    nWE = simprop.nWE
    
    
    
    Psave = np.zeros(nWT,order = 'F')
    chi_BL = np.zeros(nWI,order = 'F')
    chi_dry = np.zeros(nWI,order = 'F')
    Qinj = np.zeros(nWI,order = 'F')
    Qext = np.zeros(simprop.nWE,order = 'F')
    Rwell = np.zeros((nWT, nWT),order = 'F')
    rT1 = np.zeros((nWT, 1),order = 'F')
    rL1 = np.zeros((nWT, 1),order = 'F')
    Qsave2 = np.zeros((nWT, 1),order = 'F')
    Psave2 = np.zeros((nWT, 1),order = 'F')
    Qsave2_output = np.zeros((nWT, 1),order = 'F')
    Psave2_output = np.zeros((nWT, 1),order = 'F')
    
    simprop.nWXmin=1;
    simprop.nWXmax=1;
    
    # Beginning of Distance matrix
    
    if BCValue == 1:
        rED1 = simprop.rE / simprop.rW  # closed
        rED2 = simprop.rE / simprop.rW  # closed
    elif BCValue == 2:
        rED1 = 1e38  # open
        rED2 = 1e38  # open
    else:
        raise ValueError('Error: invalid BCValue.') # or use an appropriate exception type and error message
    
    
    #for i in range(nWI + simprop.nWE):
    #    if Xwell[i] < 0 or Xwell[i] > rsvrprop.BXL or Ywell[i] < 0 or Ywell[i] > rsvrprop.BYL:
    #        st.warning('Error: Wells should be located inside the reservoir.')
    #        sys.exit()
    
    
    for i in range (nWT):
        X=Xwell[i]
        Y=Ywell[i]
        for j in range (nWT):
            Rwell[i,j]= math.sqrt((X-Xwell[j])**2+(Y-Ywell[j])**2)
    
    # End of Distance matrix 
    
    P_ave = rsvrprop.P0
    err = 1
    it = 1
    nWX = 1
    
    while err > 1e-2:
    
        
        Fluid = Fluid_Prop(P_ave, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)
    
        Fluid.Spycher()
    
        mug, cg, rhog, rhoc = Fluid.CO2Prop()
    
        cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
    
        tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
        
        Qinj[:nWI]=InjRate[:nWI]/86400/rhoc*1000 # rm^3/s
        Qext[:simprop.nWE]=ExtRate[:simprop.nWE]/86400/rhob*rhobstd # rm^3/s
        epsilon=np.multiply(Qinj,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
        chi_BL[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgL)
        chi_dry[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgT)
        
    
        A = A_Geometry(nWX, simprop.nWE, nWI, simprop.rW, chi_BL, chi_dry, Lamda_g,
                            Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue, Rwell)
        
        P, P_ave, Psave, err, it = B_Geometry(A, simprop.nWE, nWI, Qinj,Qext,mug,mub,rsvrprop.Thickness,rsvrprop.k,relakprop.kra0,relakprop.krg0,rsvrprop.P0,Psave,it)
#        print(err)
        if it>200:
            break
        
    
    for i in range(nWI):
        rT1[i, nWX-1] = np.sqrt(zT * Qinj[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
        rL1[i, nWX-1] = np.sqrt(zL * Qinj[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
        Qsave2[i, nWX-1] = InjRate[i] # ton per day
        Psave2[i, nWX-1] = P[i]
        
    for i in range(nWI, nWT):
#        Qsave2[i, nWX-1] = -1 * ExtRate[i-nWI] * rhobstd / 1000 # ton per day
        Qsave2[i, nWX-1] = -1 * ExtRate[i-nWI]  # m^3 per day
        Psave2[i, nWX-1] = P[i]
        
    L = 0    
    for i in range(nWI):
        L += InjRate[i]
        
    TotalInjRate[nWX-1] = L # ton per day
    Capacity[nWX-1] = TotalInjRate[nWX-1] * simprop.SimTime / 86400 / 1000000 # Mton

    NPV[nWX - 1] = 0
    for years in range(1, int(simprop.SimTime / 365 / 24 / 60 / 60) + 1):
        cashflow = (Capacity[nWX - 1] / int(
            simprop.SimTime / 365 / 24 / 60 / 60) * 1000000 * npvprop.txcr - npvprop.opcostup*1000000 - npvprop.opcostdwn * 1000000*(
                                simprop.nWE + nWI))
        annual_npv = cashflow / (1 + npvprop.disrate) ** years
        NPV[nWX - 1] = NPV[nWX - 1] + annual_npv
    NPV[nWX - 1] = (NPV[nWX - 1] - npvprop.initinv*1000000) / 1000000  # Million $


    # NPV[nWX-1] = (Capacity[nWX-1] * 1000000 * npvprop.txcr - npvprop.drco * 1000000 * nWI - npvprop.drcoExt * 1000000 * simprop.nWE
    #         - npvprop.maco * 1000 * nWI * simprop.SimTime / (24 * 60**2 * 365) - npvprop.macoExt * 1000 * simprop.nWE * simprop.SimTime
    #         / (24 * 60**2 * 365) - npvprop.moco * 1000 * simprop.SimTime / (24 * 60**2 * 365) * Area / 1000000) / 1000000 # million $
    


    nWX = 50
    nWY = 50
    rW = simprop.rW
    nNodes = nWX*nWY   
    
    Pnode = np.zeros(nNodes,order = 'F')
    #XContour = np.zeros(nNodes + nWI + nWE, 'F')
    #YContour = np.zeros(nNodes + nWI + nWE, 'F')
    #PContour = np.zeros(nNodes + nWI + nWE, 'F')
    chi_BL = np.zeros(nWI,order = 'F')
    chi_dry = np.zeros(nWI,order = 'F')
    
    X = []
    Y = []
    
    for i in range(nWY):
        for j in range(nWX):
            X.append(rsvrprop.BXL / nWX / 2 + j * rsvrprop.BXL / nWX)
            Y.append(rsvrprop.BYL / nWX / 2 + i * rsvrprop.BYL / nWY)
            
    Rwell_comp = np.zeros((nNodes, nWI + nWE),'F')
    A = np.zeros((nNodes,nWI + nWE),'F')
    

    for i in range(nNodes):
        for j in range(nWI + nWE):
            Rwell_comp[i, j] = np.sqrt((X[i] - Xwell[j]) ** 2 + (Y[i] - Ywell[j]) ** 2)
    Rwell = np.real(Rwell_comp) # In case complex number occurs
            
    Qinj = np.zeros(nWI,order = 'F')
    Qext = np.zeros(simprop.nWE,order = 'F')
    Sa = np.zeros(nWI)
    
    Qinj[:nWI]=InjRate[:nWI]/86400/rhoc*1000 # rm^3/s
    Qext[:simprop.nWE]=ExtRate[:simprop.nWE]/86400/rhob*rhobstd # rm^3/s

    epsilon=np.multiply(Qinj,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
    chi_BL[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgL)
    chi_dry[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgT)
    
    # Injectors to Injectors
    
    for i in range(nNodes):
        for j in range(nWI):
            rD = Rwell[i,j]/rW
            
            #Sa = np.zeros(nWI)
            Sa[j] = 0.5*np.log(chi_dry[j]) + (0.5/F_Lg)*np.log(chi_BL[j]/chi_dry[j]) - 0.5*(Lamda_g/Lamda_w)*np.log(chi_BL[j]/eta_D3) + 0.5*0.577215*(1-Lamda_g/Lamda_w)
            
            if tD <= rD ** 2 / 4 / np.mean(chi_BL):
                A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(
                    rD ** 2 / 4 / tD / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                    rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                    -rED2 ** 2 / 4 / tD / eta_D3)
            elif tD >= rD ** 2 / 4 / np.mean(chi_dry):
                A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD) - 0.5 * exp1(chi_dry[j]) + 1 / 2 / F_Lg * exp1(
                    chi_dry[j] / F_Lg / eta_D2) - 1 / 2 / F_Lg * exp1(
                    chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                    chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                    rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                    -rED2 ** 2 / 4 / tD / eta_D3)
            elif tD < rD ** 2 / 4 / np.mean(chi_dry) and tD > rD ** 2 / 4 / np.mean(chi_BL):
                A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD / F_Lg / eta_D2) - 0.5 / F_Lg * exp1(
                    chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                    chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                    rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                    -rED2 ** 2 / 4 / tD / eta_D3)

  
        for j in range(nWI, nWI + nWE):
            rDE = Rwell[i, j] / rW
            if BCValue == 2:
                A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
            elif BCValue == 1:
                A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4
                
    real_A = np.real(A)
    real_A_nonzero = real_A[real_A.imag != 0]  
    
    def Pressure_contour_map():
    
        #1. Pressure Contour map
        if real_A_nonzero.size != 0:
            raise ValueError("Matrix A contains complex number. Exception raised.")
        else:
            A = np.copy(real_A)
            
        for i in range(nNodes):
            BInj = 0
            for j in range(nWI):
                BInj += A[i,j] * Qinj[j]
            BExt = 0
            for j in range(nWI, nWI + nWE):
                BExt += A[i,j] * Qext[j - nWI]
            Pnode[i] = mug / (2 * math.pi * rsvrprop.Thickness * rsvrprop.k * 
                              relakprop.krg0) * BInj / 1000000 - mub / (2 * math.pi * rsvrprop.Thickness 
                            * rsvrprop.k * relakprop.kra0) * BExt / 1000000 + rsvrprop.P0 # Krg0 to Krs
    
        # Prepare contour data for nodes
        XContour = [x-shift_x  for x in X]
        YContour = [y-shift_y  for y in Y]
        PContour = Pnode.tolist()
        PContour_output = PContour.copy()

        # Prepare contour data for wells
        for i in range(nWI + nWE):
            XContour.append(Xwell_input[i])
            YContour.append(Ywell_input[i])
            PContour.append(Psave2[i].tolist())
        
            
        #Interpolate data onto a grid
        #x1 = np.arange(0, rsvrprop.BXL/lencoef, 1/200*rsvrprop.BXL/lencoef)
        #x2 = np.arange(0, rsvrprop.BYL/lencoef, 1/200*rsvrprop.BYL/lencoef)
        x1 = np.arange(X_center-rsvrprop.BXL/2, X_center+rsvrprop.BXL/2, 1/200*rsvrprop.BXL)
        x2 = np.arange(Y_center-rsvrprop.BYL/2, Y_center+rsvrprop.BYL/2, 1/200*rsvrprop.BYL)
        
        X1, X2 = np.meshgrid(x1, x2)
        PContour_flatten = flatten_list(PContour) # Make PContour in a list format
        highlight_AOR = rsvrprop.P0 + simprop.CP  # Specify the contour level to highlight
        

            
        # UNIT CONVERSION
        if UnitSys == 1: # field unit, then output pressure is in psi
            PContour_flatten_unit = [x * 145.03768 for x in PContour_flatten] # Convert from MPa to psi
            Psave2_output = Psave2 * 145.03768 # Convert from MPa to psi
            highlight_AOR_unit = highlight_AOR * 145.03768 # Convert from MPa to psi
            Qsave2_output[:nWI] = Qsave2[:nWI] * 0.000365 #Injection well convert from ton/day to MMT/year
            Qsave2_output[nWI:nWI+nWE] = Qsave2[nWI:nWI+nWE] * 6.2898107 # Extraction well convert from m^3/day to bbl/day
            
        else: # SI unit
            PContour_flatten_unit = [x * 1 for x in PContour_flatten] #Output pressure still in MPa
            Psave2_output = Psave2  # Output pressure still in MPa
            highlight_AOR_unit = highlight_AOR  # in MPa
            Qsave2_output[:nWI] = Qsave2[:nWI] * 0.000365 #Injection well convert from ton/day to MMT/year
            Qsave2_output[nWI:nWI+nWE] = Qsave2[nWI:nWI+nWE] # Extraction well in m^3/day
        
        # Capture warning messages
        #if min(PContour_flatten_unit) > highlight_AOR_unit or max(PContour_flatten_unit) < highlight_AOR_unit:
        #       warnings.warn('Critical Pressure is not suitable! AOR cannot be visualized!')
        
        foo = griddata((XContour, YContour), PContour_flatten_unit, (X1, X2), method='linear') #200 by 200 pressure matrix
        
        # Now foo is the interpolated data on the grid defined by X1 and X2
        X3 = foo  
        
        # Create a contour plot with 10 levels
        fig_contour, ax = plt.subplots()
        
        contour = ax.contourf(X1, X2, X3, 10)
        
        highlight_contour = ax.contour(X1, X2, X3, levels=[highlight_AOR_unit], colors='red', linewidths=2) # Mark AOR
        
        highlight_locations = []
        for path_collection in highlight_contour.collections:
            paths = path_collection.get_paths()
            for path in paths:
                vertices = path.vertices
                for x, y in vertices:
                    highlight_locations.append((x, y))
    
        # Set plot aspect ratio, title, and axis labels
        plt.gca().set_aspect('equal', adjustable='box')
        cbar = plt.colorbar(contour)
        
        # Add labels to the highlighted contour lines
        try:
            plt.clabel(highlight_contour, inline=True, fmt='AOR', colors='red', fontsize=10)
        except Exception as e:
            # Handle any exception that occurs
            st.warning("Warning: Cannot plot the desired AOR within the reservoir boundaries, the estimated AOR is too large (or too small).")
    
        if UnitSys == 2:
            plt.title('Pressure Contour, MPa')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
            cbar.set_label('Pressure, MPa')
        elif UnitSys ==1:
            plt.title('Pressure Contour, psi')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
            cbar.set_label('Pressure, psi')
        
            
        plt.axis('equal')
        return highlight_locations, fig_contour, Psave2_output, PContour_output
    highlight_locations, fig_contour, Psave2_output, PContour_output =  Pressure_contour_map()
    

    def Plumes_map():
        #2. Plumes map
        
        
        fig_plume, ax = plt.subplots()
        # Injectors
        for w in range(nWI):
            circle = plt.Circle((Xwell_input[w], Ywell_input[w]), radius = rL1[w], color='r', fill = False)
            ax.add_patch(circle)  # Use add_patch() instead of add_artist() for circles
        
        # Extractors   
        for w in range(nWE):
            ax.scatter(Xwell_input[w+nWI], Ywell_input[w+nWI], c='b', marker='^')
            
        # Set limits for figure
        ax.set_xlim(X_center-rsvrprop.BXL/2, X_center+rsvrprop.BXL/2)
        ax.set_ylim(Y_center-rsvrprop.BYL/2, Y_center+rsvrprop.BYL/2)
        
        ax.add_patch(plt.Rectangle((X_center-rsvrprop.BXL/2, Y_center-rsvrprop.BYL/2), 
                                   rsvrprop.BXL, rsvrprop.BYL, fill=False, linestyle='--', linewidth=2, edgecolor='b'))
        
        
        for i in range(nReservoirs):
            x = (AllReservoirs_m[:, 2*i]) 
            y = (AllReservoirs_m[:, 2*i+1])
            plt.plot(x, y, c='k')
        if UnitSys == 2:  
            plt.title('CO₂ Plume Extension')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
        else:
            plt.title('CO₂ Plume Extension')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
        fig_plume.set_size_inches(8, 6)
        plt.axis('equal')
        
        df_aor = pd.DataFrame(highlight_locations, columns=['x', 'y'])
        #df_aor['x_UTM'] = df_aor['x'] - shift_x/lencoef
        #df_aor['y_UTM'] = df_aor['y'] - shift_y/lencoef
        
        return fig_plume, df_aor
    fig_plume, df_aor = Plumes_map()
    
    
    
    def Geometry_output():

        PressureCheck = [' '] * (nWE+nWI)  # Initialize PressureCheck with 'P' for all elements
        FailureCounter = 0

        for i in range(nWI):
            if Psave2[i] > MaxAllowedP[i] :
                PressureCheck[i] ='F'
                FailureCounter += 1
            else:
                PressureCheck[i] =  'P'
        for i in range(nWI, nWI + nWE):
            if Psave2[i] < MinAllowedP[i-nWI]:
                PressureCheck[i] = 'F'
                FailureCounter += 1
            else:
                PressureCheck[i] = 'P'

        PressureCheck = np.reshape(PressureCheck, (-1, 1))  # Convert to vertical column

        
        with open('EASiToolOutput_Geometry.txt', 'w') as file:        
            file.write('*********************************************************************\n')
            file.write('*********************************************************************\n')
            file.write('%70s\n' % ('EASiTool Output File for ' + proj_name).center(70))
            file.write('*********************************************************************\n\n\n\n')
            file.write('*********************************************************************\n')
            file.write('%43s\n' % 'EASiTool Input Data')
            file.write('*********************************************************************\n')  
            for column in output_print.columns:
                file.write(
                    '%43s %25s %20s\n' % (column, str(output_print.loc[0, column]), str(output_print.loc[1, column])))
              
            
            if BCValue == 1:
                file.write('%43s %25s\n\n' % ('Boundary Condition', 'Closed'))
            elif BCValue == 2:
                file.write('%43s %25s\n\n' % ('Boundary Condition', 'Open'))
            
            
            if UnitSys == 1:
                file.write('%43s %25s\n\n' % ('Unit System', 'Field Units'))
            elif UnitSys ==2:
                file.write('%43s %25s\n\n' % ('Unit System', 'SI Units'))       
                
            file.write('%43s %25s\n\n' % ('Injection/Extraction Type', 'General Geometry/Pattern'))
            
                
            file.write('\n')
            file.write('*********************************************************************\n') 
            file.write('%43s\n' % 'EASiTool Output Data')
            file.write('*********************************************************************\n') 
            file.write('%20s%20s\n' % ('Capacity (MMT):', 'NPV ($M):'))
            file.write('%20.3e%20.3e\n' % (Capacity[0], NPV[0]))
            file.write('*********************************************************************\n')
                
      
            file.write('**************************************************************************************************\n')
            file.write('%41s %5u\n' % ('Number of Injectors =', nWI))
            file.write('%41s %5u\n' % ('Number of Extractors =', nWE))
            file.write('**************************************************************************************************\n')
            if UnitSys == 2: # SI units 
                file.write('%-5s %-7s %-12s %-12s %-15s %-15s %-10s %s %s\n' % ('Well', 'Type', 'X', 'Y', 'Inj. Rate', 'Ext. Rate', 'Pressure', 'Plume radius', 'Pres. Check'))
                file.write('%-5s %-7s %-12s %-12s %-15s %-15s %-10s %s %s\n' % ('#', '', '(UTM)', '(UTM)', '(MMT/yr)', '(m^3/day)', '(MPa)', '(m)', ''))
                file.write('%s\n' % '*************************************************************************************************')
                for nW in range(nWI):
                    file.write('%-5u %-7s %-12.3e %-12.3e %-15.3e %-15s %-10.3e %-10.3e %s\n' % ((nW+1), 'Inj.', Xwell_input[nW], Ywell_input[nW], Qsave2_output[nW], 'N/A', Psave2_output[nW], rL1[nW][0], PressureCheck[nW][0]))
                    
                for nW in range(nWI, nWI+nWE):
                    file.write('%-5u %-7s %-12.3e %-12.3e %-15s %-15.3e %-10.3e %-10s %s\n' % ((nW+1), 'Ext.', Xwell_input[nW], Ywell_input[nW], 'N/A', Qsave2_output[nW], Psave2_output[nW], 'N/A', PressureCheck[nW][0]))
                file.write('\n')
                file.write('*********************************************************************\n') 
                file.write('%43s\n' % 'Basin Pressure (MPa) Data (50 by 50)')
                file.write('*********************************************************************\n') 
                # Reshape the list into a 50x50 matrix
                P_matrix = np.reshape(PContour_output, (nWX, nWY))
                for row in P_matrix:
                    rounded_row = [round(value, 2) for value in row]
                    file.write(" ".join(map(str, rounded_row)) + "\n")
                    
            elif UnitSys == 1:  # Field Units
                file.write('%-5s %-7s %-13s %-13s %-16s %-14s %-11s %-13s %-13s %-11s\n' % ('Well', 'Type', 'X', 'Y', 'Inj. Rate', 'Ext. Rate', 'Pressure', 'Plume radius', 'Plume radius', 'Pres. Check'))
                file.write('%-5s %-7s %-13s %-13s %-16s %-14s %-11s %-13s %-13s %-11s\n' % ('#', '', '(UTM)', '(UTM)', '(MMT/yr)', '(bbl/day)', '(psi)', '(m)', '(ft)', ''))
                file.write('%s\n' % '**************************************************************************************************')
                
                for nW in range(nWI):
                    file.write('%-5u %-7s %-13.3e %-13.3e %-16.3e %-14s %-11.3e %-13.3e %-13.3e %-11s\n' % ((nW + 1), 'Inj.', Xwell_input[nW], Ywell_input[nW], Qsave2_output[nW], 'N/A', Psave2_output[nW], rL1[nW][0], rL1[nW][0], PressureCheck[nW][0]))
                    
                for nW in range(nWI, nWI+nWE):
                    file.write('%-5u %-7s %-13.3e %-13.3e %-16s %-14.3e %-11.3e %-13s %-13s %-11s\n' % ((nW + 1), 'Ext.', Xwell_input[nW], Ywell_input[nW], 'N/A', Qsave2_output[nW], Psave2_output[nW], 'N/A', 'N/A', PressureCheck[nW][0]))
            
                file.write('\n')
                file.write('*********************************************************************\n') 
                file.write('%43s\n' % 'Basin Pressure (psi) Data (50 by 50)')
                file.write('*********************************************************************\n') 
                # Reshape the list into a 50x50 matrix
                P_matrix = np.reshape((np.array(PContour_output) * 145.03768), (nWX, nWY))
                for row in P_matrix:
                    rounded_row = [round(value, 2) for value in row]
                    file.write(" ".join(map(str, rounded_row)) + "\n")
            
            else: 
                raise Exception('Wrong Unit System Selection!')
             

                file.write('%s\n\n\n\n\n' % ' ')
            
           
               
        file.close()
             
    Geometry_output()
    
    return Psave2_output, Qsave2_output, Capacity[0], NPV[0], fig_contour, fig_plume, highlight_locations, shift_x, shift_y, df_aor

def flatten_list(lst):
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened


def S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi):
    # This function is the core algorithm of S3, and return the final pressure distribution and NPV prediction for sensitivity analysis
    
    TotalInjRate = [0] * simprop.nWXmax
    Capacity = [0] * simprop.nWXmax
    NPV: Any= [0] * simprop.nWXmax
    
    if UnitSys == 2:
        lencoef = 1000
    elif UnitSys == 1:
        lencoef = 1609.344
    else: 
        raise Exception ("Wrong Unit System selection!")

    df_inj = pd.read_excel(uploaded_file,sheet_name='Injectors')
    df_ext = pd.read_excel(uploaded_file,sheet_name='Extractors')
    
    WellNumber_inj = df_inj["Well Number"]
    WellNumber_ext = df_ext["Well Number"]
    
    #Xwell_UTM = df_inj["UTM_Well_X"]
    #Ywell_UTM = df_inj["UTM_Well_Y"]
    
    Xwell = df_inj["UTM_Well_X"]
    Xwell = pd.concat([df_inj["UTM_Well_X"], df_ext["UTM_Well_X"]], ignore_index=True)
    Ywell = df_inj["UTM_Well_Y"]
    Ywell = pd.concat([df_inj["UTM_Well_Y"], df_ext["UTM_Well_Y"]], ignore_index=True)      
    
    Xwell_input = Xwell.copy()
    Ywell_input = Ywell.copy()

    AllReservoirs_orig = pd.read_excel(uploaded_file, sheet_name='Geometry')
    
    # Convert AllReservoirs values to km
    AllReservoirs_nonan = AllReservoirs_orig.dropna(axis=1, how='all')
    AllReservoirs_m = AllReservoirs_nonan.values  #in meter
    AllReservoirs = AllReservoirs_m / lencoef   # in km or mile
    # Get the number of reservoirs
    nReservoirs = AllReservoirs.shape[1] // 2
    
    # To store all reservoirs' geometry locations
    Xres = []
    Yres = []
    for i in range(nReservoirs):
        x = (AllReservoirs_m[:, 2*i]) 
        y = (AllReservoirs_m[:, 2*i+1])
        Xres.extend(x)
        Yres.extend(y)
    Xres = np.array(Xres)
    Yres = np.array(Yres)
    Xres = Xres[~np.isnan(Xres)]
    Yres = Yres[~np.isnan(Yres)]                   
    
    try: 
        if UnitSys == 2:
            InjRate = df_inj["Injection Rate (Ton/day)"]
            ExtRate = df_ext["Extraction Rate (m^3/day)"]
            #MaxAllowedP = df_inj["Max Injection Pressure (Mpa)"]
            #MinAllowedP = df_ext["Min Extraction Pressure (Mpa)"]
        elif UnitSys == 1:
            InjRate = df_inj["Injection Rate (MMT/yr)"] * 2739.726027
            ExtRate = df_ext["Extraction Rate (bbl/day)"] * 0.158987
            #MaxAllowedP = df_inj["Max Injection Pressure (psi)"] * 0.00689476
            #MinAllowedP = df_ext["Min Extraction Pressure (psi)"] * 0.00689476
    except (KeyError, ValueError, ZeroDivisionError) as e:
        st.warning('Error in column names or data. Please check your input units. :exclamation:')
        sys.exit()
        
    simprop.rE = math.sqrt(rsvrprop.BXL*rsvrprop.BYL/math.pi)

    Area = rsvrprop.XL * rsvrprop.YL # Project Area
    
    # Shift well and reservoir locations based on the new center, from previous rectangle's center to average well loc coordiates in UTM. 
    
    combined_X = np.concatenate([Xwell.values, Xres])
    combined_Y = np.concatenate([Ywell.values, Yres])
    
    X_center = sts.mean(combined_X)
    Y_center = sts.mean(combined_Y)
    
    # To get the Basin boundary based on the well locations
    try:
       GeoRatio = (max(combined_X) - min(combined_X)) / (max(combined_Y) - min(combined_Y))
    except ZeroDivisionError:
       GeoRatio = 1  

       
   # Check if Xwell and Ywell have only 1 element
    if len(combined_X) == 1 or len(combined_Y) == 1:
        GeoRatio = 1
       
    # Make sure the shape of reservoir is not too strange. 
    if GeoRatio < 0.5: 
       GeoRatio = 0.5
    if GeoRatio > 2:
       GeoRatio = 2

        
    Area_Basin = rsvrprop.BXL*rsvrprop.BYL # Based on inputs, both reservoir and basin are set to be square first, then adjusted into rectangles
    
    rsvrprop.BXL = math.sqrt(GeoRatio*Area_Basin)
    rsvrprop.BYL = rsvrprop.BXL/GeoRatio
    
    # In case the basin cannot cover all the project area. 
    long_x_half = max(X_center-min(combined_X),max(combined_X)-X_center)
    long_y_half = max(Y_center-min(combined_Y),max(combined_Y)-Y_center) 
    if long_x_half>rsvrprop.BXL/2 or long_y_half>rsvrprop.BYL/2:
        rsvrprop.BXL = max(combined_X)-min(combined_X)
        rsvrprop.BYL = max(combined_Y)-min(combined_Y)
        X_center = (max(combined_X)+min(combined_X))/2
        Y_center = (max(combined_Y)+min(combined_Y))/2
        #Area_Basin_new = rsvrprop.BXL*rsvrprop.BYL/lencoef/lencoef
        #st.warning('The Sensitivity Analysis is based on the adjusted reservoir area shown in Output Figures.')
        
    
    shift_x = rsvrprop.BXL/2 - X_center
    shift_y = rsvrprop.BYL/2 - Y_center
    

    Xwell_shift = Xwell + shift_x
    Ywell_shift = Ywell + shift_y
    
    Xwell = Xwell_shift.copy()
    Ywell = Ywell_shift.copy()
    
   
    
    nWT = len(WellNumber_inj)+len(WellNumber_ext)
    nWI = len(WellNumber_inj)
    
    
    
    simprop.nWE = len(WellNumber_ext)
    nWE = simprop.nWE
    
    
    
    Psave = np.zeros(nWT,order = 'F')
    chi_BL = np.zeros(nWI,order = 'F')
    chi_dry = np.zeros(nWI,order = 'F')
    Qinj = np.zeros(nWI,order = 'F')
    Qext = np.zeros(simprop.nWE,order = 'F')
    Rwell = np.zeros((nWT, nWT),order = 'F')
    rT1 = np.zeros((nWT, 1),order = 'F')
    rL1 = np.zeros((nWT, 1),order = 'F')
    Qsave2 = np.zeros((nWT, 1),order = 'F')
    Psave2 = np.zeros((nWT, 1),order = 'F')
    Qsave2_output = np.zeros((nWT, 1),order = 'F')
    Psave2_output = np.zeros((nWT, 1),order = 'F')
    
    simprop.nWXmin=1;
    simprop.nWXmax=1;
    
    # Beginning of Distance matrix
    
    if BCValue == 1:
        rED1 = simprop.rE / simprop.rW  # closed
        rED2 = simprop.rE / simprop.rW  # closed
    elif BCValue == 2:
        rED1 = 1e38  # open
        rED2 = 1e38  # open
    else:
        raise ValueError('Error: invalid BCValue.') # or use an appropriate exception type and error message
    
    
    #for i in range(nWI + simprop.nWE):
    #    if Xwell[i] < 0 or Xwell[i] > rsvrprop.BXL or Ywell[i] < 0 or Ywell[i] > rsvrprop.BYL:
    #        st.warning('Error: Wells should be located inside the reservoir.')
    #        sys.exit()
    
    
    for i in range (nWT):
        X=Xwell[i]
        Y=Ywell[i]
        for j in range (nWT):
            Rwell[i,j]= math.sqrt((X-Xwell[j])**2+(Y-Ywell[j])**2)
    
    # End of Distance matrix 
    
    P_ave = rsvrprop.P0
    err = 1
    it = 1
    nWX = 1
    
    while err > 1e-2:
    
        
        Fluid = Fluid_Prop(P_ave, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)
    
        Fluid.Spycher()
    
        mug, cg, rhog, rhoc = Fluid.CO2Prop()
    
        cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
    
        tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
        
        Qinj[:nWI]=InjRate[:nWI]/86400/rhoc*1000 # rm^3/s
        Qext[:simprop.nWE]=ExtRate[:simprop.nWE]/86400/rhob*rhobstd # rm^3/s
        epsilon=np.multiply(Qinj,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
        chi_BL[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgL)
        chi_dry[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgT)
        
    
        A = A_Geometry(nWX, simprop.nWE, nWI, simprop.rW, chi_BL, chi_dry, Lamda_g,
                            Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue, Rwell)
        
        P, P_ave, Psave, err, it = B_Geometry(A, simprop.nWE, nWI, Qinj,Qext,mug,mub,rsvrprop.Thickness,rsvrprop.k,relakprop.kra0,relakprop.krg0,rsvrprop.P0,Psave,it)
#        print(err)
        if it>200:
            break
        
    
    for i in range(nWI):
        rT1[i, nWX-1] = np.sqrt(zT * Qinj[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
        rL1[i, nWX-1] = np.sqrt(zL * Qinj[i] * simprop.SimTime / np.pi / rsvrprop.Porosity / rsvrprop.Thickness)
        Qsave2[i, nWX-1] = InjRate[i] # ton per day
        Psave2[i, nWX-1] = P[i]
        
    Psave2_ave = np.mean(Psave2[0:nWI, nWX-1])
        
    for i in range(nWI, nWT):
#        Qsave2[i, nWX-1] = -1 * ExtRate[i-nWI] * rhobstd / 1000 # ton per day
        Qsave2[i, nWX-1] = -1 * ExtRate[i-nWI]  # m^3 per day
        Psave2[i, nWX-1] = P[i]
        
    
    L = 0    
    for i in range(nWI):
        L += InjRate[i]
        
    TotalInjRate[nWX-1] = L # ton per day
    Capacity[nWX-1] = TotalInjRate[nWX-1] * simprop.SimTime / 86400 / 1000000 # Mton

    NPV[nWX - 1] = 0
    for years in range(1, int(simprop.SimTime / 365 / 24 / 60 / 60) + 1):
        cashflow = (Capacity[nWX - 1] / int(
            simprop.SimTime / 365 / 24 / 60 / 60) * 1000000 * npvprop.txcr - npvprop.opcostup*1000000 - npvprop.opcostdwn * 1000000*(
                                simprop.nWE + nWI))
        annual_npv = cashflow / (1 + npvprop.disrate) ** years
        NPV[nWX - 1] = NPV[nWX - 1] + annual_npv
    NPV[nWX - 1] = (NPV[nWX - 1] - npvprop.initinv*1000000) / 1000000  # Million $
    
    def plot_multi_AOR_geo():
        nWX = 50
        nWY = 50
        rW = simprop.rW
        nNodes = nWX*nWY   
        
        Pnode = np.zeros(nNodes,order = 'F')
        #XContour = np.zeros(nNodes + nWI + nWE, 'F')
        #YContour = np.zeros(nNodes + nWI + nWE, 'F')
        #PContour = np.zeros(nNodes + nWI + nWE, 'F')
        chi_BL = np.zeros(nWI,order = 'F')
        chi_dry = np.zeros(nWI,order = 'F')
        
        X = []
        Y = []
        
        for i in range(nWY):
            for j in range(nWX):
                X.append(rsvrprop.BXL / nWX / 2 + j * rsvrprop.BXL / nWX)
                Y.append(rsvrprop.BYL / nWX / 2 + i * rsvrprop.BYL / nWY)
                
        Rwell_comp = np.zeros((nNodes, nWI + nWE),'F')
        A = np.zeros((nNodes,nWI + nWE),'F')
        
    
        for i in range(nNodes):
            for j in range(nWI + nWE):
                Rwell_comp[i, j] = np.sqrt((X[i] - Xwell[j]) ** 2 + (Y[i] - Ywell[j]) ** 2)
        Rwell = np.real(Rwell_comp) # In case complex number occurs
                
        Qinj = np.zeros(nWI,order = 'F')
        Qext = np.zeros(simprop.nWE,order = 'F')
        Sa = np.zeros(nWI)
        
        Qinj[:nWI]=InjRate[:nWI]/86400/rhoc*1000 # rm^3/s
        Qext[:simprop.nWE]=ExtRate[:simprop.nWE]/86400/rhob*rhobstd # rm^3/s
    
        epsilon=np.multiply(Qinj,(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
        chi_BL[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgL)
        chi_dry[:nWI]=np.full(nWI,(1/4)*epsilon*dfgdSgT)
        
        # Injectors to Injectors
        
        for i in range(nNodes):
            for j in range(nWI):
                rD = Rwell[i,j]/rW
                
                #Sa = np.zeros(nWI)
                Sa[j] = 0.5*np.log(chi_dry[j]) + (0.5/F_Lg)*np.log(chi_BL[j]/chi_dry[j]) - 0.5*(Lamda_g/Lamda_w)*np.log(chi_BL[j]/eta_D3) + 0.5*0.577215*(1-Lamda_g/Lamda_w)
                
                if tD <= rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(
                        rD ** 2 / 4 / tD / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD >= rD ** 2 / 4 / np.mean(chi_dry):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD) - 0.5 * exp1(chi_dry[j]) + 1 / 2 / F_Lg * exp1(
                        chi_dry[j] / F_Lg / eta_D2) - 1 / 2 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD < rD ** 2 / 4 / np.mean(chi_dry) and tD > rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD / F_Lg / eta_D2) - 0.5 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
    
      
            for j in range(nWI, nWI + nWE):
                rDE = Rwell[i, j] / rW
                if BCValue == 2:
                    A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
                elif BCValue == 1:
                    A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4
                    
        real_A = np.real(A)
        real_A_nonzero = real_A[real_A.imag != 0]  
        
    
    
        #1. Pressure Contour map
        if real_A_nonzero.size != 0:
            raise ValueError("Matrix A contains complex number. Exception raised.")
        else:
            A = np.copy(real_A)
            
        for i in range(nNodes):
            BInj = 0
            for j in range(nWI):
                BInj += A[i,j] * Qinj[j]
            BExt = 0
            for j in range(nWI, nWI + nWE):
                BExt += A[i,j] * Qext[j - nWI]
            Pnode[i] = mug / (2 * math.pi * rsvrprop.Thickness * rsvrprop.k * 
                              relakprop.krg0) * BInj / 1000000 - mub / (2 * math.pi * rsvrprop.Thickness 
                            * rsvrprop.k * relakprop.kra0) * BExt / 1000000 + rsvrprop.P0 # Krg0 to Krs
    
        # Prepare contour data for nodes
        XContour = [x-shift_x  for x in X]
        YContour = [y-shift_y  for y in Y]
        PContour = Pnode.tolist()
        PContour_output = PContour.copy()

        # Prepare contour data for wells
        for i in range(nWI + nWE):
            XContour.append(Xwell_input[i])
            YContour.append(Ywell_input[i])
            PContour.append(Psave2[i].tolist())
        
            
        #Interpolate data onto a grid
        #x1 = np.arange(0, rsvrprop.BXL/lencoef, 1/200*rsvrprop.BXL/lencoef)
        #x2 = np.arange(0, rsvrprop.BYL/lencoef, 1/200*rsvrprop.BYL/lencoef)
        x1 = np.arange(X_center-rsvrprop.BXL/2, X_center+rsvrprop.BXL/2, 1/200*rsvrprop.BXL)
        x2 = np.arange(Y_center-rsvrprop.BYL/2, Y_center+rsvrprop.BYL/2, 1/200*rsvrprop.BYL)
        
        X1, X2 = np.meshgrid(x1, x2)
        PContour_flatten = flatten_list(PContour) # Make PContour in a list format
        highlight_AOR = rsvrprop.P0 + simprop.CP  # Specify the contour level to highlight
        

            
        # UNIT CONVERSION
        if UnitSys == 1: # field unit, then output pressure is in psi
            PContour_flatten_unit = [x * 145.03768 for x in PContour_flatten] # Convert from MPa to psi
            Psave2_output = Psave2 * 145.03768 # Convert from MPa to psi
            highlight_AOR_unit = highlight_AOR * 145.03768 # Convert from MPa to psi
            Qsave2_output[:nWI] = Qsave2[:nWI] * 0.000365 #Injection well convert from ton/day to MMT/year
            Qsave2_output[nWI:nWI+nWE] = Qsave2[nWI:nWI+nWE] * 6.2898107 # Extraction well convert from m^3/day to bbl/day
            
        else: # SI unit
            PContour_flatten_unit = [x * 1 for x in PContour_flatten] #Output pressure still in MPa
            Psave2_output = Psave2  # Output pressure still in MPa
            highlight_AOR_unit = highlight_AOR  # in MPa
            Qsave2_output[:nWI] = Qsave2[:nWI] * 0.000365 #Injection well convert from ton/day to MMT/year
            Qsave2_output[nWI:nWI+nWE] = Qsave2[nWI:nWI+nWE] # Extraction well in m^3/day
        
        # Capture warning messages
        #if min(PContour_flatten_unit) > highlight_AOR_unit or max(PContour_flatten_unit) < highlight_AOR_unit:
        #       warnings.warn('Critical Pressure is not suitable! AOR cannot be visualized!')
        
        foo = griddata((XContour, YContour), PContour_flatten_unit, (X1, X2), method='linear') #200 by 200 pressure matrix
        
        # Now foo is the interpolated data on the grid defined by X1 and X2
        X3 = foo  
        
        # Create a contour plot with 10 levels
        fig_contour, ax = plt.subplots()
        
        contour = ax.contourf(X1, X2, X3, 10)
        
        highlight_contour = ax.contour(X1, X2, X3, levels=[highlight_AOR_unit], colors='red', linewidths=2) # Mark AOR
        
        highlight_locations = []
        for path_collection in highlight_contour.collections:
            paths = path_collection.get_paths()
            for path in paths:
                vertices = path.vertices
                for x, y in vertices:
                    highlight_locations.append((x, y))

    
        # Set plot aspect ratio, title, and axis labels
        plt.gca().set_aspect('equal', adjustable='box')
        cbar = plt.colorbar(contour)
        
        # Add labels to the highlighted contour lines
        try:
            plt.clabel(highlight_contour, inline=True, fmt='AOR', colors='red', fontsize=10)
        except Exception as e:
            # Handle any exception that occurs
            st.warning("Warning: Cannot plot the desired AOR within the reservoir boundaries in Sensitivity Analysis, the estimated AOR is too large (or too small).")
        
        if UnitSys == 2:
            plt.title('Pressure Contour, MPa')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
            cbar.set_label('Pressure, MPa')
        elif UnitSys ==1:
            plt.title('Pressure Contour, psi')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
            cbar.set_label('Pressure, psi')
        
        plt.axis('equal')
        
        return highlight_locations

          
    if Sensi == 1:
        highlight_locations = plot_multi_AOR_geo()
    else:
        highlight_locations = []
    
    return Psave2_ave, NPV[0], highlight_locations, X_center, Y_center, rL1, rT1


'''
## TESTING
from properties import Unit_conversion, read_rsvr_input, userexcel_conversion
import pandas as pd
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_SI_testing.xlsx'
UnitSys = 2
BCValue = 2
Sensi = 0

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
userinput_inj = pd.read_excel(uploaded_file, sheet_name='Injectors', header=None)
userinput_ext = pd.read_excel(uploaded_file, sheet_name='Extractors', header=None)
input_df, proj_name = userexcel_conversion(userinput_df, userinput_inj, userinput_ext)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)


Psave2_output, Qsave2_output, Capacity, NPV, fig_contour, fig_plume, highlight_locations, shift_x, shift_y, df_aor= Scenario3(rsvrprop, simprop, relakprop, npvprop, BCValue,uploaded_file,UnitSys, output_print, proj_name) 
'''