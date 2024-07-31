
"""
Created on Mon May  1 16:28:25 2023

@author: wangz
"""

# --Import libraries--

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import streamlit as st
import warnings
import math
from Matcal import Distance_Mat_Builder
from reservoir_funcs import Fluid_Prop
from scipy.special import exp1
import pandas as pd
import plotly.graph_objects as go

def warning_pressure_plume(nWXmin,nWXmax, rL1, Xwell, Ywell, Rwell, ConRate, Psave2, Qsave2, rsvrprop, nWE, nWX_ctrl):
    WarningCountPlume = np.zeros(nWXmax,order = 'F')
    WarningPlume = 0
    BXL = rsvrprop.BXL
    BYL = rsvrprop.BYL
    XL = rsvrprop.XL
    YL = rsvrprop.YL
    for nWX in range(1, nWXmax + 1):
        nWY = nWX
        for i in range(1, nWX + 1):           
            if rL1[i - 1][nWX - 1] > (Ywell[i - 1][nWX - 1] - BYL / 2 + YL / 2):
                WarningCountPlume[nWX - 1] += 1
                WarningPlume += 1
            if rL1[(nWX) * (nWY - 1) + i - 1][nWX - 1] > (BYL / 2 + YL / 2 - Ywell[(nWX) * (nWY - 1) + i - 1][nWX - 1]):
                WarningCountPlume[nWX - 1] += 1
                WarningPlume += 1
    
        for j in range(1, nWY):
            for i in range(1, nWX + 1):
                if (rL1[(j - 1) * nWX + i - 1][nWX - 1] + rL1[j * nWX + i - 1][nWX - 1]) > Rwell[(j - 1) * nWX + i - 1][(j * nWX + i - 1)][nWX - 1]:
                    WarningCountPlume[nWX - 1] += 1
                    WarningPlume += 1
    
        for i in range(1, nWY + 1):
            if rL1[(i - 1) * nWX][nWY - 1] > (Xwell[(i - 1) * nWX][nWY - 1] - BXL / 2 + XL / 2):
                WarningCountPlume[nWY - 1] += 1
                WarningPlume += 1
            if rL1[i * nWX - 1][nWX - 1] > (BXL / 2 + XL / 2 - Xwell[i * nWX - 1][nWY - 1]):
                WarningCountPlume[nWY - 1] += 1
                WarningPlume += 1
    
        for j in range(1, nWY + 1):
            for i in range(1, nWX):
                if (rL1[(j - 1) * nWX + i - 1][nWY - 1] + rL1[(j - 1) * nWX + i][nWX - 1]) > Rwell[(j - 1) * nWX + i - 1][(j - 1) * nWX + i][nWX - 1]:
                    WarningCountPlume[nWY - 1] += 1
                    WarningPlume += 1
    #return WarningCountPlume, WarningPlume

#def calculate_warning_pressure(nWXmin, nWXmax, ConRate, Psave2, rsvrprop, Qsave2, nWE):

    WarningCountPressure = np.zeros(nWXmax)
    WarningPressure = 0
    FracP = rsvrprop.FracP
    P0 = rsvrprop.P0
    
    if ConRate == 1:
        for j in range(nWXmin, nWXmax+1):
            for i in range(j*j):
                if Psave2[i,j-1] > FracP:
                    WarningCountPressure[j-1] += 1
                    WarningPressure += 1
            for i in range(j*j, j*j+nWE):
                if Psave2[i,j-1] < 0.5*P0:
                    WarningCountPressure[j-1] += 1
                    WarningPressure += 1
    else:
        for j in range(nWXmin, nWXmax+1):
            for i in range(j*j):
                if Qsave2[i,j-1] > 2000: #ton/day of CO₂, adjusted accordingly
                    WarningCountPressure[j-1] += 1
                    WarningPressure += 1
            for i in range(j*j, j*j+nWE):
                if Qsave2[i,j-1] > 2000: #m^3/day of brine, adjusted accordingly
                    WarningCountPressure[j-1] += 1
                    WarningPressure += 1
    #return WarningCountPressure, WarningPressure

#def plume_pressure_warning(WarningCountPlume,WarningCountPressure,WarningPressure,WarningPlume, nWX_ctrl):
    if WarningCountPressure[nWX_ctrl-1] == 0 and WarningCountPlume[nWX_ctrl-1] > 0:
        st.warning(f'Design Considerations for {nWX_ctrl**2} injectors {nWE} extractors :\n\n'
                   '- CO₂ plume of at least one injector crosses the reservoir boundary or overlaps another plume in the selected well patterns. \n\n'
                   'Consider adjusting these parameters:\n'
                   '- Maximum Allowable injection pressure - Minimum extraction pressure - Injection rate - Extraction rate.')
    elif ConRate == 1 and  WarningCountPressure[nWX_ctrl-1] > 0 and WarningCountPlume[nWX_ctrl-1] == 0:
        st.warning(f'Design Considerations for {nWX_ctrl**2} injectors {nWE} extractors :\n\n'
                   '- Bottomhole pressure of at least one injector exceeds the frac pressure in this selected well patterns.\n'
                   ' and/or Bottomhole pressure of at least one extractor drops below 50% of the initial pressure in the selected well patterns.\n\n'
                   'Consider adjusting these parameters:\n'
                   '- Injection rate - Extraction rate.')
    elif ConRate == 1 and WarningCountPressure[nWX_ctrl-1] > 0 and WarningCountPlume[nWX_ctrl-1] > 0:
        st.warning(f'Design Considerations for {nWX_ctrl**2} injectors {nWE} extractors :\n\n'
                   '- Bottomhole pressure of at least one injector exceeds the frac pressure in this selected well patterns.\n'
                   ' and/or Bottomhole pressure of at least one extractor drops below 50% of the initial pressure in the selected well patterns.\n'
                   '- Also, CO₂ plume of at least one injector crosses the reservoir boundary or overlaps another plume in the selected well patterns.\n\n'
                   'Consider adjusting these parameters:\n'
                   '- Injection rate - Extraction rate.')
    elif ConRate == 0 and WarningCountPressure[nWX_ctrl-1] > 0 and WarningCountPlume[nWX_ctrl-1] == 0:
        st.warning(f'Design Considerations for {nWX_ctrl**2} injectors {nWE} extractors :\n\n'
                   '- Injection rate of at least one injector exceeds 2000 ton/day (0.73MMT/yr) in this selected well patterns.\n'
                   ' and/or Extraction rate of at least one extractor exceeds 2000 m^3/day (12580bbl/day) in the selected well patterns.\n\n'
                   'Consider adjusting these parameters:\n'
                   '- Maximum Allowable injection pressure - Minimum extraction pressure.')
    elif ConRate == 0 and WarningCountPressure[nWX_ctrl-1] > 0 and WarningCountPlume[nWX_ctrl-1] > 0:
        st.warning(f'Design Considerations for {nWX_ctrl**2} injectors {nWE} extractors :\n\n'
                   '- Injection rate of at least one injector exceeds 2000 ton/day (0.73MMT/yr) in this selected well patterns.\n'
                   ' and/or Extraction rate of at least one extractor exceeds 2000 m^3/day (12580bbl/day) in the selected well patterns.\n'
                   '- Also, CO₂ plume of at least one injector crosses the reservoir boundary or overlaps another plume in the selected well patterns.\n\n'
                   'Consider adjusting these parameters:\n'
                   '- Maximum Allowable injection pressure - Minimum extraction pressure.')
    else: 
        st.success(f"😄 The design in the selected well pattern ({nWX_ctrl**2} injectors and {nWE} extractors) is theoretically feasible and operational.")
            

def plot_npv_capacity(simprop, Capacity, NPV):
    
    nWXmin = simprop.nWXmin
    nWXmax = simprop.nWXmax
    ### Format the NPV array as a table for display
    npv_table = pd.DataFrame({'nWX': range(nWXmin, (nWXmax + 1)), 'NPV $M': NPV[(nWXmin - 1):(nWXmax)]})
    npv_table['Number of Injection Wells'] = npv_table['nWX'] ** 2
    Capacity_table = pd.DataFrame(
        {'nWX': range(nWXmin, (nWXmax + 1)), 'Capacity MMT of CO₂': Capacity[(nWXmin - 1):(nWXmax)]})
    Capacity_table['Number of Injection Wells'] = Capacity_table['nWX'] ** 2
    # st.header('Output Results')
    col1, col2 = st.columns(2)

    with col1:
        #st.header('Capacity')
        # Convert the DataFrame to a Plotly figure
        capacity_fig = go.Figure(data=go.Scatter(x=Capacity_table['Number of Injection Wells'],
                                                 y=Capacity_table['Capacity MMT of CO₂'],
                                                 mode='lines+markers'))
        capacity_fig.update_layout(title="Capacity vs Number of Injection Wells",
                                   xaxis_title="Number of Injection Wells",
                                   yaxis_title="Capacity (MMT of CO₂)")
        st.plotly_chart(capacity_fig)     

    with col2: 
        #st.header('NPV')
        # Convert the DataFrame to a Plotly figure
        npv_fig = go.Figure(
            data=go.Scatter(x=npv_table['Number of Injection Wells'], y=npv_table['NPV $M'],
                            mode='lines+markers'))
        npv_fig.update_layout(title="NPV vs Number of Injection Wells",
                              xaxis_title="Number of Injection Wells", yaxis_title="NPV ($M)")
        st.plotly_chart(npv_fig)
    return npv_table, Capacity_table



def plot_plume_extension(nWX, Xwell, Ywell, rL1, nWE, rsvrprop, UnitSys):
    
    BXL = rsvrprop.BXL
    BYL = rsvrprop.BYL
    XL = rsvrprop.XL
    YL = rsvrprop.YL
    
    
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    
   
    if UnitSys == 2:
        lencoef = 1000
    elif UnitSys == 1:
        lencoef = 1609.34
    else: 
        raise Exception ("Wrong Unit System selection!")
    # Draw circles     
    for w in range(nWX*nWX):
        x = Xwell[w,nWX-1]/lencoef
        y = Ywell[w,nWX-1]/lencoef
        r = rL1[w,nWX-1]/lencoef

        circle = plt.Circle((x, y), radius=r, color='r', fill=False)
        ax.add_artist(circle)

    # Draw scatter points
    if nWE > 0: 
        for w in range(nWE):
            x = Xwell[w+nWX*nWX,nWX-1]/lencoef
            y = Ywell[w+nWX*nWX,nWX-1]/lencoef
            ax.scatter(x, y, color='b', marker='^')

    # Draw rectangles
    rect1 = plt.Rectangle(((BXL-XL)/2/lencoef, (BYL-YL)/2/lencoef), XL/lencoef, YL/lencoef, linewidth=2, linestyle='--', edgecolor='g', facecolor='none', label='Project Area')
    ax.add_patch(rect1)
    rect2 = plt.Rectangle((0, 0), BXL/lencoef, BYL/lencoef, linewidth=2, linestyle='--', edgecolor='b', facecolor='none', label='Reservoir Area')
    ax.add_patch(rect2)
    plt.legend(loc='upper right')
    # Set title, axis labels, and background color
    ax.set_title('CO₂ Plume Extension')
    if UnitSys ==2:
        ax.set_xlabel('X , km')
        ax.set_ylabel('Y , km')
    else: 
        ax.set_xlabel('X , mile')
        ax.set_ylabel('Y , mile')
    fig.set_facecolor('w')
    ax.set_aspect('equal')
    ax.set_xlim([-BXL/lencoef*0.05, BXL/lencoef*1.05])
    ax.set_ylim([-BYL/lencoef*0.05, BYL/lencoef*1.05])

    return fig


def plot_well_property(Xwell, Ywell, nWX, nWE, ConRate, Qsave2, Psave2, UnitSys, rsvrprop):
    
    BXL = rsvrprop.BXL
    BYL = rsvrprop.BYL
    
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    if UnitSys == 2:
        lencoef = 1000
        pcoef = 1
    elif UnitSys == 1:
        lencoef = 1609.34
        pcoef = 145.038 # convert from MPa to psi
    else: 
        raise Exception ("Wrong Unit System selection!")
    if ConRate == 0:
        if UnitSys == 2: # SI units
            #vmin = min(np.min(Qsave2[0:nWX*nWX, nWX-1]), np.min(Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
            #vmax = max(np.max(Qsave2[0:nWX*nWX, nWX-1]), np.max(Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))     
            p = ax.scatter(Xwell[0:nWX*nWX, nWX-1]/1000, Ywell[0:nWX*nWX, nWX-1]/1000, c=Qsave2[0:nWX*nWX, nWX-1] , cmap='cool', alpha=0.8, edgecolors='none', s=50)
            #ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, c=Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50)
            ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, c='red', cmap='cool', alpha=0.8, edgecolors='none', s=50)
            ax.set_title('Well Injection Rate')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, BXL/1000])
            ax.set_ylim([0, BYL/1000])
            ax.set_xlabel('X , km')
            ax.set_ylabel('Y , km')
            ax.set_facecolor('w')
            cbar = plt.colorbar(p)
            cbar.set_label('Flow Rate (ton/day)')
        elif UnitSys == 1: # Field units
            #Qsave2[0:nWX*nWX,nWX-1] = Qsave2[0:nWX*nWX,nWX-1] * 0.000365 #Injection well convert from ton/day to MMT/year
            #Qsave2[nWX*nWX:nWX*nWX+nWE,nWX-1] = Qsave2[nWX*nWX:nWX*nWX+nWE,nWX-1] * 6.2898107 # Extraction well convert from m^3/day to bbl/day

            #vmin = min(np.min(Qsave2[0:nWX*nWX, nWX-1]), np.min(Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
            #vmax = max(np.max(Qsave2[0:nWX*nWX, nWX-1]), np.max(Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))     
            p = ax.scatter(Xwell[0:nWX*nWX, nWX-1]/1609.34, Ywell[0:nWX*nWX, nWX-1]/1609.34, c=Qsave2[0:nWX*nWX, nWX-1] * 0.000365, cmap='cool', alpha=0.8, edgecolors='none', s=50)
            #ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, c=Qsave2[nWX*nWX:nWX*nWX+nWE, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50)
            ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, c='red', cmap='cool', alpha=0.8, edgecolors='none', s=50)
            ax.set_title('Well Rate')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, BXL/1609.34])
            ax.set_ylim([0, BYL/1609.34])
            ax.set_xlabel('X , mile')
            ax.set_ylabel('Y , mile')
            ax.set_facecolor('w')
            cbar = plt.colorbar(p)
            cbar.set_label('Flow Rate (MMT/year)')
        else: 
            raise Exception ("Wrong Unit System selection!")
        
    else:
        if UnitSys == 2: # SI units
            if nWE > 0:
                vmin = min(np.min(Psave2[0:nWX*nWX, nWX-1]), np.min(Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
                vmax = max(np.max(Psave2[0:nWX*nWX, nWX-1]), np.max(Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
            else:
                vmin = np.min(Psave2[0:nWX*nWX, nWX-1])
                vmax = np.max(Psave2[0:nWX*nWX, nWX-1])
            p = ax.scatter(Xwell[0:nWX*nWX, nWX-1]/1000, Ywell[0:nWX*nWX, nWX-1]/1000, c=Psave2[0:nWX*nWX, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50, vmin=vmin, vmax=vmax)
            ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1000, c=Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50, vmin=vmin, vmax=vmax)
            ax.set_title('Well Pressure')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, BXL/1000])
            ax.set_ylim([0, BYL/1000])
            ax.set_xlabel('X , km')
            ax.set_ylabel('Y , km')
            ax.set_facecolor('w')
            cbar = plt.colorbar(p)
            cbar.set_label('Well Pressure (MPa)')
        elif UnitSys == 1: # Field Units
            Psave2 = Psave2 * 145.038 # Convert from MPa to psi
            if nWE > 0:
                vmin = min(np.min(Psave2[0:nWX*nWX, nWX-1]), np.min(Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
                vmax = max(np.max(Psave2[0:nWX*nWX, nWX-1]), np.max(Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1]))
            else:
                vmin = np.min(Psave2[0:nWX*nWX, nWX-1])
                vmax = np.max(Psave2[0:nWX*nWX, nWX-1])    
            p = ax.scatter(Xwell[0:nWX*nWX, nWX-1]/1609.34, Ywell[0:nWX*nWX, nWX-1]/1609.34, c=Psave2[0:nWX*nWX, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50, vmin=vmin, vmax=vmax)
            ax.scatter(Xwell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, Ywell[nWX*nWX:nWX*nWX+nWE, nWX-1]/1609.34, c=Psave2[nWX*nWX:nWX*nWX+nWE, nWX-1], cmap='cool', alpha=0.8, edgecolors='none', s=50, vmin=vmin, vmax=vmax)
            ax.set_title('Well Pressure')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, BXL/1609.34])
            ax.set_ylim([0, BYL/1609.34])
            ax.set_xlabel('X , mile')
            ax.set_ylabel('Y , mile')
            ax.set_facecolor('w')
            cbar = plt.colorbar(p)
            cbar.set_label('Well Pressure (Psi)')
        else: 
            raise Exception ("Wrong Unit System selection!")
    # Draw rectangles
    rect1 = plt.Rectangle(((rsvrprop.BXL-rsvrprop.XL)/2/lencoef, (rsvrprop.BYL-rsvrprop.YL)/2/lencoef), rsvrprop.XL/lencoef, rsvrprop.YL/lencoef, linewidth=2, linestyle='--', edgecolor='g', facecolor='none', label='Project Area')
    ax.add_patch(rect1)

    rect2 = plt.Rectangle((0, 0), rsvrprop.BXL/lencoef, rsvrprop.BYL/lencoef, linewidth=2, linestyle='--', edgecolor='b', facecolor='none', label='Reservoir Area')
    ax.add_patch(rect2)
    plt.legend(loc='upper right')
    ax.set_xlim([-rsvrprop.BXL/lencoef*0.05, rsvrprop.BXL/lencoef*1.05])
    ax.set_ylim([-rsvrprop.BYL/lencoef*0.05, rsvrprop.BYL/lencoef*1.05])
    fig.set_facecolor('w')
    ax.set_aspect('equal')
        
    return fig


def plot_AOR(Xwell, Ywell, Psave2, nWX, rsvrprop, simprop, UnitSys):
 
    if UnitSys == 2:
        lencoef = 1000
        pcoef = 1
    elif UnitSys == 1:
        lencoef = 1609.34
        pcoef = 145.038 # convert from MPa to psi
    else: 
        raise Exception ("Wrong Unit System selection!")
        
   
    # Prepare contour data for nodes
    XContour = [x / lencoef for x in Xwell[0:nWX**2+simprop.nWE,nWX-1]]
    YContour = [y / lencoef for y in Ywell[0:nWX**2+simprop.nWE,nWX-1]]
    PContour = [p * pcoef for p in Psave2[0:nWX**2+simprop.nWE,nWX-1]]

   
    
    #Interpolate data onto a grid
    x1 = np.arange(0, rsvrprop.BXL/lencoef, 1/200*rsvrprop.BXL/lencoef)
    x2 = np.arange(0, rsvrprop.BYL/lencoef, 1/200*rsvrprop.BYL/lencoef)
    X1, X2 = np.meshgrid(x1, x2)
    #PContour_flatten = flatten_list(PContour) # Make PContour in a list format
    
    highlight_AOR = rsvrprop.P0 + simprop.CP  # Specify the contour level to highlight
    highlight_AOR_unit = highlight_AOR * pcoef
    
    # Capture warning messages
    #AOR_warn = 0
    #if min(PContour) > highlight_AOR_unit or max(PContour) < highlight_AOR_unit:
    #    AOR_warn = 1
    #    warnings.warn('Critical Pressure is not suitable! AOR cannot be visualized!')
    foo = griddata((XContour, YContour), PContour, (X1, X2), method='linear')
   
    
    # Now foo is the interpolated data on the grid defined by X1 and X2
    X3 = foo  # You can use X3 for further processing
    
    # Create a contour plot with 9 levels
    fig_contour, ax = plt.subplots()
    contour = ax.contourf(X1, X2, X3, 9)
    
        # Set the x-axis and y-axis limits
    ax.set_xlim([0, rsvrprop.BXL/lencoef])
    ax.set_ylim([0, rsvrprop.BYL/lencoef])
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # Raise UserWarning as an error

    try:
        highlight_contour = ax.contour(X1, X2, X3, levels=[highlight_AOR_unit], colors='green', linewidths=2)
    except UserWarning:
        raise Exception("AOR cannot be visualized in the domain area.")
        st.warning("AOR cannot be visualized in the domain area.")

    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(contour)
    if UnitSys == 2:
        plt.title('Pressure Contour, MPa')
        plt.xlabel('X , km')
        plt.ylabel('Y , km')
        cbar.set_label('Pressure, MPa')
    else:
        plt.title('Pressure Contour, psi')
        plt.xlabel('X , mile')
        plt.ylabel('Y , mile')
        cbar.set_label('Pressure, psi')
        
    plt.axis('equal')
    
    # Draw rectangles
    rect1 = plt.Rectangle(((rsvrprop.BXL-rsvrprop.XL)/2/lencoef, (rsvrprop.BYL-rsvrprop.YL)/2/lencoef), rsvrprop.XL/lencoef, rsvrprop.YL/lencoef, linewidth=2, linestyle='--', edgecolor='g', facecolor='none', label='Project Area')
    ax.add_patch(rect1)

    rect2 = plt.Rectangle((0, 0), rsvrprop.BXL/lencoef, rsvrprop.BYL/lencoef, linewidth=2, linestyle='--', edgecolor='b', facecolor='none', label='Reservoir Area')
    ax.add_patch(rect2)
    plt.legend(loc='upper right')
    ax.set_xlim([-rsvrprop.BXL/lencoef*0.05, rsvrprop.BXL/lencoef*1.05])
    ax.set_ylim([-rsvrprop.BYL/lencoef*0.05, rsvrprop.BYL/lencoef*1.05])
    fig_contour.set_facecolor('w')
    ax.set_aspect('equal')
    return fig_contour, X3


def S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2, Qsave_res, ConRate):
    nWX = 50
    nWY = 50
    
    if UnitSys == 2:
        lencoef = 1000
        pcoef = 1
    elif UnitSys == 1:
        lencoef = 1609.34
        pcoef = 145.038 # convert from MPa to psi
    else: 
        raise Exception ("Wrong Unit System selection!")
        
    if BCValue == 1:
        rED1 = simprop.rE / simprop.rW  # closed
        rED2 = simprop.rE / simprop.rW  # closed
    elif BCValue == 2:
        rED1 = 1e38  # open
        rED2 = 1e38  # open
    else:
        raise ValueError('Error: invalid BCValue.') # or use an appropriate exception type and error message
    
    rW = simprop.rW
    nNodes = nWX*nWY   
    
    nWI = nWX_ctrl**2
    nWE = simprop.nWE
    
    Pnode = np.zeros(nNodes,order = 'F')
    chi_BL = np.zeros(nWI,order = 'F')
    chi_dry = np.zeros(nWI,order = 'F')
    Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
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
            Rwell_comp[i, j] = np.sqrt((X[i] - Xwell[j,nWX_ctrl-1]) ** 2 + (Y[i] - Ywell[j, nWX_ctrl-1]) ** 2)
    Rwell = np.real(Rwell_comp) # In case complex number occurs
            
    
    Sa = np.zeros(nWI)
    
    Fluid = Fluid_Prop(rsvrprop.P0, rsvrprop.temp, rsvrprop.salinity, relakprop.Sar, relakprop.Sgc,relakprop.m, relakprop.n, relakprop.kra0, relakprop.krg0)

    Fluid.Spycher()

    mug, cg, rhog, rhoc = Fluid.CO2Prop()

    cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()

    tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(rsvrprop.k,rsvrprop.Porosity,rsvrprop.cr,simprop.SimTime,simprop.rW)
    
    if ConRate == 1:
        Qinj=simprop.InjRate/86400/rhoc*1000 # rm^3/s # in reservoir condition
        Qext=simprop.ExtRate/86400/rhob*rhobstd # rm^3/s # in reservoir condition
    
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
                BInj += A[i,j] * Qinj 
            BExt = 0
            for j in range(nWI, nWI + nWE):
                BExt += A[i,j] * Qext
            Pnode[i] = mug / (2 * math.pi * rsvrprop.Thickness * rsvrprop.k * 
                              relakprop.krg0) * BInj / 1000000 - mub / (2 * math.pi * rsvrprop.Thickness 
                            * rsvrprop.k * relakprop.kra0) * BExt / 1000000 + rsvrprop.P0 # Krg0 to Krs
    elif ConRate == 0:
        Qinj = [0] * nWI
        Qext = [0] * simprop.nWE 
        for w in range(nWI):
             Qinj[w] = Qsave_res[w,nWX_ctrl-1] #rm^3/s # in reservoir condition    
             epsilon=np.multiply(Qinj[w],(cg+rsvrprop.cr)*mug/(4*math.pi*rsvrprop.Thickness*rsvrprop.k*relakprop.krg0))
             chi_BL[w] = (1/4)*epsilon*dfgdSgL
             chi_dry[w] = (1/4)*epsilon*dfgdSgT
        for w in range(nWE):
             Qext[w] = -1 * Qsave_res[w + nWI, nWX_ctrl-1] # rm^3/s # in reservoir condition
        
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
                BExt += A[i,j] * Qext[j-nWI]
            Pnode[i] = mug / (2 * math.pi * rsvrprop.Thickness * rsvrprop.k * 
                              relakprop.krg0) * BInj / 1000000 - mub / (2 * math.pi * rsvrprop.Thickness 
                            * rsvrprop.k * relakprop.kra0) * BExt / 1000000 + rsvrprop.P0 # Krg0 to Krs
    else: 
        raise ValueError("Inconrrect Flow Rate Condition Selection.")
    
                                                                    
    # Prepare contour data for nodes
    XContour = np.array(X)/lencoef
    YContour = np.array(Y)/lencoef
    PContour = np.array(Pnode)
    
    # Prepare contour data for wells
    for i in range(nWI + nWE):
        XContour = np.append(XContour, Xwell[i, nWX_ctrl-1]/lencoef)
        YContour = np.append(YContour, Ywell[i, nWX_ctrl-1]/lencoef)    
        PContour = np.append(PContour, Psave2[i, nWX_ctrl-1])
    
    # Interpolate data onto a grid
    x1 = np.linspace(0, rsvrprop.BXL/lencoef, 200)
    x2 = np.linspace(0, rsvrprop.BYL/lencoef, 200)
    
    X1, X2 = np.meshgrid(x1, x2)
    
    PContour_flatten_unit = np.array(PContour)  # No need for flatten_list() conversion
    
    highlight_AOR = rsvrprop.P0 + simprop.CP
    
    # UNIT CONVERSION
    if UnitSys == 1:  # field unit, then output pressure is in psi
        PContour_flatten_unit *= 145.03768  # Convert from MPa to psi
        Psave2_output = Psave2 * 145.03768  # Convert from MPa to psi
        highlight_AOR_unit = highlight_AOR * 145.03768  # Convert from MPa to psi
    else:  # SI unit
        Psave2_output = Psave2  # Output pressure still in MPa
        highlight_AOR_unit = highlight_AOR  # in MPa
   

    # Interpolate data onto a grid
    foo = griddata((XContour, YContour), PContour_flatten_unit, (X1, X2), method='linear')
    
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
        st.warning("Warning: Cannot plot AOR in the pressure contour map, either AOR it too large or too small.")

    plt.gca().set_aspect('equal', adjustable='box')
    if UnitSys == 2:
        plt.title('Pressure Contour, MPa')
        plt.xlabel('X , km')
        plt.ylabel('Y , km')
        cbar.set_label('Pressure, MPa')
    else:
        plt.title('Pressure Contour, psi')
        plt.xlabel('X , mile')
        plt.ylabel('Y , mile')
        cbar.set_label('Pressure, psi')
        
    plt.axis('equal')
    
    # Draw rectangles
    rect1 = plt.Rectangle(((rsvrprop.BXL-rsvrprop.XL)/2/lencoef, (rsvrprop.BYL-rsvrprop.YL)/2/lencoef), rsvrprop.XL/lencoef, rsvrprop.YL/lencoef, linewidth=2, linestyle='--', edgecolor='g', facecolor='none', label='Project Area')
    ax.add_patch(rect1)
 
    rect2 = plt.Rectangle((0, 0), rsvrprop.BXL/lencoef, rsvrprop.BYL/lencoef, linewidth=2, linestyle='--', edgecolor='b', facecolor='none', label='Reservoir Area')
    ax.add_patch(rect2)
    plt.legend(loc='upper right')
    ax.set_xlim([-rsvrprop.BXL/lencoef*0.05, rsvrprop.BXL/lencoef*1.05])
    ax.set_ylim([-rsvrprop.BYL/lencoef*0.05, rsvrprop.BYL/lencoef*1.05])
    fig_contour.set_facecolor('w')
    ax.set_aspect('equal')
    

    # Draw scatter points
    # Create a set to store unique labels
    legend_labels = set()
    for w in range(nWX_ctrl**2):
        x = Xwell[w, nWX_ctrl-1] / lencoef
        y = Ywell[w, nWX_ctrl-1] / lencoef
        ax.scatter(x, y, color='r', marker='o', s=20)  # Adjust the 's' parameter as desired
        legend_labels.add('Injectors')
    if nWE > 0: 
        for w in range(nWE):
            x = Xwell[w + nWX_ctrl**2, nWX_ctrl-1] / lencoef
            y = Ywell[w + nWX_ctrl**2, nWX_ctrl-1] / lencoef
            ax.scatter(x, y, color='b', marker='^', s=20)  # Adjust the 's' parameter as desired
            legend_labels.add('Extractors')
    ax.legend()
    return highlight_locations, fig_contour, Pnode, highlight_contour
    
'''
#### TESTING
# Read input properties from the uploaded Excel file, for testing purposes
from properties import read_rsvr_input, Unit_conversion, userexcel_conversion
import pandas as pd
from S2 import Scenario2
from S1 import Scenario1
from Matcal import Distance_Mat_Builder

sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_SI.xlsx'
UnitSys = 2
BCValue = 2

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
input_df, proj_name = userexcel_conversion(userinput_df)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)

Psave2, Qsave2, rL1, Capacity, NPV = Scenario1(rsvrprop, simprop, relakprop, npvprop, BCValue)
         

nWXmax = simprop.nWXmax 
nWXmin = simprop.nWXmin
nWE = simprop.nWE
XL = rsvrprop.XL
YL = rsvrprop.YL
BXL = rsvrprop.BXL
BYL = rsvrprop.BYL
ConRate = 1  


Xwell, Ywell, Rwell = Distance_Mat_Builder(nWXmax, nWE, XL, YL, BXL, BYL)
nWX_ctrl = 9

#fig1 = plot_plume_extension(nWX_ctrl, Xwell, Ywell, rL1, nWE, BXL, XL, BYL, YL, UnitSys)
#fig2 = plot_well_property(Xwell, Ywell, nWX_ctrl,nWE, ConRate, BXL, BYL, Qsave2, Psave2, UnitSys)
# Clear the figure
#plt.clf()
#plot_capacity_vs_NumInj(NumInj, Capacity, WarningCountPlume, WarningCountPressure, nWXmin, nWXmax)
#plot_npv_vs_NumInj(NumInj, NPV, nWXmin, nWXmax)
WarningCountPlume, WarningPlume = calculate_warning_plume(nWXmax, rL1, Xwell, Ywell, Rwell, BXL, BYL, XL, YL)
WarningCountPressure, WarningPressure = calculate_warning_pressure(nWXmin, nWXmax, ConRate, Psave2, rsvrprop, Qsave2, nWE, UnitSys)

#fig3, X3 =plot_AOR(Xwell, Ywell, Psave2, nWX_ctrl, rsvrprop, simprop, UnitSys)
'''