"""
Created on Sun May 21 17:40:40 2023

@author: wangz
"""
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
#from tqdm import tqdm
from S2 import Scenario2
from S3 import S3_core
import streamlit as st
import matplotlib.patches as mpatches
from plotting import S1S2_PContour
from Matcal import Distance_Mat_Builder
import pandas as pd

class Inputclass:
    def __init__(self, real_num, means, sds):
        self.real_num = real_num
        self.attributes = ['P0', 'temp', 'Thickness', 'k', 'salinity', 'm','n', 'cr', 'FracP', 'kra0', 'krg0', 'Porosity','Sar','Sgc']
        self.means = means
        self.sds = sds

        self.generate_values()


    def generate_values(self):
        for attr, mean, sd in zip(self.attributes, self.means, self.sds):
            setattr(self, f'{attr}_values', np.random.normal(mean, sd, self.real_num))
            setattr(self, f'{attr}_mean', mean)
            setattr(self, f'{attr}_sd', sd)



def sensi_ana(self, Sensi, rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, mins, maxes, geometry):
    self.real_num = self.real_num

    attributes = ['P0', 'temp', 'Thickness', 'k', 'salinity', 'm', 'n', 'cr', 'FracP', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc']
    means = [self.P0_mean, self.temp_mean, self.Thickness_mean, self.k_mean, self.salinity_mean, self.m_mean, self.n_mean, 
             self.cr_mean, self.FracP_mean, self.kra0_mean, self.krg0_mean, self.Porosity_mean, self.Sar_mean, self.Sgc_mean]
    sds = [self.P0_sd, self.temp_sd, self.Thickness_sd, self.k_sd, self.salinity_sd, self.m_sd, self.n_sd, self.cr_sd, 
           self.FracP_sd, self.kra0_sd, self.krg0_sd, self.Porosity_sd, self.Sar_sd, self.Sgc_sd]
    
    # Random sample only the selected parameter, others para are taken as mean values:
    #for attr, mean, sd in zip(attributes, means, sds):
    #    if attr == Sensi_para:
    #        setattr(self, f'{attr}_values', np.random.normal(mean, sd, self.real_num))
    #    else:
    #        setattr(self, f'{attr}_values', np.full(self.real_num, mean))
    
    # Random sample all parameters based on sd
    #for attr, mean, sd in zip(attributes, means, sds):
    #    setattr(self, f'{attr}_values', np.random.normal(mean, sd, self.real_num))
        
    # Random sample from min and max
    for attr, min_val, max_val in zip(attributes, mins, maxes):
        setattr(self, f'{attr}_values', np.random.uniform(min_val, max_val, self.real_num))
    
    SA_values = np.zeros(self.real_num, order='F')
    NPV_values = np.zeros(self.real_num, order='F')
    
   
    
    ## Create a progress bar
    #progress_bar = tqdm(total=self.real_num, bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')
    
    for i in range(self.real_num):
        
        ## Update the progress bar
        #progress_bar.update(1)
        
        rsvrprop.P0 = self.P0_values[i]
        rsvrprop.temp = self.temp_values[i]
        rsvrprop.Thickness = self.Thickness_values[i]
        rsvrprop.k = self.k_values[i]
        rsvrprop.salinity = self.salinity_values[i]
        rsvrprop.FracP = self.FracP_values[i]
        rsvrprop.Porosity = self.Porosity_values[i]
        rsvrprop.cr = self.cr_values[i]
        relakprop.m = self.m_values[i]
        relakprop.n = self.n_values[i]
        relakprop.Sar = self.Sar_values[i]
        relakprop.Sgc = self.Sgc_values[i]
        relakprop.kra0 = self.kra0_values[i]
        relakprop.krg0 = self.krg0_values[i]
        

        if geometry == 0: 
            SA_values[i], NPV_values[i] = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi) # Here SA_values stors Capacity esitmation
        elif geometry == 1:
            Sensi = 0
            SA_values[i], NPV_values[i], _, _,_,_,_ = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi) # Here SA_values stors P_BHP ave for inj esitmation
            if UnitSys == 1:
                SA_values[i] = SA_values[i] * 145.038 # As output is average Inj. BHP pressure
                
        else:
            raise Exception('Wrong Geometry Function Selection!')


    return SA_values, NPV_values
    
    


def plot_density_curve(hist_data, geometry, UnitSys):
    # Create the KDE object using the histogram data
    fig_hist, ax = plt.subplots()
    
    hist_data = hist_data[~np.isnan(hist_data) & np.isfinite(hist_data)] #remove nan and inf
    kde = gaussian_kde(hist_data)

    # Generate the smooth curve by evaluating the KDE over a range of values
    x_vals = np.linspace(hist_data.min(), hist_data.max(), 1000)
    smooth_curve = kde.evaluate(x_vals)

    # Calculate p5, p50 and p95 values
    p90 = np.percentile(hist_data, 10)  # 5th percentile
    p10 = np.percentile(hist_data, 90)  # 95th percentile
    p50 = np.percentile(hist_data,50) # 50th percentile, median
    #mean_value = np.mean(hist_data)  # Mean

    # Plot the histogram and the smooth curve
    n, bins, patches = plt.hist(hist_data, bins='auto', density=True, alpha=0.5)


    #plt.plot(x_vals, smooth_curve, 'r-', label = 'smooth curve' )
    plt.plot(x_vals, smooth_curve, 'r-' )
    plt.axvline(p10, color='blue', linestyle='--', label=f'P10 = {p10:.2f}')
    plt.axvline(p50, color='red', linestyle='--', label=f'P50 = {p50:.2f}')
    plt.axvline(p90, color='green', linestyle='--', label=f'P90 = {p90:.2f}')
    plt.legend()
    
    if geometry == 0:
        plt.xlabel('Capacity, MMT')
    elif geometry ==1:
        if UnitSys == 1:
            plt.xlabel('Avg. BHP of Injectors, psi')
        else: 
            plt.xlabel('Avg. BHP of Injectors, MPa')
    
    plt.ylabel('Histogram Density')
    plt.title('Monte Carlo Simulation for Sensitivity Analysis')

    # Annotate 90% confidence range on the graph
    
    #plt.annotate(f'{p10:.2f}', xy=(p10, 0), xytext=(p10-0.05, 0.01), ha='center', color='blue')
    #plt.annotate(f'{p90:.2f}', xy=(p90, 0), xytext=(p90-0.05, 0.01), ha='center', color='green')
    #plt.annotate(f'{p50:.2f}', xy=(p50, 0), xytext=(p50-0.05, 0.01), ha='center', color='red')
  
    #plt.annotate(f'{mean_value:.2f}', xy=(mean_value, 0), xytext=(mean_value, 0.01), ha='center', color='yellow')


    return fig_hist


def generate_histdata(Sensi_para_rel_num, rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, geometry, output_print):
    
    means = [rsvrprop.P0, rsvrprop.temp, rsvrprop.Thickness, rsvrprop.k, rsvrprop.salinity, relakprop.m, relakprop.n, rsvrprop.cr, rsvrprop.FracP,
            relakprop.kra0, relakprop.krg0, rsvrprop.Porosity, relakprop.Sar, relakprop.Sgc]
    sds = [mean * 0.01 for mean in means]
    
    # Samples number
    real_num = Sensi_para_rel_num
    
    SA_input_data, max_values_df, min_values_df = read_SA_input(rsvrprop, relakprop, output_print, UnitSys, uploaded_file)
    attributes = ['P0', 'temp', 'Thickness', 'k', 'salinity', 'm', 'n', 'cr', 'FracP', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc']
    act_attributes = ['Pressure', 'Temperature', 'Thickness', 'Permeability', 'Salinity', 'm', 'n', 'Rock Comp.', 'Max Inj. Pressure', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc']
    
    attribute_mapping = {
    'Pressure': 'P0',
    'Temperature': 'temp',
    'Thickness': 'Thickness',
    'Permeability': 'k',
    'Salinity': 'salinity',
    'Rock Comp.': 'cr',
    'Max Inj. Pressure': 'FracP',
    'kra0': 'kra0',
    'krg0': 'krg0',
    'Porosity':  'Porosity',
    'Sar': 'Sar',
    'Sgc': 'Sgc',
    'm':'m',
    'n':'n'}

    maxes = []

    for attr in attributes:
        max_value = max_values_df[max_values_df['Parameter'] == attr]['Max'].values[0]
        maxes.append(max_value)
        
    mins = []

    for attr in attributes:
        min_value = min_values_df[min_values_df['Parameter'] == attr]['Min'].values[0]
        mins.append(min_value)
        
    
    input_data = Inputclass(real_num, means, sds)
    Sensi = 1
    SA_values, NPV_values = sensi_ana(input_data, Sensi, rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, mins, maxes, geometry)
    hist_data = SA_values
    return hist_data




def sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop, BCValue, output_print, UnitSys, geometry, uploaded_file, proj_name):
    ''' 
    Tornado chart to run sensitivity analysis
    
    
    '''
    SA_input_data, max_values_df, min_values_df = read_SA_input(rsvrprop, relakprop, output_print, UnitSys, uploaded_file)
    
    means = [rsvrprop.P0, rsvrprop.temp, rsvrprop.Thickness, rsvrprop.k, rsvrprop.salinity, relakprop.m, relakprop.n, rsvrprop.cr, rsvrprop.FracP,
             relakprop.kra0, relakprop.krg0, rsvrprop.Porosity, relakprop.Sar, relakprop.Sgc]
    
    sds = [mean * 0.01 for mean in means] # Meaningless here, just for input
     
    # Samples number
    real_num = 1 # Meaning less here
     
    input_data = Inputclass(real_num, means, sds)
     
    attributes = ['P0', 'temp', 'Thickness', 'k', 'salinity', 'm', 'n', 'cr', 'FracP', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc']
    act_attributes = ['Pressure', 'Temperature', 'Thickness', 'Permeability', 'Salinity', 'm', 'n', 'Rock Comp.', 'Max Inj. Pressure', 'kra0', 'krg0', 'Porosity', 'Sar', 'Sgc']
    
    attribute_mapping = {
    'Pressure': 'P0',
    'Temperature': 'temp',
    'Thickness': 'Thickness',
    'Permeability': 'k',
    'Salinity': 'salinity',
    'Rock Comp.': 'cr',
    'Max Inj. Pressure': 'FracP',
    'kra0': 'kra0',
    'krg0': 'krg0',
    'Porosity':  'Porosity',
    'Sar': 'Sar',
    'Sgc': 'Sgc',
    'm':'m',
    'n':'n'}

    maxes = []

    for attr in attributes:
        max_value = max_values_df[max_values_df['Parameter'] == attr]['Max'].values[0]
        maxes.append(max_value)
        
    mins = []

    for attr in attributes:
        min_value = min_values_df[min_values_df['Parameter'] == attr]['Min'].values[0]
        mins.append(min_value)
    
    
    means = [input_data.P0_mean, input_data.temp_mean, input_data.Thickness_mean, input_data.k_mean, input_data.salinity_mean, input_data.m_mean, input_data.n_mean, 
             input_data.cr_mean, input_data.FracP_mean, input_data.kra0_mean, input_data.krg0_mean, input_data.Porosity_mean, input_data.Sar_mean, input_data.Sgc_mean]
    sds = [input_data.P0_sd, input_data.temp_sd, input_data.Thickness_sd, input_data.k_sd, input_data.salinity_sd, input_data.m_sd, input_data.n_sd, input_data.cr_sd, 
           input_data.FracP_sd, input_data.kra0_sd, input_data.krg0_sd, input_data.Porosity_sd, input_data.Sar_sd, input_data.Sgc_sd]
    

    count = len(attributes) * 3
    #attr_values = [[] for _ in range(count)]
    
    #for i, (attr, mean, sd) in enumerate(zip(attributes, means, sds)):
    #    values_sd_mean = [mean - 2 * sd] + [mean] + [mean + 2 * sd] + [mean] * (count - 3)
   
    #    for j in range(count):
    #        attr_values[j].append(values_sd_mean[(j - 3*i) % (count)]) # shift element by 3*i when enumerating
    
    
    attr_values = [[] for _ in range(count)]
    for i, (attr, mean) in enumerate(zip(attributes, means)):
        values_sd_mean = [mins[i]] + [mean] + [maxes[i]] + [mean] * (count - 3)
        for j in range(count):
            attr_values[j].append(values_sd_mean[(j - 3*i) % (count)])  # shift element by 3*i when enumerating

    # Print attr_values
    #for attr, values in zip(attributes, attr_values):
    #    print(f"{attr}_values: {values}")
    Sensi = 1
        
    SA_values_tor = np.zeros(count, order='F')
    NPV_values_tor = np.zeros(count, order='F')
        
    for i in range(count):
        rsvrprop.P0 = attr_values[i][attributes.index('P0')]
        rsvrprop.temp = attr_values[i][attributes.index('temp')]
        rsvrprop.Thickness = attr_values[i][attributes.index('Thickness')]
        rsvrprop.k = attr_values[i][attributes.index('k')]
        rsvrprop.salinity = attr_values[i][attributes.index('salinity')]
        rsvrprop.FracP = attr_values[i][attributes.index('FracP')]
        rsvrprop.Porosity = attr_values[i][attributes.index('Porosity')]
        rsvrprop.cr = attr_values[i][attributes.index('cr')]
        relakprop.m = attr_values[i][attributes.index('m')]
        relakprop.n = attr_values[i][attributes.index('n')]
        relakprop.Sar = attr_values[i][attributes.index('Sar')]
        relakprop.Sgc = attr_values[i][attributes.index('Sgc')]
        relakprop.kra0 = attr_values[i][attributes.index('kra0')]
        relakprop.krg0 = attr_values[i][attributes.index('krg0')]
        
        if geometry == 0: 
            SA_values_tor[i], NPV_values_tor[i] = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi) # Here SA_values stors Capacity esitmation
        elif geometry == 1:
            Sensi = 0
            SA_values_tor[i], NPV_values_tor[i], _, _, _,_,_ = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi) # Here SA_values stors P_BHP ave for inj esitmation
            if UnitSys == 1:
                SA_values_tor[i] = SA_values_tor[i] * 145.038
                
        else:
            raise Exception('Wrong Geometry Function Selection!')    
        #Capacity_values_tor[i], NPV_values_tor[i] = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
    
    mean_values = []
    low_limit_values = []
    up_limit_values = []
    
    for i in range(len(attributes)):
        index = i * 3  # Starting index for the attribute
        mean_values.append(SA_values_tor[index + 1]) # means 
        low_limit_values.append(SA_values_tor[index]) # lower limits 
        up_limit_values.append(SA_values_tor[index + 2]) # upper limits
    
    difference_values = np.array(up_limit_values) - np.array(low_limit_values)
    sorted_indices = np.argsort(difference_values)[::-1]
    sorted_attributes = [act_attributes[i] for i in sorted_indices]
    sorted_mean_values = [mean_values[i] for i in sorted_indices]
    #sorted_sd_values = [sd_values[i] for i in sorted_indices]
    sorted_low_limit_values = [low_limit_values[i] for i in sorted_indices]
    sorted_up_limit_values = [up_limit_values[i] for i in sorted_indices]

    # Tornado Chart
    fig_tor, ax = plt.subplots()
    positions = np.arange(len(sorted_attributes))
    width = 0.5
    
    # Plotting the bars for mean values
    for sorted_mean_value in sorted_mean_values:
        ax.axvline(x=sorted_mean_value, color='gray', linestyle='--', linewidth = 0.5)

    
    # Plotting the bars for high and low limits
    # Create a list of colors based on the condition
    colors = ['gray' if up < low else 'blue' for up, low in zip(sorted_up_limit_values, sorted_low_limit_values)]
    posi_rela = [ -1 if up < low else 1 for up, low in zip(sorted_up_limit_values, sorted_low_limit_values)]
    #Positive relationship between attributes to capacity.
 
    
    # Create the bar plot
    ax.barh(positions, np.subtract(sorted_up_limit_values, sorted_low_limit_values),
        left=sorted_low_limit_values, height=width, color=colors, alpha=0.3)
    
    # Create legend patches
    negative_patch = mpatches.Patch(color='gray', label='Negative Correlation')
    positive_patch = mpatches.Patch(color='blue', label='Positive Correlation')
    
    # Add the legend
    ax.legend(handles=[negative_patch, positive_patch])
    ax.set_yticks(positions)
    ax.set_yticklabels(sorted_attributes)
    if geometry == 0:
        plt.xlabel('Capacity, MMT')
    elif geometry ==1:
        if UnitSys == 1:
            plt.xlabel('Avg. BHP of Injectors, psi')
        else: 
            plt.xlabel('Avg. BHP of Injectors, MPa')
    ax.set_title('Sensitivity Tornado Chart')
    #ax.legend()

    plt.tight_layout()
    
    if geometry == 0:
    
        with open('EASiToolOutput_Sensitivity.txt', 'w') as file:        
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
                       
                
            file.write('\n')
            file.write('*********************************************************************\n') 
            file.write('%43s\n' % 'EASiTool Output Tornado Chart Data')
            file.write('*********************************************************************\n') 
            
            file.write('{:<30}\n'.format('Mean Values Based on Input Parameters'))
            file.write('{:<20}{:<20}{:<20}\n'.format('Mean Value','Capacity (MMT):', 'NPV ($M):'))
            file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ',SA_values_tor[1], NPV_values_tor[1]))
            file.write('{:<75}\n'.format('*' * 69))    
            
            for i in range(len(attributes)):
                index = i * 3  # Starting index for the attribute
                file.write('{:<20}\n'.format(act_attributes[i]))
                file.write('{:<20}{:<20}{:<20}\n'.format('Upper Limit', 'Capacity (MMT):', 'NPV ($M):'))
                file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index], NPV_values_tor[index]))
                file.write('{:<20}{:<20}{:<20}\n'.format('Lower Limit', 'Capacity (MMT):', 'NPV ($M):'))
                file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index + 2], NPV_values_tor[index + 2]))
                file.write('{:<75}\n'.format('*' * 69))
    elif geometry == 1:
        with open('EASiToolOutput_Sensitivity.txt', 'w') as file:        
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
                       
                
            file.write('\n')
            file.write('*********************************************************************\n') 
            file.write('%43s\n' % 'EASiTool Output Tornado Chart Data')
            file.write('*********************************************************************\n') 
            
            file.write('{:<30}\n'.format('Mean Values Based on Input Parameters'))
            if UnitSys == 1:
                
                for i in range(len(attributes)):
                    index = i * 3  # Starting index for the attribute
                    file.write('{:<20}\n'.format(act_attributes[i]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Upper Limit', 'Ave Inj BHP(psi):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index], NPV_values_tor[index]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Mean', 'Ave Inj BHP(psi):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index + 1], NPV_values_tor[index + 1]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Lower Limit', 'Ave Inj BHP(psi):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index + 2], NPV_values_tor[index + 2]))
                    file.write('{:<75}\n'.format('*' * 69))
            elif UnitSys == 2:
             
                
                for i in range(len(attributes)):
                    index = i * 3  # Starting index for the attribute
                    file.write('{:<20}\n'.format(act_attributes[i]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Upper Limit', 'Ave Inj BHP(MPa):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index], NPV_values_tor[index]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Mean', 'Ave Inj BHP(MPa):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index + 1], NPV_values_tor[index + 1]))
                    file.write('{:<20}{:<20}{:<20}\n'.format('Lower Limit', 'Ave Inj BHP(MPa):', 'NPV ($M):'))
                    file.write('{:<20}{:<20.3e}{:<20.3e}\n'.format(' ', SA_values_tor[index + 2], NPV_values_tor[index + 2]))
                    file.write('{:<75}\n'.format('*' * 69))

    def plot_multi_AOR():
        
        df_aor_low = []
        df_aor_up = []
        df_aor_mean = []
        X_center = 0
        Y_center = 0

        ### Mean ###################################################################
        if geometry ==0: 
            nWX_ctrl = simprop.nWX_npv_max
            ConRate = 0
            Sensi = 0
            Psave2, _, _, _, _, _, Qsave_res = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
            highlight_locations, _,_, highlight_contour = S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2, Qsave_res, ConRate)
            if 'highlight_locations' in locals() and highlight_locations:
                x_mean, y_mean = zip(*highlight_locations)
            else:
                x_mean, y_mean = [], []
            
        elif geometry == 1:
            Sensi = 1
            _,_, highlight_locations, X_center, Y_center,_,_ = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi) 
            if 'highlight_locations' in locals() and highlight_locations:
                x_mean, y_mean = zip(*highlight_locations)
                df_aor_mean = pd.DataFrame(highlight_locations, columns=['x', 'y'])
            else:
                x_mean, y_mean = [], []
                df_aor_mean = pd.DataFrame(columns=['x', 'y'])

       
        ### Up #######################################################################
        #rsvrprop.P0 = input_data.P0_mean + 2 *input_data.P0_sd * posi_rela[sorted_attributes.index('Pressure')]
        #rsvrprop.temp = input_data.temp_mean + 2 *input_data.temp_sd * posi_rela[sorted_attributes.index('Temperature')]
        #rsvrprop.Thickness = input_data.Thickness_mean + 2 *input_data.Thickness_sd * posi_rela[sorted_attributes.index('Thickness')]
        #rsvrprop.k = input_data.k_mean + 2 *input_data.k_sd * posi_rela[sorted_attributes.index('Permeability')]
        #rsvrprop.salinity = input_data.salinity_mean + 2 *input_data.salinity_sd * posi_rela[sorted_attributes.index('Salinity')]
        #rsvrprop.FracP = input_data.FracP_mean + 2 *input_data.FracP_sd * posi_rela[sorted_attributes.index('Max Inj. Pressure')]
        #rsvrprop.Porosity = input_data.Porosity_mean + 2 *input_data.Porosity_sd * posi_rela[sorted_attributes.index('Porosity')]
        #rsvrprop.cr = input_data.cr_mean + 2 *input_data.cr_sd * posi_rela[sorted_attributes.index('Rock Comp.')]
        #relakprop.m = input_data.m_mean + 2 *input_data.m_sd * posi_rela[sorted_attributes.index('m')]
        #relakprop.n = input_data.n_mean + 2 *input_data.n_sd * posi_rela[sorted_attributes.index('n')]
        #relakprop.Sar = input_data.Sar_mean + 2 *input_data.Sar_sd * posi_rela[sorted_attributes.index('Sar')]
        #$relakprop.Sgc = input_data.Sgc_mean + 2 *input_data.Sgc_sd * posi_rela[sorted_attributes.index('Sgc')]
        #relakprop.kra0 = input_data.kra0_mean + 2 *input_data.kra0_sd * posi_rela[sorted_attributes.index('kra0')]
        #relakprop.krg0 = input_data.krg0_mean + 2 *input_data.krg0_sd * posi_rela[sorted_attributes.index('krg0')]
        
        parameters_to_check_rsvr = ['Pressure', 'Temperature', 'Thickness', 'Permeability', 'Salinity', 'Max Inj. Pressure', 'Porosity', 'Rock Comp.']
        parameters_to_check_relak = ['m','n','Sar','Sgc','kra0','krg0']
        for parameter in parameters_to_check_rsvr:
            attr_name = attribute_mapping.get(parameter)            
            if attr_name:
                attr_index = sorted_attributes.index(parameter)
                if posi_rela[attr_index] > 0:
                    max_value = max_values_df[max_values_df['Parameter'] == attr_name]['Max'].values[0]
                    setattr(rsvrprop, attr_name, max_value)
                else:
                    min_value = min_values_df[min_values_df['Parameter'] == attr_name]['Min'].values[0]
                    setattr(rsvrprop, attr_name, min_value)
            else:
                print(f"No mapping found for parameter: {parameter}")
                
        for parameter in parameters_to_check_relak:
            attr_name = attribute_mapping.get(parameter)            
            if attr_name:
                attr_index = sorted_attributes.index(parameter)
                if posi_rela[attr_index] > 0:
                    max_value = max_values_df[max_values_df['Parameter'] == attr_name]['Max'].values[0]
                    setattr(relakprop, attr_name, max_value)
                else:
                    min_value = min_values_df[min_values_df['Parameter'] == attr_name]['Min'].values[0]
                    setattr(relakprop, attr_name, min_value)
            else:
                print(f"No mapping found for parameter: {parameter}")
        
        
        
        if geometry == 0:
            nWX_ctrl = simprop.nWX_npv_max
            ConRate = 0
            Sensi = 0
            Psave2_up, _, _, _, _, _, Qsave_res_up = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
            highlight_locations_up, _,_, highlight_contour_up = S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2_up, Qsave_res_up, ConRate)
            if 'highlight_locations_up' in locals() and highlight_locations_up:
                x_up, y_up = zip(*highlight_locations_up)
            else:
                x_up, y_up = [], []
                
        elif geometry == 1:
            Sensi = 1
            _,_, highlight_locations_up, _, _,_,_ = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi) # Here SA_values stors P_BHP ave for inj esitmation
            if 'highlight_locations_up' in locals() and highlight_locations_up:
                x_up, y_up = zip(*highlight_locations_up)
                df_aor_up = pd.DataFrame(highlight_locations_up, columns=['x', 'y'])
            else:
                x_up, y_up = [], []
                df_aor_up = pd.DataFrame(columns=['x', 'y'])
        
        ### low ###########################################################################################
        #rsvrprop.P0 = input_data.P0_mean - 2 *input_data.P0_sd * posi_rela[sorted_attributes.index('Pressure')]
        #rsvrprop.temp = input_data.temp_mean - 2 *input_data.temp_sd * posi_rela[sorted_attributes.index('Temperature')]
        #rsvrprop.Thickness = input_data.Thickness_mean - 2 *input_data.Thickness_sd * posi_rela[sorted_attributes.index('Thickness')]
        #rsvrprop.k = input_data.k_mean - 2 *input_data.k_sd * posi_rela[sorted_attributes.index('Permeability')]
        #rsvrprop.salinity = input_data.salinity_mean - 2 *input_data.salinity_sd * posi_rela[sorted_attributes.index('Salinity')]
        #rsvrprop.FracP = input_data.FracP_mean - 2 *input_data.FracP_sd * posi_rela[sorted_attributes.index('Max Inj. Pressure')]
        #rsvrprop.Porosity = input_data.Porosity_mean - 2 *input_data.Porosity_sd * posi_rela[sorted_attributes.index('Porosity')]
        #rsvrprop.cr = input_data.cr_mean - 2 *input_data.cr_sd * posi_rela[sorted_attributes.index('Rock Comp.')]
        #relakprop.m = input_data.m_mean - 2 *input_data.m_sd * posi_rela[sorted_attributes.index('m')]
        #relakprop.n = input_data.n_mean - 2 *input_data.n_sd * posi_rela[sorted_attributes.index('n')]
        #relakprop.Sar = input_data.Sar_mean - 2 *input_data.Sar_sd * posi_rela[sorted_attributes.index('Sar')]
        #relakprop.Sgc = input_data.Sgc_mean - 2 *input_data.Sgc_sd * posi_rela[sorted_attributes.index('Sgc')]
        #relakprop.kra0 = input_data.kra0_mean - 2 *input_data.kra0_sd * posi_rela[sorted_attributes.index('kra0')]
        #relakprop.krg0 = input_data.krg0_mean - 2 *input_data.krg0_sd * posi_rela[sorted_attributes.index('krg0')]
        
        for parameter in parameters_to_check_rsvr:
            attr_name = attribute_mapping.get(parameter)            
            if attr_name:
                attr_index = sorted_attributes.index(parameter)
                if posi_rela[attr_index] > 0:
                    min_value = min_values_df[min_values_df['Parameter'] == attr_name]['Min'].values[0]
                    setattr(rsvrprop, attr_name, min_value)
                else:
                    max_value = max_values_df[max_values_df['Parameter'] == attr_name]['Max'].values[0]
                    setattr(rsvrprop, attr_name, max_value)
            else:
                print(f"No mapping found for parameter: {parameter}")
                
        for parameter in parameters_to_check_relak:
            attr_name = attribute_mapping.get(parameter)            
            if attr_name:
                attr_index = sorted_attributes.index(parameter)
                if posi_rela[attr_index] > 0:
                    min_value = min_values_df[min_values_df['Parameter'] == attr_name]['Min'].values[0]
                    setattr(relakprop, attr_name, min_value)
                else:
                    max_value = max_values_df[max_values_df['Parameter'] == attr_name]['Max'].values[0]
                    setattr(relakprop, attr_name, max_value)
            else:
                print(f"No mapping found for parameter: {parameter}")


        if geometry == 0:
            nWX_ctrl = simprop.nWX_npv_max
            ConRate = 0
            Sensi = 0
            Psave2_low, _, _, _, _, _, Qsave_res_low = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
            highlight_locations_low, _, _, highlight_contour_low = S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2_low, Qsave_res_low, ConRate)
            if 'highlight_locations_low' in locals() and highlight_locations_low:
                x_low, y_low = zip(*highlight_locations_low)
            else:
                x_low, y_low = [], []
        elif geometry ==1: 
            Sensi = 1
            _,_, highlight_locations_low, _, _,_,_ = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi) # Here SA_values stors P_BHP ave for inj esitmation
            if 'highlight_locations_low' in locals() and highlight_locations_low:
                x_low, y_low = zip(*highlight_locations_low)
                df_aor_low = pd.DataFrame(highlight_locations_low, columns=['x', 'y'])
            else:
                x_low, y_low = [], []
                df_aor_low = pd.DataFrame(columns=['x', 'y'])

        
        if UnitSys == 2:
            lencoef = 1000
            pcoef = 1
        elif UnitSys == 1:
            lencoef = 1609.34
            pcoef = 145.038 # convert from MPa to psi
        else: 
            raise Exception ("Wrong Unit System selection!")
       
        # Create a contour plot with 10 levels
        fig_multi, ax = plt.subplots()
        plt.gca().set_aspect('equal', adjustable='box')
            
        min_aor_pt_len = 10
        if len(x_up) > min_aor_pt_len and len(y_up) > min_aor_pt_len: 
            plt.scatter(x_up, y_up, color='indigo', label='Max AOR Estimation', s=2)
        if len(x_mean) > min_aor_pt_len and len(y_mean) > min_aor_pt_len:     
            plt.scatter(x_mean, y_mean, color='red', label='Mean AOR Estimation',s=2)
        if len(x_low) > min_aor_pt_len and len(y_low) > min_aor_pt_len: 
            plt.scatter(x_low,y_low,color='tan',label='Min AOR Estimation',s=2)
        
        
        
        plt.gca().set_aspect('equal', adjustable='box')
        if geometry == 0:
            if UnitSys == 2:
                plt.title('AORs with Probability')
                plt.xlabel('X , km')
                plt.ylabel('Y , km')
                
            else:
                plt.title('AORs with Probability')
                plt.xlabel('X , mile')
                plt.ylabel('Y , mile')
        elif geometry == 1:
            plt.title('AORs with Probability')
            plt.xlabel('X , UTM')
            plt.ylabel('Y , UTM')
            
            
            
        plt.axis('equal')
        if geometry == 0:
            # Draw rectangles
            rect1 = plt.Rectangle(((rsvrprop.BXL-rsvrprop.XL)/2/lencoef, (rsvrprop.BYL-rsvrprop.YL)/2/lencoef), rsvrprop.XL/lencoef, rsvrprop.YL/lencoef, linewidth=2, linestyle='--', edgecolor='g', facecolor='none', label='Project Area')
            ax.add_patch(rect1)
             
            rect2 = plt.Rectangle((0, 0), rsvrprop.BXL/lencoef, rsvrprop.BYL/lencoef, linewidth=2, linestyle='--', edgecolor='b', facecolor='none', label='Reservoir Area')
            ax.add_patch(rect2)
            plt.legend(loc='upper right')
            ax.set_xlim([-rsvrprop.BXL/lencoef*0.05, rsvrprop.BXL/lencoef*1.05])
            ax.set_ylim([-rsvrprop.BYL/lencoef*0.05, rsvrprop.BYL/lencoef*1.05])
            #fig_contour.set_facecolor('w')
            ax.set_aspect('equal')
            
            Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
            nWE = simprop.nWE
            # Draw scatter points
            # Create a set to store unique labels
            legend_labels = set()
            for w in range(nWX_ctrl**2):
                x = Xwell[w, nWX_ctrl-1] / lencoef
                y = Ywell[w, nWX_ctrl-1] / lencoef
                ax.scatter(x, y, color='orange', marker='o', s=20)  # Adjust the 's' parameter as desired
                legend_labels.add('Injectors')
            if nWE > 0: 
                for w in range(nWE):
                    x = Xwell[w + nWX_ctrl**2, nWX_ctrl-1] / lencoef
                    y = Ywell[w + nWX_ctrl**2, nWX_ctrl-1] / lencoef
                    ax.scatter(x, y, color='b', marker='^', s=20)  # Adjust the 's' parameter as desired
                    legend_labels.add('Extractors')
            ax.legend()
        
        elif geometry ==1:
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
            nWI = len(WellNumber_inj)
            nWE = len(WellNumber_ext)
            # Get the number of reservoirs
            nReservoirs = AllReservoirs.shape[1] // 2
             # Injectors
            for w in range(nWI):
                ax.scatter(Xwell_input[w], Ywell_input[w], c='orange', s=20)
            
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
           
            plt.axis('equal')
            ax.legend()
        plt.legend(loc='upper right')
        return fig_multi, df_aor_up, df_aor_low, X_center, Y_center
    
    fig_multi, df_aor_up, df_aor_low, X_center, Y_center= plot_multi_AOR()
    
    return fig_tor, fig_multi, SA_values_tor, NPV_values_tor, df_aor_up, df_aor_low, X_center, Y_center



class SAInputStore:
    def __init__(self):
        self.parameter_data = pd.DataFrame(columns=['Parameter', 'Value', 'Min', 'Max'])

    def add_parameter(self, parameter_name, value, min_value, max_value):
        new_row = pd.DataFrame({
          'Parameter': [parameter_name],
          'Value': [value],
          'Min': [min_value],
          'Max': [max_value]
        })
        self.parameter_data = pd.concat([self.parameter_data, new_row], ignore_index=True)


    def get_parameter_data(self):
        return self.parameter_data


def read_SA_input(rsvrprop, relakprop, input_data, UnitSys, uploaded_file):
    #means = [input_data.P0_mean, input_data.temp_mean, input_data.Thickness_mean, input_data.k_mean, input_data.salinity_mean, input_data.m_mean, input_data.n_mean, 
    #         input_data.cr_mean, input_data.FracP_mean, input_data.kra0_mean, input_data.krg0_mean, input_data.Porosity_mean, input_data.Sar_mean, input_data.Sgc_mean]
    input_SA = pd.read_excel(uploaded_file, sheet_name='SensitivityAnalysis', header=None)
    Para_names_excel_list = ['Initial Pressure', 'Temperature', 'Thickness', 'Salinity', 'Porosity', 'Permeability', 'Rock Compressibility', 'Max Allowable Injection Pressure',
                             'Water Corey Exponent (m)', 'Gas Corey Exponent (n)', 'Endpoint Water Relative Permeability (kra0)', 'Endpoint Gas Relative Permeability (krg0)',
                             'Residual Water Saturation (Sar)', 'Critical Gas Saturation (Sgc)']
    attributes = ['P0', 'temp', 'Thickness', 'salinity', 'Porosity', 'k', 'cr', 'FracP', 'm', 'n', 'kra0', 'krg0', 'Sar', 'Sgc']
    coef_SA_SI = [1, 1, 1, 1, 1, 9.87E-16, 1, 1, 1, 1, 1, 1, 1, 1]
    coef_SA_field = [0.00689476, 1, 0.3048, 1.71607E-05, 1, 9.87E-16, 0.000145038, 0.00689476, 1, 1, 1, 1, 1, 1]
    
    SA_input_store = SAInputStore()
    for para_name, attr in zip(Para_names_excel_list, attributes):
        row = input_SA[input_SA.iloc[:, 0] == para_name]
        if not row.empty:
            if UnitSys == 1:
                value = row.iloc[0, 1] * coef_SA_field[attributes.index(attr)]
                min_value = row.iloc[0, 2] * coef_SA_field[attributes.index(attr)]
                max_value = row.iloc[0, 3] * coef_SA_field[attributes.index(attr)]
             # Convert 'Temperature' from F to C if it's the 'Temperature' attribute
                if attr == 'temp':
                    value = (value - 32) * 5/9
                    min_value = (min_value - 32) * 5/9
                    max_value = (max_value - 32) * 5/9
                 
            if UnitSys == 2:
                value = row.iloc[0, 1] * coef_SA_SI[attributes.index(attr)]
                min_value = row.iloc[0, 2] * coef_SA_SI[attributes.index(attr)]
                max_value = row.iloc[0, 3] * coef_SA_SI[attributes.index(attr)]
            
            
            SA_input_store.add_parameter(attr, value, min_value, max_value)
            
    SA_input_data = SA_input_store.get_parameter_data()
    
    
    # Create a list of dictionaries containing parameter names and 'Max' values
    max_values_list = [{'Parameter': attr, 'Max': SA_input_data.loc[SA_input_data['Parameter'] == attr, 'Max'].values[0]} for attr in attributes]
    max_values_df = pd.DataFrame(max_values_list)
    
    # Create a list of dictionaries containing parameter names and 'Min' values
    min_values_list = [{'Parameter': attr, 'Min': SA_input_data.loc[SA_input_data['Parameter'] == attr, 'Min'].values[0]} for attr in attributes]
    min_values_df = pd.DataFrame(min_values_list)

    
    return SA_input_data, max_values_df, min_values_df


### TESTING for density histogram plot
# Read input properties from the uploaded Excel file, for testing purposes

'''
from properties import Unit_conversion, read_rsvr_input, userexcel_conversion
import pandas as pd
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_SI.xlsx'
UnitSys = 2
BCValue = 1
geometry = 0
userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
userinput_inj = pd.read_excel(uploaded_file, sheet_name='Injectors', header=None)
userinput_ext = pd.read_excel(uploaded_file, sheet_name='Extractors', header=None)
input_df, proj_name = userexcel_conversion(userinput_df, userinput_inj, userinput_ext)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)
SA_input_data, max_values_df, min_values_df = read_SA_input(rsvrprop, relakprop, output_print, UnitSys, uploaded_file)
## Test 1

Sensi_para_rel_num = 100
hist_data = generate_histdata(Sensi_para_rel_num, rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, geometry, output_print)
fig_hist = plot_density_curve(hist_data, geometry, UnitSys)

## TEST2 for tornado chart
#real_num = 100
#Sensi = 1

#means = [rsvrprop.P0, rsvrprop.temp, rsvrprop.Thickness, rsvrprop.k, rsvrprop.salinity, relakprop.m, relakprop.n, rsvrprop.cr, rsvrprop.FracP,
#        relakprop.kra0, relakprop.krg0, rsvrprop.Porosity, relakprop.Sar, relakprop.Sgc]
#sds = [mean * 0.01 for mean in means]
#input_data = Inputclass(real_num, means, sds)
#fig_tor, fig_multi, SA_values_tor, NPV_values_tor, df_aor_up, df_aor_low, X_center, Y_center = sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop, BCValue, output_print, UnitSys, geometry, uploaded_file,proj_name)
'''