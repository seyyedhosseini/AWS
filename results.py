# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:44:46 2023

@author: wangz
"""

from Matcal import Distance_Mat_Builder
import numpy as np

#def results_output(output_print,BCValue, ConRate, NPV, Capacity, simprop,rsvrprop, UnitSys, Qsave2, Psave2, X3):
def results_output(output_print,BCValue, ConRate, NPV, Capacity, simprop,rsvrprop, UnitSys, Qsave2, Psave2, proj_name, 
                   Pnode, rL1):   
    nWXmin = simprop.nWXmin
    nWXmax = simprop.nWXmax
    nWE = simprop.nWE
    Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
    Xnode = 50
    Ynode = 50
    with open('EASiToolOutput.txt', 'w') as file:
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
        
        if ConRate == 1:
            file.write('%43s %25s\n\n' % ('Injection/Extraction Type', 'Uniform'))
        elif ConRate == 0:
            file.write('%43s %25s\n\n' % ('Injection/Extraction Type', 'Optimal'))

        if UnitSys == 1:
            file.write('%43s %25s\n\n' % ('Unit System', 'Field Units'))
        elif UnitSys ==2:
            file.write('%43s %25s\n\n' % ('Unit System', 'SI Units'))       
            
        file.write('\n')
        file.write('*********************************************************************\n') 
        file.write('%43s\n' % 'EASiTool Output Data')
        file.write('*********************************************************************\n') 
        file.write('{:>70}\r'.format('     Number of        Capacity              NPV                       '))
        file.write('{:>70}\r'.format('     Injectors    (M metric ton)           ($M)                       '))
        file.write('*********************************************************************\n')  
        
        for nWX in range(nWXmin, nWXmax + 1):
            file.write('%10u' % (nWX * nWX))
            file.write('%20.3e' % Capacity[nWX-1])
            file.write('%20.3e\n' % NPV[nWX-1])
        file.write('*********************************************************************\n')
            
        for nWX in range(nWXmin, nWXmax + 1):
            file.write('*******************************************************************************\n')
            file.write('%41s %5u\n' % ('Number of Injectors =', (nWX * nWX)))
            file.write('%41s %5u\n' % ('Number of Extractors =', nWE))
            file.write('*******************************************************************************\n')
            if UnitSys == 2: # SI units 
                file.write('%80s\n' % ' Well   Well     X        Y      Inj. Rate      Ext. Rate   Pressure    Plume radius')
                file.write('%80s\n' % ' #      Type    (m)      (m)     (MMT/yr)      (m^3/day)      (MPa)        (m)')
                file.write('%s\n' % '*****************************************************************************')
                for nW in range(0, (nWX * nWX)):
                    file.write('%4u' % (nW+1))
                    file.write('%10s' % 'Inj.')
                    file.write('%11.3e' % Xwell[nW, nWX-1])
                    file.write('%11.3e' % Ywell[nW, nWX-1])
                    #file.write('%11.3e' % (Qsave2[nW, nWX-1])) #ton/day
                    file.write('%11.3e' % (Qsave2[nW, nWX-1]*0.000365)) #MMT/yr
                    file.write('%11s' % 'N/A')
                    file.write('%11.3e' % (Psave2[nW, nWX-1])) #MPa
                    file.write('%11.3e\n' % (rL1[nW, nWX-1])) # m

                
                for nW in range((nWX * nWX), (nWX * nWX) + nWE):
                    file.write('%4u' % (nW+1))
                    file.write('%10s' % 'Ext.')
                    file.write('%11.3e' % Xwell[nW, nWX-1])
                    file.write('%11.3e' % Ywell[nW, nWX-1])
                    file.write('%11s' % 'N/A')
                    file.write('%11.3e' % Qsave2[nW, nWX-1])  # m^3/day
                    file.write('%11.3e' % Psave2[nW, nWX-1])
                    file.write('%11s\n' % 'N/A')
                    
                file.write('\n')
                file.write('*********************************************************************\n') 
                    
                file.write('%43s\n' % 'Basin Pressure (MPa) Data (50 by 50)')
                file.write('*********************************************************************\n') 
                # Reshape the list into a 50x50 matrix
                P_matrix = np.reshape(Pnode, (Xnode, Ynode))
                for row in P_matrix:
                    rounded_row = [round(value, 2) for value in row]
                    file.write(" ".join(map(str, rounded_row)) + "\n")
                    
                    
            elif UnitSys == 1: # Field Units
                file.write('%70s\n' % ' Well   Well     X        Y      Inj. Rate      Ext. Rate   Pressure    Plume radius ')
                file.write('%70s\n' % ' #      Type   (mile)   (mile)    (MMT/yr)     (bbl/day)      (psi)         (mile)')
                file.write('%s\n' % '*******************************************************************************')
                for nW in range(0, (nWX * nWX)):
                    file.write('%4u' % (nW+1))
                    file.write('%10s' % 'Inj.')
                    file.write('%11.3e' % (Xwell[nW, nWX-1]/1609.34)) # mile
                    file.write('%11.3e' % (Ywell[nW, nWX-1]/1609.34)) # mile
                    #file.write('%11.3e' % (Qsave2[nW, nWX-1])) #ton/day
                    file.write('%11.3e' % (Qsave2[nW, nWX-1]*0.000365)) # MMT/year
                    file.write('%11s' % 'N/A')
                    file.write('%11.3e' % (Psave2[nW, nWX-1] *145.038)) #psi
                    file.write('%11.3e\n' % (rL1[nW, nWX-1]/1609.34)) # mile
  
                
                for nW in range((nWX * nWX), (nWX * nWX) + nWE):
                    file.write('%4u' % (nW+1))
                    file.write('%10s' % 'Ext.')
                    file.write('%11.3e' % (Xwell[nW, nWX-1]/1609.34)) # mile
                    file.write('%11.3e' % (Ywell[nW, nWX-1]/1609.34)) # mile
                    file.write('%11s' % 'N/A')
                    file.write('%11.3e' % (Qsave2[nW, nWX-1] * 6.2898107))  # bbl/day
                    file.write('%11.3e' % (Psave2[nW, nWX-1] *145.038)) #psi
                    file.write('%11s\n' % 'N/A')
                
                file.write('\n')
                file.write('*********************************************************************\n') 
                file.write('%43s\n' % 'Basin Pressure (psi) Data (50 by 50)')
                file.write('*********************************************************************\n') 
                # Reshape the list into a 50x50 matrix
                P_matrix = np.reshape((np.array(Pnode) * 145.03768), (Xnode, Ynode))
                for row in P_matrix:
                    rounded_row = [round(value, 2) for value in row]
                    file.write(" ".join(map(str, rounded_row)) + "\n")
            else: 
                raise Exception('Wrong Unit System Selection!')
                        

            file.write('%s\n\n\n\n\n' % ' ')
        
        # if need to ouptut pressure contour     
        #for row in X3:
        #    for element in row:
        #        file.write(str(element) + ' ')
        #    file.write('\n')
                            
    file.close()



