from shapely.geometry import Polygon
from math import atan2
from plotting import S1S2_PContour
from S2 import Scenario2
import numpy as np
from S3 import S3_core
import streamlit as st

def calculate_small_polygon_area(small_loc, large_loc):
    x_loc_sm = [x for x,y in small_loc]
    y_loc_sm = [y for x,y in small_loc]
    x_loc_lg = [x for x,y in large_loc]
    y_loc_lg = [y for x,y in large_loc]
    
    polygon_lg = Polygon(zip(x_loc_lg, y_loc_lg))
    polygon_sm = Polygon(zip(x_loc_sm,  y_loc_sm))
    
    intersection = polygon_lg.intersection(polygon_sm)
    
    return intersection.area

def calculate_S2_storeff(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, ave_sat):
    Sensi = 0
    ConRate = 0
    
    if UnitSys == 2:
        lencoef = 1000
    elif UnitSys == 1:
        lencoef = 1609.34
    else: 
        raise Exception ("Wrong Unit System selection!")
        
    Psave2, _, rL1, rT1, _, _, Qsave_res = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
    
    highlight_locations, _,_, highlight_contour = S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2, Qsave_res, ConRate)

    large_loc = [(0,0), (rsvrprop.BXL/lencoef,0),(rsvrprop.BXL/lencoef,rsvrprop.BYL/lencoef), (0, rsvrprop.BYL/lencoef)]
    total_plum_area = np.sum((rL1[:, nWX_ctrl-1]/lencoef) ** 2 * np.pi)
    
    if 'highlight_locations' in locals() and len(highlight_locations)>=4: 
        aor_area = calculate_small_polygon_area(highlight_locations, large_loc)
    else:
        aor_area = rsvrprop.BXL/lencoef*rsvrprop.BYL/lencoef
        st.warning('Cannot visualize the AOR inside the reservoir boundaries. Total AOR is taken as the reservoir area.')
    stor_eff = total_plum_area*ave_sat/aor_area
    return stor_eff


def calculate_S3_storeff(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, ave_sat, uploaded_file):
    # Make Sensi = 1 here 
    Sensi = 1
    ConRate = 0
    
    Psave2_ave, _ , highlight_locations, X_center, Y_center , rL1, rT1 = S3_core(rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, Sensi)
    
    x_start = X_center-rsvrprop.BXL/2
    y_start = Y_center-rsvrprop.BYL/2
    x_end = X_center+rsvrprop.BXL/2
    y_end = Y_center+rsvrprop.BYL/2
    large_loc = [(x_start,y_start), (x_end,y_start),(x_end,y_end), (x_start, y_end)]
    
    total_plum_area = np.sum((rL1[:]) ** 2 * np.pi)
    if 'highlight_locations' in locals() and len(highlight_locations)>=4: 
        aor_area = calculate_small_polygon_area(highlight_locations, large_loc)
    else:
        aor_area = rsvrprop.BXL*rsvrprop.BYL
        st.warning('Cannot visualize the AOR inside the reservoir boundaries. Total AOR is taken as the reservoir area.')
        
    stor_eff = total_plum_area*ave_sat/aor_area
    return stor_eff

'''
# Example coordinates of the larger polygon (domain limits)
larger_polygon_x = [0, 2, 2, 1, 0]
larger_polygon_y = [0, 1, 2, 2, 1]

# Example coordinates of the small polygon
small_polygon_x = [1, 2, 2, 1]
small_polygon_y = [0, 1, 2, 2]


small_polygon_x = [1, 1, 2, 2]
small_polygon_y = [0, 2, 2, 1]

small_vertices = [(x, y) for x, y in zip(small_polygon_x, small_polygon_y)]

polygon = Polygon(small_vertices)

is_ccw = polygon.exterior.is_ccw

if not is_ccw:
    centroid = polygon.centroid
    sorted_vertices = sorted(small_vertices, key=lambda vertex: atan2(vertex[1] - centroid.y, vertex[0] - centroid.x))
    reordered_polygon = Polygon(sorted_vertices)
else:
    reordered_polygon = polygon

'''



from properties import Unit_conversion, read_rsvr_input, userexcel_conversion
import pandas as pd
import numpy as np
sheet_name = 'Parameters'
uploaded_file = 'ET_input_temp_field.xlsx'
UnitSys = 1

# currently only closed boundary, BCValue = 1
BCValue = 1

userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
userinput_inj = pd.read_excel(uploaded_file, sheet_name='Injectors', header=None)
userinput_ext = pd.read_excel(uploaded_file, sheet_name='Extractors', header=None)
input_df, proj_name = userexcel_conversion(userinput_df, userinput_inj, userinput_ext)

output_df, output_print = Unit_conversion(input_df, UnitSys)
rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)


nWX_ctrl = simprop.nWX_npv_max

ave_sat = 0.5


stor_eff_S2 = calculate_S2_storeff(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, ave_sat)
stor_eff_S3 = calculate_S3_storeff(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, ave_sat, uploaded_file)