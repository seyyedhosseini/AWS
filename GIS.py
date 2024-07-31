
import folium
from folium.plugins import MiniMap, MeasureControl, BeautifyIcon
import pandas as pd
import pyproj
import streamlit as st

import numpy as np
import utm
#from folium import raster_layers
#from io import BytesIO

def GIS_map(filename, df_aor_loc, df_aor_up, df_aor_low, zone_number, UnitSys, X_center, Y_center, rsvrprop):

    
    
    if UnitSys == 2:
        lencoef = 1000
    elif UnitSys == 1:
        lencoef = 1609.34
    else: 
        raise Exception ("Wrong Unit System selection!")
    
    def wells_rsvr_GIS():
        
        # 1. Convert UTM into coordinates
        
        
        df_well_int = pd.read_excel(filename,sheet_name='Injectors')
        df_well_ext = pd.read_excel(filename,sheet_name='Extractors')
        df_well = pd.concat([df_well_int, df_well_ext], ignore_index=True)
        df_geo_orig = pd.read_excel(filename,sheet_name='Geometry')
        df_geo = df_geo_orig.dropna(axis=1, how='all')
        
        #map = folium.Map(location=[ref_lat, ref_lon], zoom_start=12)
        
        # 1.Mark wells locations

        well_utm_x = df_well_int['UTM_Well_X']
        well_utm_y = df_well_int['UTM_Well_Y']
    
        # injectors plotting
        well_lat, well_lon = convert_to_latlon(well_utm_x, well_utm_y, zone_number)
        
        map = folium.Map(location=[well_lat[0], well_lon[0]], zoom_start=12)
        
        for lat,lon,well_num in zip(well_lat, well_lon, df_well_int["Well Number"]):
            folium.Marker(location=[lat,lon],popup=f"Injector {well_num}",
            icon=folium.Icon(color="green")).add_to(map)
        # extractors plotting    

        well_utm_x = df_well_ext['UTM_Well_X']
        well_utm_y = df_well_ext['UTM_Well_Y']    
        well_lat, well_lon = convert_to_latlon(well_utm_x, well_utm_y, zone_number)
        for lat,lon,well_num in zip(well_lat, well_lon, df_well_ext["Well Number"]):
            folium.Marker(location=[lat,lon],popup=f"Extractor {well_num}",
            icon=folium.Icon(color="blue")).add_to(map)

        
        df_geo_utm = df_geo.copy()
        lat_geo = []
        lon_geo = []

        
        for i in range(1, (df_geo.shape)[1] // 2 + 1):  
            column_name_x = f'ProjX{i}_UTM'
            column_name_y = f'ProjY{i}_UTM'
            
            #df_geo_utm is a 10 columns df
            
            df_geo_utm_wonan = df_geo_utm.dropna(subset=[column_name_x, column_name_y]) # drop nan
            lat, lon = convert_to_latlon(df_geo_utm_wonan[column_name_x], df_geo_utm_wonan[column_name_y], zone_number)
            lat_geo.append(lat)
            lon_geo.append(lon)
        
        
        
        
        # 2. Reservoir boundaries 
        for lat_line, lon_line in zip(lat_geo, lon_geo):
            # Create a polyline object for the boundary line
            polyline = folium.PolyLine(
                locations=list(zip(lat_line, lon_line)),popup="Reservoir boundary",
                color='black',
                weight=2
            )
            
            # Add the polyline to the map
            polyline.add_to(map)
        return map
            
    map = wells_rsvr_GIS()
            
    
    def AOR_GIS():
        # 3. AOR boundaries, based on UTM coordiantes
        
        lat_aor = []
        lon_aor = []
        #df_aor_loc = pd.read_excel('highlight_locations.xlsx')
        lat, lon = convert_to_latlon(df_aor_loc['x'], df_aor_loc['y'], zone_number)
        lat_aor.append(lat)
        lon_aor.append(lon)
        # AOR boundaries 
        for lat_line, lon_line in zip(lat_aor, lon_aor):
            ## Create a polyline object for the boundary line
            
            #polyline = folium.PolyLine(
            #    locations=list(zip(lat_line, lon_line)),popup = "Area of Review",
            #    color='red',
            #    weight=2
            #)
            ## Add the polyline to the map
            #polyline.add_child(folium.Popup('AOR')).add_to(map)
            
            for lat, lon in zip(lat_line, lon_line):
                # Create a scatter marker for each point
                marker = folium.CircleMarker(
                    location=[lat, lon],
                    radius = 0.5,
                    color = 'red',
                    fill =True,
                    fill_color = 'red',
                    popup='AOR Point',
                    #icon=folium.Icon(color='red', icon='circle', icon_size=(0.5,0.5))
                )
                # Add the marker to the map
                #circle_marker.add_to(map)
                marker.add_child(folium.Popup('Mean AOR Estimation')).add_to(map)
        
        lat_aor = []
        lon_aor = []
        lat, lon = convert_to_latlon(df_aor_up['x'], df_aor_up['y'], zone_number)
        lat_aor.append(lat)
        lon_aor.append(lon)
        # AOR boundaries 
        for lat_line, lon_line in zip(lat_aor, lon_aor): 
            
            for lat, lon in zip(lat_line, lon_line):
                # Create a scatter marker for each point
                marker = folium.CircleMarker(
                    location=[lat, lon],
                    radius = 0.5,
                    color = 'indigo',
                    fill =True,
                    fill_color = 'indigo',
                    popup='AOR Point',
                    #icon=folium.Icon(color='red', icon='circle', icon_size=(0.5,0.5))
                )
                # Add the marker to the map
                #circle_marker.add_to(map)
                marker.add_child(folium.Popup('Max AOR Estimation')).add_to(map)
        lat_aor = []
        lon_aor = []        
        lat, lon = convert_to_latlon(df_aor_low['x'], df_aor_low['y'], zone_number)
        lat_aor.append(lat)
        lon_aor.append(lon)
        # AOR boundaries 
        for lat_line, lon_line in zip(lat_aor, lon_aor):
           
            
            for lat, lon in zip(lat_line, lon_line):
                # Create a scatter marker for each point
                marker = folium.CircleMarker(
                    location=[lat, lon],
                    radius = 0.5,
                    color = 'gray',
                    fill =True,
                    fill_color = 'gray',
                    popup='AOR Point',
                    #icon=folium.Icon(color='red', icon='circle', icon_size=(0.5,0.5))
                )
                # Add the marker to the map
                #circle_marker.add_to(map)
                marker.add_child(folium.Popup('Min AOR Estimation')).add_to(map)
                
        # Rectangle properties
        
        rectangle_center = [X_center, Y_center]
        rectangle_width = rsvrprop.BXL
        rectangle_height = rsvrprop.BYL
        
        # Define the rectangle's vertices
        vertices = [
            (rectangle_center[0] - rectangle_width/2, rectangle_center[1] - rectangle_height/2),
            (rectangle_center[0] + rectangle_width/2, rectangle_center[1] - rectangle_height/2),
            (rectangle_center[0] + rectangle_width/2, rectangle_center[1] + rectangle_height/2),
            (rectangle_center[0] - rectangle_width/2, rectangle_center[1] + rectangle_height/2)
        ]
        vertices_df = pd.DataFrame(vertices, columns=['x', 'y'])
       
        lat, lon = convert_to_latlon(vertices_df['x'],vertices_df['y'], zone_number)
        coordinates_list = [(lat_val, lon_val) for lat_val, lon_val in zip(lat, lon)]
        # Create a Folium polygon representing the rectangle
        rectangle = folium.Polygon(
            locations=coordinates_list,
            fill=False,
            color='blue',
            line_opacity=1,
            line_weight=2,
            line_dasharray='10, 10'  # Dash pattern
        )
        
        # Add the rectangle to the map
        rectangle.add_child(folium.Popup('Reservoir Boundary')).add_to(map)
        
        return map
    map = AOR_GIS()
        
        
    '''
    # 4. Add Basin Bounds
    
    # Define the UTM coordinates for the rectangular bound
    utm_coords = [convert_to_latlon(utm_x, utm_y, zone_number, zone_letter), convert_to_latlon(utm_x + 20000, utm_y + 10000, zone_number, zone_letter)]
    
    
    # Create a rectangle overlay using the UTM coordinates
    rectangle = folium.Rectangle(
        bounds=utm_coords,
        color='red',
        fill=False,
        weight=2
    )
    
    # Add the rectangle overlay to the map
    rectangle.add_to(map)
    '''
    #map.save('map.html')
    # Add the MeasureControl plugin for the scale bar
    measure_control = MeasureControl(position='bottomleft', primary_length_unit='meters')
    map.add_child(measure_control)
    # Add a scale bar to the map
    minimap = MiniMap(toggle_display=True, position="bottomright", metric=True)
    map.add_child(minimap)
    # Add search capability
    draw_options = folium.plugins.Draw(export=False, filename='data.geojson', position='topleft', show_geometry_on_click=True, draw_options=None, edit_options=None)
    map.add_child(draw_options)
    full_screen = folium.plugins.Fullscreen(position='topright', title='Full Screen', title_cancel='Exit Full Screen', force_separate_button=False)
    map.add_child(full_screen)
    mouse_position = folium.plugins.MousePosition(position='bottomleft', separator=' : ', empty_string='Unavailable', lng_first=False,
                           num_digits=5, prefix='', lat_formatter=None, lng_formatter=None)
    map.add_child(mouse_position)

    # # Define baselayer options
    # baselayer_options = {
    #     "Open Street Map": folium.TileLayer("openstreetmap"),
    #     "Stamen Terrain": folium.TileLayer("stamenterrain"),
    #     "Stamen Toner": folium.TileLayer("stamentoner"),
    #     # Add more baselayers as needed
    # }
    # map.add_child(baselayer_options["Open Street Map"])
    # # Add a dropdown for baselayer selection
    # selected_baselayer = st.selectbox("Select Baselayer", list(baselayer_options.keys()))
    # # Update the map with the selected baselayer
    # map.add_child(baselayer_options[selected_baselayer])

    return map


def convert_to_utm(lat, lon, zone_number):
    # Define the UTM projection
    utm_proj = pyproj.Proj(proj="utm", zone=zone_number, ellps="WGS84", datum="WGS84", units="m")

    # Convert the geographical coordinates to UTM
    utm_x, utm_y = utm_proj(lon, lat)

    return utm_x, utm_y


def convert_to_latlon(utm_x, utm_y, zone_number):
    # Define the UTM projection
    utm_proj = pyproj.Proj(proj="utm", zone=zone_number, 
                           ellps="WGS84", datum="WGS84", units="m")

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(utm_x, utm_y, inverse=True)

    return lat, lon

def get_utm_zone(lat, lon):
    utm_zone = utm.from_latlon(lat, lon)
    zone_number = utm_zone[2]
    zone_letter = utm_zone[3]
    return zone_number, zone_letter


'''
## TESTING
from properties import Unit_conversion, read_rsvr_input
from S3 import Scenario3
sheet_name = 'Parameters'
uploaded_file = 'Inputs_SI _UTM.xlsx'
UnitSys = 2
BCValue = 1
ref_lat = 29.7299
ref_lon = -94.2913


geometry = 1
_, _, _, _, df_aor_up, df_aor_low, X_center, Y_center = sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop,
                                              BCValue,    output_print, UnitSys, geometry, uploaded_file)
map = GIS_map(uploaded_file, df_aor, df_aor_up, df_aor_low, zone_number, UnitSys, X_center, Y_center, rsvrprop)
# Render the map as an HTML string
map_html = map._repr_html_()
map.width = '100%'
map.height = '600px'
st.components.v1.html(map_html, width=1350, height=1000, scrolling=True)
'''