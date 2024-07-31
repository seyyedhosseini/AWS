# --Import libraries--
import streamlit as st
import streamlit_authenticator as stauth
import folium
import pandas as pd
import sys
import logging
import datetime
from S1 import Scenario1
from S2 import Scenario2
from S3 import Scenario3
from properties import Unit_conversion, read_rsvr_input, check_parameter_ranges, userexcel_conversion
from plotting import plot_plume_extension, plot_well_property, plot_AOR, plot_npv_capacity, S1S2_PContour, warning_pressure_plume
from Matcal import Distance_Mat_Builder
from Sensitivity import generate_histdata, plot_density_curve, sensi_ana_tor
import base64
from GIS import GIS_map
from S3 import Scenario3
import warnings
from results import results_output
import re
import numpy as np
import yaml
from yaml.loader import SafeLoader
import plotly.graph_objects as go
# Ignore runtime warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Define the file path

file_path_1 = "ET_input_temp_SI.xlsx"
file_path_2 = 'ET_input_temp_field.xlsx'


# Define a function to download the Excel file
def download_excel(file_path):
    with open(file_path, "rb") as file:
        b64 = base64.b64encode(file.read()).decode()  # Convert file to bytes
        href = f'<a href="data:application/octet-stream;base64,{b64}" download="{file_path}"><font color="#0000FF"><font size="+5">Download Input File Template</font></a>'
        return href


def download_output(file_path_output):
    with open(file_path_output, "rb") as file:
        b64 = base64.b64encode(file.read()).decode()  # Convert file to bytes
        href = f'<a href="data:application/octet-stream;base64,{b64}" download="{file_path_output}"><font color="#0000FF"><font size="+5">Download Simulation Output File</font></a>'
        return href

def download_SA_output(file_path_output):
    with open(file_path_output, "rb") as file:
        b64 = base64.b64encode(file.read()).decode()  # Convert file to bytes
        href = f'<a href="data:application/octet-stream;base64,{b64}" download="{file_path_output}"><font color="#0000FF"><font size="+5">Download Sensitivity Analysis Output File</font></a>'
        return href

def download_uplodad_file():
    # Display the download button
    st.markdown(download_excel(), unsafe_allow_html=True)
    st.markdown(""":red[Please download the Excel input file and fill in the ] **:blue[input parameters.]**""")

    # --Add file uploader to select an input Excel file--
    uploaded_file = st.file_uploader('Select Input Excel file', type=['xlsx'])
    return uploaded_file


def reset_app():
    # Reset session state variables or perform any other necessary actions
    st.session_state["Run EASiTool"] = False
    st.session_state["Run Sensitivity Analysis"] = False
    #st.session_state["Generate GIS"] = False
    # Clear all display and cache
    # st.experimental_rerun()


def show_license_window():
    license_text = "Copyright <2023> <Gulf Coast Carbon Center>\n\n"\
                    "Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n\n"\
                    "1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n"\
                    "2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n"\
                    "3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n\n"\
                    "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, " \
                    "INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. " \
                    "IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, " \
                    "PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED " \
                    "AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT " \
                    "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." \

    st.write('### License Agreement')
    accepted_terms = st.text_area("License Agreement", value=license_text, height=350,
                                  key="license_agreement",label_visibility='hidden')
    accept_terms = st.checkbox("I accept the terms and conditions.")
    continue_button = st.button("Continue")



    # Check if the continue button is clicked and the terms are accepted
    if continue_button and accept_terms:
        st.session_state.show_license = False
        st.experimental_rerun()  # Trigger a rerun to execute the main program


# --Define the main function--


global df

def main():
    # Sidebar
    st.set_page_config(layout="wide")
    st.image('E.jpg', width=200)
    st.title('EASiTool 5.0 Alpha Version')

    # Set the CSS styles for center alignment and larger font
    centered_style = "<style> .centered { text-align: center; }</style>"
    larger_font_style = "<style> .larger-font { font-size: 24px; }</style>"

    # Add content to the sidebar
    st.sidebar.title("Welcome to EASiTool 5.0")
    st.sidebar.write(
        "This tool incorporates analytical-based solutions to provide CO₂ storage capacity estimation for Geological Carbon Storage (GCS) and offers a range of powerful features to support efficient CO₂ storage evaluation and decision-making.")
    st.sidebar.write(
        "This app is developed by the [Gulf Coast Carbon Center](https://gccc.beg.utexas.edu/) at the [Bureau of Economic Geology](https://www.beg.utexas.edu/).")
    st.sidebar.title("How to Use EASiTool 5.0?")

    # Step 1
    st.sidebar.title("Step 1")
    st.sidebar.write("Download the template input file and edit as needed. "
                     " (Data in the template is just an example and does not represent an actual project.) "
                     )
    st.sidebar.write("Select the module from the input file and type in the required parameters.")

    # Module Descriptions
    module_text1 = "**User Given Inputs (Given Geometry, Constant Injection/Extraction Rate)**"
    module_text1 += " - Choose this module if you have project and well locations, and injection/extraction rates."
    st.sidebar.markdown(module_text1, unsafe_allow_html=True)

    module_text2 = "**Maximum Storage Capacity (Fixed Bottomhole Pressure)**"
    module_text2 += " - Start with this module for a general project without specific location and well information."
    st.sidebar.markdown(module_text2, unsafe_allow_html=True)

    # Step 2
    st.sidebar.title('Step 2')
    st.sidebar.write("Upload the input file, then select the reservoir boundary condition.")

    st.sidebar.title("Step 3")
    st.sidebar.write("Double check your inputs. Run the tool, and EASiTool 5.0 generates the following results "
                     "for your GCS project:")
    st.sidebar.write("1. Pressure contour map\n"
                     "2. CO₂ plume extension map\n"
                     "3. Well pressure/flow rates map*\n"
                     "4. Area of Review (AoR) evaluations\n"
                     "5. Sensitivity Analysis\n"
                     "6. Geographic Information System (GIS) map*\n"
                     "7. NPV and Capacity vs. number of injection optimization*\n")
    st.sidebar.write("*: For the certain module")

    # Sensitivity Analysis
    st.sidebar.title("Step 4")
    st.sidebar.write("If you need to run Sensitivity Analysis, EASiTool displays:")
    st.sidebar.write("1. Tornado chart\n"
                     "2. Sensitivity analysis of capacity estimation using the Monte Carlo Simulation ")

    # Additional Information
    st.sidebar.write("(MMT means million metric tons)")
    st.sidebar.write(
        "(If using the user provided inputs, in the geometry tab, please make sure that the first input point and the last input point of each project area are the same point.)")

    # Feedback
    st.sidebar.write("**[Send us your feedback.](mailto:your-email@example.com)**")

    #add color to bottons

    m = st.markdown("""
    <style>
    div.stButton > button:first-child {
        background-color: #0099ff;
        color:#ffffff;
    }
    div.stButton > button:hover {
        background-color: #00ff00;
        color:#ff0000;
        }
    </style>""", unsafe_allow_html=True)

    css = '''
    <style>
        .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
        font-size:1.5rem;
        }
    </style>
    '''

    st.markdown(css, unsafe_allow_html=True)

    #Add feedback area to sidebar
    with st.sidebar:
        contact_form = """
        <form action="https://formsubmit.co/easitool@beg.utexas.edu" method="POST">
             <input type="hidden" name="_captcha" value="false">
             <input type="text" name="name" placeholder="Your name" required>
             <input type="email" name="email" placeholder="Your email" required>
             <textarea name="message" placeholder="Your feedback"></textarea>
             <button style="color:white; background-color:MediumSeaGreen;border-color:Gainsboro" type="submit">Send</button>
        </form>    
        """
        st.markdown(contact_form, unsafe_allow_html=True)

        # Use Local CSS File
        def local_css(file_name):
            with open(file_name) as f:
                st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

        local_css("style/style.css")

    with st.sidebar:
        url = 'https://gccc.beg.utexas.edu/easitool'

        st.markdown(f'''
        <a href={url}><button style="color:white; background-color:MediumSeaGreen;border-color:Gainsboro;padding: 10px;" type="submit">Visit our website</button></a>
        ''',
                    unsafe_allow_html=True)

    # Initialize session state variable
    if "show_license" not in st.session_state:
        st.session_state.show_license = True

    if st.session_state.show_license:
        show_license_window()
    else:
        # Start the app
        ##username/password
        ##these two lines are to generate new passwords for config.YAML file
        #hashed_passwords = stauth.Hasher(['password', 'training']).generate()
        #print(hashed_passwords)

        with open('config.yaml') as file:
            config = yaml.load(file, Loader=SafeLoader)
            authenticator = stauth.Authenticate(
                config['credentials'],
                config['cookie']['name'],
                config['cookie']['key'],
                config['cookie']['expiry_days'],
                config['preauthorized']
            )
        name, authentication_status, username = authenticator.login('Login', 'main')

        # # track user login
        # logger = logging.getLogger('ET')
        # logger.setLevel(logging.INFO)
        # logger.propagate = False
        # formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
        #                               datefmt='%Y-%m-%d %H:%M:%S')
        # fh = logging.FileHandler('ET.log')
        # fh.setFormatter(formatter)
        # if (logger.hasHandlers()):
        #     logger.handlers.clear()
        # logger.addHandler(fh)
        # if authentication_status:
        #     authenticator.logout('Logout', 'main')
        #     st.write(f"Logged in as: {name} ")
            # logger.info(f"User '{name}' logged in.")
        ## track user login
        #logger = logging.getLogger('ET')
        #logger.setLevel(logging.INFO)
        #logger.propagate = False
        #formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
        #                              datefmt='%Y-%m-%d %H:%M:%S')
        #fh = logging.FileHandler('ET.log')
        #fh.setFormatter(formatter)
        #if (logger.hasHandlers()):
            #logger.handlers.clear()
        #logger.addHandler(fh)

        if authentication_status:
            authenticator.logout('Logout', 'main')
            st.write(f"Logged in as: {name} ")
            #logger.info(f"User '{name}' logged in.")
        elif authentication_status == False:
            st.error('Username/password is incorrect')
            st.stop()
        elif authentication_status == None:
            st.warning('Please enter your username and password')
            st.stop()

        #Select Units

        Unit_col1, unit_col2 = st.columns(2)
        with Unit_col1:
            st.markdown('<h3 style="margin-bottom: 0;">Select unit system:</h3>', unsafe_allow_html=True)
            function3 = st.selectbox(
                ' ',
                ('SI Unit', 'Field Unit')
            )

        with unit_col2:
            # st.write('Please select your input files based on your unit system:')
            # Display the download button
            st.write('')
            st.write('')
            st.write('')
            st.write('')
            st.write('')
            if function3 == 'SI Unit':
                st.markdown(download_excel(file_path_1), unsafe_allow_html=True)
            elif function3 == 'Field Unit':
                st.markdown(download_excel(file_path_2), unsafe_allow_html=True)
            else:
                st.write('Please select the correct unit system.')


        # Display the CSS styles
        st.markdown(centered_style, unsafe_allow_html=True)
        st.markdown(larger_font_style, unsafe_allow_html=True)
        #st.write(
        #    """:red[Please download the Excel input file and edit the input values. Data in the template is just an example and is not representing an actual project. ]""")


        # --Add file uploader to select an input Excel file--
        uploader_col1, uploader_col2 = st.columns(2)
        with uploader_col1:
            st.markdown('<h3 style="margin-bottom: 0;">Upload the input file here:</h3>', unsafe_allow_html=True)
            uploaded_file = st.file_uploader(' ', type=['xlsx'])
        with uploader_col2:
            if uploaded_file is not None:
                st.markdown('<h3 style="margin-bottom: 0;">Select reservoir boundary condition:</h3>',
                            unsafe_allow_html=True)
                function2 = st.selectbox(
                    ' ',
                    ('Closed Boundary', 'Open Boundary')
                )

        if uploaded_file is not None:
            # Read input properties from the uploaded Excel file
            Sensi = 0
            file_name = uploaded_file.name
            sheet_name = 'Parameters'
            try:
                #input_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
                userinput_df = pd.read_excel(uploaded_file, sheet_name=sheet_name, header=None)
                userinput_inj = pd.read_excel(uploaded_file, sheet_name='Injectors', header=None)
                userinput_ext = pd.read_excel(uploaded_file, sheet_name='Extractors', header=None)
                input_df, proj_name = userexcel_conversion(userinput_df, userinput_inj, userinput_ext)
            except:
                print('this is the wrong file:', file_name)
                st.write('This is the wrong file. Please upload the correct input Excel file and fill in the parameters.')
                sys.exit()



            if input_df.iloc[1, 28] == 1:
                function1 = 'User Given Inputs'
            elif input_df.iloc[1, 28] == 2:
                function1 = 'Maximum Storage Capacity'
            else:
                # You can add a default value here in case the condition doesn't match
                function1 = 'Maximum Storage Capacity'  # Change 'Default Value' to your desired default
                st.write(' There is an error in your module selection. Maximum Storage Capacity Module is selected here.')


            if function3 == 'Field Unit':
                UnitSys = 1
            elif function3 == 'SI Unit':
                UnitSys = 2
            else:
                # st.write('Invalid boundary condition selected')
                UnitSys = 3
            if function2 == 'Closed Boundary':
                BCValue = 1
            elif function2 == 'Open Boundary':
                BCValue = 2
            else:
                # st.write('Invalid boundary condition selected')
                BCValue = 3

            output_df, output_print = Unit_conversion(input_df, UnitSys)
            rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df)

            st.write(' ')
            st.write(f'<span class="larger-font">Your input parameters for module <b>{function1}</b> are:</span>',
                     unsafe_allow_html=True)


            # Create a function to format the numbers to two decimal places
            def format_number(num):
                if isinstance(num, (int, float)):
                    return "{:.2f}".format(num)
                return num  # Return the original value for non-numeric elements


            # Apply the formatting function to the DataFrame
            formatted_output = output_print.applymap(format_number)
            if function1 == 'Maximum Storage Capacity':
                st.table(formatted_output.iloc[0:2, 0:9])
                st.table(formatted_output.iloc[0:2, 9:18])
                st.table(formatted_output.iloc[0:2, 18:])
            else:
                st.table(formatted_output.iloc[0:2, 0:14])
                st.table(formatted_output.iloc[0:2, 17:])
                
            check_parameter_ranges(output_df, UnitSys)

            ConRate = 0

            nWXmax = simprop.nWXmax
            nWXmin = simprop.nWXmin
            nWE = simprop.nWE

            #        Add a button to calculate the NPV array
            st.write(larger_font_style, unsafe_allow_html=True)
            st.write('<span class="larger-font">Please click "Clear Results" before a new run.</span>',
                     unsafe_allow_html=True)

            # Generate the clear/restart button
            if st.button("Clear Results"):
                reset_app()

            if "Run EASiTool" not in st.session_state:
                st.session_state["Run EASiTool"] = False

            if "Run Sensitivity Analysis" not in st.session_state:
                st.session_state["Run Sensitivity Analysis"] = False

            #if "Generate GIS" not in st.session_state:
            #    st.session_state["Generate GIS"] = False

            if st.button("Run EASiTool"):
                st.session_state["Run EASiTool"] = not st.session_state["Run EASiTool"]

            if st.session_state["Run EASiTool"]:
                # Calculates results and display output figures.
                #if function1 == 'Uniform Injection/Extraction Rate':
                #    Psave2, Qsave2, rL1, Capacity, NPV = Scenario1(rsvrprop, simprop, relakprop, npvprop, BCValue)
                #    ConRate = 1

                if function1 == 'Maximum Storage Capacity':
                    #if output_df['simmodule'][1] == 1:
                    #    st.warning('Please check your simulation module selection in your excel input file', icon="⚠️")
                    #    st.stop()
                    Psave2, Qsave2, rL1, rT1, Capacity, NPV, Qsave_res = Scenario2(rsvrprop, simprop, relakprop, npvprop, BCValue, Sensi)
                    ConRate = 0
                elif function1 == 'User Given Inputs':
                    #if output_df['simmodule'][1] == 2:
                    #    st.warning('Please check your simulation module selection in your excel input file', icon="⚠️")
                    #    st.stop()
                    Psave2_output, Qsave2_output, Capacity_S3, NPV_S3, fig_contour_geo, fig_plume_geo, highlight_locations, shift_x, shift_y, df_aor = Scenario3(
                        rsvrprop, simprop, relakprop, npvprop, BCValue, uploaded_file, UnitSys, output_print, proj_name)


                else:
                    st.write('Invalid function selected')
                    NPV = None
                    Capacity = None

                if function1 == 'Maximum Storage Capacity':
                    if NPV is None or not all(isinstance(x, (int, float)) for x in NPV):
                        st.write('Invalid NPV array')
                    if Capacity is None or not all(isinstance(x, (int, float)) for x in Capacity):
                        st.write('Invalid Capacity array')
                    else:
                        tab1_1, tab1_2 = st.tabs(['Results', 'Sensitivity Analysis'])
                        with tab1_1:
                            tab1_1_1, tab1_1_2 = st.tabs(['Capacity & NPV vs Number of Injection Wells', '2D Figures'])
                            with tab1_1_1:
                                npv_table, capacity_table = plot_npv_capacity(simprop, Capacity, NPV)
                            with tab1_1_2:
                                # Extract data from tables
                                #NPV = npv_table.iloc[:]['NPV $M']
                                #capacity = capacity_table.iloc[:]['Capacity MMT of CO₂']
                                ctrl_col1, ctrl_col2 = st.columns([0.52, 0.48])

                                with ctrl_col1:
                                    # Generate the options for the select box
                                    options = [i ** 2 for i in range(nWXmin, nWXmax + 1)]
                                    # Find the row index where 'NPV $M' reaches its maximum and make it the default selection
                                    max_npv_index = npv_table['NPV $M'].idxmax()
                                    simprop. nWX_npv_max = npv_table.loc[max_npv_index, 'nWX']
                                    default_slct= simprop.nWX_npv_max ** 2
                                    # Get the selected value in result figures
                                    selected_value = st.selectbox("Select Number of Injection Wells in the Result Figures:",
                                                                  options, index=options.index(default_slct))

                                    #st.write("Selected Number of Injection wells:", selected_value)
                                    # Compute the square root of the selected value
                                    nWX_ctrl = round(int(selected_value ** 0.5))

                                with ctrl_col2:
                                    st.write(' ')
                                    highlight_locations, fig_contour, Pnode, _ = S1S2_PContour(rsvrprop, simprop, relakprop, npvprop, BCValue, nWX_ctrl, UnitSys, Psave2, Qsave_res, ConRate)

                                results_output(output_print, BCValue, ConRate, NPV, Capacity, simprop, rsvrprop,
                                               UnitSys, Qsave2, Psave2, proj_name, Pnode, rL1)
                                output_file_path = 'EASiToolOutput.txt'
                                if output_file_path is not None:
                                    st.markdown(download_output(output_file_path), unsafe_allow_html=True)

                                ### Plume Extension and Well Property plots
                                Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
                                fig_plume = plot_plume_extension(nWX_ctrl, Xwell, Ywell, rL1, nWE, rsvrprop, UnitSys)
                                fig_well = plot_well_property(Xwell, Ywell, nWX_ctrl, nWE, ConRate, Qsave2, Psave2,
                                                              UnitSys,rsvrprop)


                                st.subheader('Output Figures')
                                fig_col1, fig_col2 = st.columns([0.52, 0.48])
                                if ConRate == 0:
                                    with fig_col1:
                                        st.pyplot(fig_contour)
                                        st.pyplot(fig_plume)
                                        st.markdown(
                                            f"<p class='centered larger-font'>Diameter of red circles represent the actual CO₂ plume size. </p>",
                                            unsafe_allow_html=True)

                                    with fig_col2:
                                        st.pyplot(fig_well)
                                        if nWE > 0:
                                            st.markdown(
                                                f"<p class='centered larger-font'>Red dots: Extraction wells </p>",
                                                unsafe_allow_html=True)
                                        st.write(' ')
                                        st.write(' ')
                                        warning_pressure_plume(nWXmin, nWXmax, rL1, Xwell, Ywell, Rwell, ConRate, Psave2,
                                                               Qsave2, rsvrprop, nWE, nWX_ctrl)

                                elif ConRate == 1:
                                    with fig_col1:
                                        #fig_contour, X3 = plot_AOR(Xwell, Ywell, Psave2, nWX_ctrl, rsvrprop, simprop, UnitSys)
                                        st.pyplot(fig_contour)
                                        st.pyplot(fig_plume)

                                    with fig_col2:
                                        st.pyplot(fig_well)
                                        st.write(' ')
                                        st.write(' ')
                                        warning_pressure_plume(nWXmin, nWXmax, rL1, Xwell, Ywell, Rwell, ConRate, Psave2,
                                                               Qsave2, rsvrprop, nWE, nWX_ctrl)
                        with tab1_2:
                            if st.button("Run Sensitivity Analysis"):
                                st.session_state["Run Sensitivity Analysis"] = not st.session_state["Run Sensitivity Analysis"]
                                geometry = 0
                                with st.spinner('Please wait...'):
                                    # Create a progress bar widget
                                    # progress_bar = st.progress(0)
                                    fig_col1, fig_col2 = st.columns(2)
                                    fig_tor, fig_multi, _, _, _, _, _, _ = sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop, BCValue, output_print, UnitSys, geometry, uploaded_file, proj_name)
                                    with fig_col1:
                                        output_file_path = 'EASiToolOutput_Sensitivity.txt'
                                        if output_file_path is not None:
                                            st.markdown(download_SA_output(output_file_path), unsafe_allow_html=True)
                                        st.pyplot(fig_tor)
                                    with fig_col2:
                                        st.write('')
                                        st.write('')
                                        st.write('')
                                        st.pyplot(fig_multi)
                                st.success('Sensitivity analysis (tornado chart) completed.')

                                selection_col1, selection_col2 = st.columns(2)
                                with selection_col1:
                                    st.markdown(
                                        '<h3 style="margin-bottom: 0;">Change number of samples used in the Monte Carlo Simulation:</h3>',
                                        unsafe_allow_html=True)
                                    function7 = st.selectbox(
                                        ' ',
                                        ('100, <1 minute', '500, <3 minutes')
                                    )
                                    drawdown_options_num_rel = ['100, <1 minute', '500, <3 minutes']
                                with selection_col2:
                                    pass

                                # sensitivity_parameters = ['P0', 'temp', 'Thickness', 'salinity', 'Porosity', 'k', 'cr', 'FracP', 'Sar',
                                #                          'Sgc', 'm', 'n', 'kar0', 'krg0']
                                sensitivity_parameters_num_rel = [100, 500]
                                # Sensi_para = sensitivity_parameters[drawdown_options.index(function6)]
                                Sensi_para_rel_num = sensitivity_parameters_num_rel[drawdown_options_num_rel.index(function7)]

                                with st.spinner('Please wait...'):
                                    # Create a progress bar widget
                                    # progress_bar = st.progress(0)
                                    rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df) # Initialized the properties again
                                    hist_data = generate_histdata(Sensi_para_rel_num, rsvrprop, simprop, relakprop, npvprop,
                                                                  BCValue, uploaded_file, UnitSys, geometry, output_print)

                                    figure_hist = plot_density_curve(hist_data, geometry, UnitSys)
                                    # Reset the progress bar at the end
                                    fig2_col1, fig2_col2 = st.columns(2)
                                    with fig2_col1:
                                        st.pyplot(figure_hist)
                                    with fig2_col2:
                                        pass
                                st.success(
                                    'Sensitivity analysis completed. Random sampling of parameters space for a broader exploration of the sensitivity. '
                                    'Standard deviation is 1% of the mean value of the input parameters.')

                elif function1 == 'User Given Inputs':
                    tab2_1, tab2_3, tab2_2 = st.tabs(['Results', 'GIS Map','Sensitivity Analysis'])
                    with tab2_1:
                        output_file_path = 'EASiToolOutput_Geometry.txt'
                        if output_file_path is not None:
                            st.markdown(download_output(output_file_path), unsafe_allow_html=True)

                        st.subheader('Output Figures')
                        fig_col1, fig_col2 = st.columns([0.52, 0.48])
                        with fig_col1:
                            st.pyplot(fig_contour_geo)
                        with fig_col2:
                            st.pyplot(fig_plume_geo)
                        st.subheader('Summary of Results')
                        st.markdown(
                            f"<p class='larger-font'>NPV: {round(NPV_S3, 2)} $M, Storage Capacity: {round(Capacity_S3, 2)} MMT of CO₂</p>",
                            unsafe_allow_html=True)
                    with tab2_2:
                        if st.button("Run Sensitivity Analysis"):
                            st.session_state["Run Sensitivity Analysis"] = not st.session_state[
                                "Run Sensitivity Analysis"]
                            geometry = 1
                            with st.spinner('Please wait...'):
                                # Create a progress bar widget
                                # progress_bar = st.progress(0)
                                fig_col1, fig_col2 = st.columns(2)
                                fig_tor, fig_multi, _, _, _, _, _, _ = sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop, BCValue, output_print, UnitSys,
                                                              geometry, uploaded_file, proj_name)
                                with fig_col1:
                                    output_file_path = 'EASiToolOutput_Sensitivity.txt'
                                    if output_file_path is not None:
                                        st.markdown(download_SA_output(output_file_path), unsafe_allow_html=True)
                                    st.pyplot(fig_tor)
                                with fig_col2:
                                    st.write('')
                                    st.write('')
                                    st.write('')
                                    st.pyplot(fig_multi)
                            st.success('Sensitivity analysis (tornado chart) completed.')

                            selection_col1, selection_col2 = st.columns(2)
                            with selection_col1:
                                st.markdown(
                                    '<h3 style="margin-bottom: 0;">Change number of samples used in the Monte Carlo Simulation:</h3>',
                                    unsafe_allow_html=True)
                                function7 = st.selectbox(
                                    ' ',
                                    ('100, <1 minute', '500, <3 minutes', '1000, >5 minutes')
                                )
                                drawdown_options_num_rel = ['100, <1 minute', '500, <3 minutes', '1000, >5 minutes']
                            with selection_col2:
                                pass

                            # sensitivity_parameters = ['P0', 'temp', 'Thickness', 'salinity', 'Porosity', 'k', 'cr', 'FracP', 'Sar',
                            #                          'Sgc', 'm', 'n', 'kar0', 'krg0']
                            sensitivity_parameters_num_rel = [100, 500, 1000]
                            # Sensi_para = sensitivity_parameters[drawdown_options.index(function6)]
                            Sensi_para_rel_num = sensitivity_parameters_num_rel[drawdown_options_num_rel.index(function7)]

                            with st.spinner('Please wait...'):
                                # Create a progress bar widget
                                # progress_bar = st.progress(0)
                                rsvrprop, simprop, relakprop, npvprop = read_rsvr_input(output_df) # Initialized the properties again
                                hist_data = generate_histdata(Sensi_para_rel_num, rsvrprop, simprop, relakprop, npvprop, BCValue,
                                                              uploaded_file, UnitSys, geometry, output_print)
                                figure_hist = plot_density_curve(hist_data, geometry, UnitSys)
                                # Reset the progress bar at the end
                                fig2_col1, fig2_col2 = st.columns(2)
                                with fig2_col1:
                                    st.pyplot(figure_hist)
                                with fig2_col2:
                                    pass
                            st.success(
                                'Sensitivity analysis completed. Random sampling of parameters space for a broader exploration of the sensitivity. '
                                'Standard deviation is 1% of the mean value of the input parameters.')
                    with tab2_3:
                        ### GIS map
                        st.write(
                            "You can also generate GIS map for the reservoir and plume based on inputs' UTM coordinates, and the typed zone number and zone letter.")
                        st.write(
                            "https://www.dmap.co.uk/utmworld.htm and https://www.latlong.net/lat-long-utm.html \n"
                            "are good websites to find the UTM zone, and to convert UTM coordinates from your latitude and longitude."
                        )
                        st.subheader("GIS Map:")

                        # Create two columns for input fields
                        GIS_col1_1, GIS_col1_2 = st.columns(2)

                        # Input field for zone number
                        with GIS_col1_1:
                            zone_number = st.number_input("Enter UTM zone number:", value=14, step=1)
                        # Input field for zone letter
                        with GIS_col1_2:
                            pass

                        #if st.button("Generate GIS"):
                        #    st.session_state["Generate GIS"] = not st.session_state["Generate GIS"]
                        with st.spinner('Please wait...'):
                            ### GIS map based on UTM
                            geometry = 1
                            _, _, _, _, df_aor_up, df_aor_low, X_center, Y_center = sensi_ana_tor(rsvrprop, simprop, relakprop, npvprop,
                                                                              BCValue,
                                                                              output_print, UnitSys,
                                                                              geometry, uploaded_file, proj_name)
                            map = GIS_map(uploaded_file, df_aor, df_aor_up, df_aor_low, zone_number, UnitSys, X_center, Y_center, rsvrprop)
                            # Render the map as an HTML string
                            map_html = map._repr_html_()
                            map.width = '100%'
                            map.height = '600px'
                            st.components.v1.html(map_html, width=1350, height=1000, scrolling=True)



if __name__ == "__main__":
    main()
