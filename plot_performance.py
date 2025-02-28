import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sympy

st.title(
"Interactive Web Application to Plot the Peformance of Drones and Cubesats Using AWG Photonic Chips")
a_surface_km = st.number_input(
"Please enter the altitude of the Cubesat Orbit in kilometers",
value=500.0, placeholder="Type a number...") 
st.write("Entered cubesat orbit altitude: ", str(a_surface_km), " km")
cubesat_aperture_cm = st.number_input(
"Please enter the aperture of the primary optic of the cubesat in cm",
value=10.0, placeholder="Type a number...")
cubesat_aperture_m = cubesat_aperture_cm * 10**-2 
st.write("Entered cubesat aperture: ", str(cubesat_aperture_cm), " cm")
time_exposure = st.number_input(
"Please enter the exposure time of the Cubesat's Detector in seconds",
value=0.1, placeholder="Type a number...") 
st.write("Entered cubesat exposure time: ", str(time_exposure), " seconds")
#now convert to kilograms, meters, and seconds SI.
cubesat_aperture, a_surface = cubesat_aperture_cm * 10**-2, a_surface_km * 10**3
D, lam = sympy.symbols("D"), sympy.symbols("lambda")
lam_symbol, lam_val = sympy.symbols("lambda"), 1.662 * 10**-6
lam_def = "the center wavelength, in meters, of the AWG photonic chip"
D_symbol, D_val = sympy.symbols("D"), cubesat_aperture_m
D_def = "the diameter of the telescope's primary optic in meters"
vars_dict = {} 
vars_dict[lam], vars_dict[D] = {"Definition": lam_def, "Value": lam_val}, {"Definition": D_def, "Value": D_val}
rayleigh_criterion = 1.22 * lam / D
vars_dict[lam]["Symbol"], vars_dict[D]["Symbol"] = lam_symbol, D_symbol
#Calculate the rayleigh criterion as a float
rayleigh_criterion_val = rayleigh_criterion.subs([(x, vars_dict[x]["Value"])
    for x in rayleigh_criterion.free_symbols])
rayleigh_criterion_symbol = sympy.symbols("theta")
vars_dict[rayleigh_criterion] = {"Value": rayleigh_criterion_val, 
                                 "Symbol": rayleigh_criterion_symbol,
"Definition": "the angular resolution in radians of the telescope on the cubesat perpendicular to the motion of the cubesat"}
#print the result
radians_viewed = float(rayleigh_criterion_val)
fov_radius_meters = radians_viewed * a_surface / 2
M = 5.9721679 * 10**24 #kg
G = 6.6743 * 10**-11 #m^3 kg^-1 s^-2
R_earth = 6378100 #m
a_center_earth = a_surface + R_earth
P = np.sqrt(4 * np.pi**2 * a_center_earth**3 / (G*M))
v = np.sqrt(G*M/a_center_earth) #meters per second
orbit_fraction = time_exposure / P #meters
distance_swept_on_Earth = orbit_fraction * R_earth * 2 * np.pi
solid_angle_viewed = radians_viewed**2 #steradians
area = ((cubesat_aperture/2)**2) * np.pi #square meters
solid_angle_of_telescope_viewed_by_Earth = np.pi * (cubesat_aperture/2)**2 / (a_surface**2) #steradians
side_length_pixel_perp_to_orbit_dir = 2 * fov_radius_meters 
st.write("The side length of the pixel in the direction perpendicular to the orbit direction is " + str(np.round(
    side_length_pixel_perp_to_orbit_dir)) + " meters.")
st.write("The side length of the pixel in the direction parallel to the orbit direction is " + str(np.round(
    distance_swept_on_Earth, 1)) + " meters.")
st.write("So the pixel size is " + str(np.round(
    side_length_pixel_perp_to_orbit_dir)) + " meters by " +
    str(np.round(distance_swept_on_Earth, 1)) + " meters.")
def get_vars(expression):
    return [x for x in expression.free_symbols]
if st.button("Show work?"):
    vars = get_vars(rayleigh_criterion)
    vars.append(rayleigh_criterion)
    st.write("$$\\Large " + sympy.latex(vars_dict[rayleigh_criterion]["Symbol"])
     + " = " + sympy.latex(rayleigh_criterion) + "$$")
    for x in vars:
        val_str = str(vars_dict[x]["Value"])
        if "e" in str(vars_dict[x]["Value"]):
            floats, oom = val_str.split("e")
            floats = float(floats)
            st.write("$ "+ sympy.latex(vars_dict[x]["Symbol"]), " $ is ", vars_dict[x]["Definition"], 
                 " and has a value of $", str(np.round(floats,4)), "\\times 10^{", 
                 str(int(oom)), "}$")
        else:
            val_str = str(np.round(vars_dict[x]["Value"],2))
            st.write("$ "+ sympy.latex(x), "$ is ", vars_dict[x]["Definition"], 
                 " and has a value of $", val_str, "$")
            val_str = str(vars_dict[x]["Value"])
scale_factors = [2**n for n in range(10)]
for n in range(16):
    scale_factors.append(1.5**(n+1))
for x in [1.001, 1.005, 1.01, 1.05, 1.1, 1.25, 1.75]:
    scale_factors.append(x)
scale_factors.sort()
def get_file_path(parent_dir, scale_factor, name_start, file_type):
    """
    Generate a file path string based on the provided parameters.

    Parameters
    ----------
    parent_dir : str
        The parent directory of the file.
    scale_factor : float
        The scale factor to be converted to a string.
    name_start : str
        The starting string before the scale factor.
    file_type : str
        The file type to append to the end of the string.

    Returns
    -------
    str
        A string representing the file path.
    """
    return (parent_dir + name_start +
            str(scale_factor).replace(".", "point") + file_type)
parent_dir = 'spectralcalc_transmission/'
file_paths = [get_file_path(parent_dir, scale_factor, "cubesat_scale", ".txt") 
    for scale_factor in scale_factors]

def txt_to_float_arr(path, first_index=1):
    """
    Convert a text file to a list of floats. We only use this to import the extraterrestrial
    solar spectrum.

    Parameters
    ----------
    path : str
        The path to the text file.
    first_index : int
        The index of the first element to be converted to a float.

    Returns
    -------
    list
        A list of floats.
    """
    with open(path, 'r') as file:
        output_arr = [float(x) for x in file.read().split('\n')[first_index:]]
        all_floats = all([type(x) == float for x in output_arr])
        return output_arr, all_floats
p = 'spectralcalc_transmission/extraterrestrial_solar_spectrum_wavelengths.txt' #microns
solar_microns, all_floats = txt_to_float_arr(p)
p = 'spectralcalc_transmission/extraterrestrial_solar_irradiance.txt' #W/m2/micron
solar_radiances, all_floats = txt_to_float_arr(p)
def read_spectralcalc_output(path, start_line=25, end_line=1):
    """
    Reads the spectral calculation output from a file, extracting wavelength
    and transmittance data.

    Parameters
    ----------
    path : str
        The path to the file containing spectral data.
    start_line : int, optional
        The number of initial lines to skip before reading data (default is 25).
    end_line : int, optional
        The number of lines at the end of the file to exclude from reading (default is 1).

    Returns
    -------
    tuple
        A tuple containing two lists:
        - wavelen_arr: List of wavelengths in microns extracted from the file.
        - transmi_arr: List of transmittance values (dimensionless fraction between 0 and 1) 
            corresponding to the wavelengths.
    """
    with open(path, 'r') as file:
        file_lines_arr = file.read().split('\n')
    wavelen_arr, transmi_arr = [], []
    for i in range(len(file_lines_arr)-(start_line+end_line)):
        wavelength, transmittance = [float(x) for x in 
                                     file_lines_arr[i+start_line].split(" ")]
        wavelen_arr.append(wavelength)
        transmi_arr.append(transmittance)
    return wavelen_arr, transmi_arr

#Now we need arrays describing the boundaries and centers of each spectral channel in the AWG
channel_edges, channel_centers = np.linspace(1665.15, 1666.25, 21), np.linspace(1665.2, 1666.2, 20)
def get_smallest_neighbor(channel_wavelen, wavelen_arr): 
    """Returns the index and wavelength of the nearest neighbor 
    of channel_wavelen to the channel_wavelength in wavelen_arr"""
    for i in range(len(wavelen_arr)):
        if wavelen_arr[i] > channel_wavelen/1000:
            return (i, wavelen_arr[i])
solar_smallest_neighbors = [get_smallest_neighbor(c, solar_microns) for c in channel_edges]
p = file_paths[0]
scale_1_wavelen_arr, scale_1_transmi_arr = read_spectralcalc_output(p)
smallest_neighbors = [get_smallest_neighbor(c, scale_1_wavelen_arr) for c in channel_edges]
wavelen_arrs_dict, transmi_arrs_dict, equal_lengths = {}, {}, {}
#Now we will make dictionaries. The first dictionary is callend wavelen_arrs_dict 
#and each key is a scale_factor and each value is be the corresponding wavelengths used in the 
#output of the spectralcalc simulation. The second dictionary is callend transmi_arrs_dict and each 
#key is a scale_factor and each value is the corresponding transmission as a dimensionless number 
#between 0 and 1 used in the output of the spectralcalc simulation. 
for i in range(len(scale_factors)):
    sf = scale_factors[i]
    wavelen_arrs_dict[sf], transmi_arrs_dict[sf] = read_spectralcalc_output(file_paths[i])
    equal_lengths[sf] = len(wavelen_arrs_dict[sf]) == len(transmi_arrs_dict[sf])
def get_channel_avg_transmi_squared(transmi_arr_input):
    return [np.average(transmi_arr_input[smallest_neighbors[i][0]:
    smallest_neighbors[i+1][0]]) **2 for i in range(20)]
awg_channels_avg_transmi_arrs_dict = {}
for sf in scale_factors:
    awg_channels_avg_transmi_arrs_dict[sf] = (
    get_channel_avg_transmi_squared(transmi_arrs_dict[sf]))
def get_avg_awg_radiance(awg_channels_avg_transmi, solar_radiances):
    #that coefficient is because the spectral width of each spectral channel is 0.05 nm
    #and the solar radiances are given in W/m2/micron
    return [ (0.05/1000) *
        awg_channels_avg_transmi[i] * solar_radiances[solar_smallest_neighbors[i][0]] 
        for i in range(len(awg_channels_avg_transmi)) 
    ]
awg_channels_avg_irr_arrs_dict = {}
for sf in scale_factors:
    awg_channels_avg_irr_arrs_dict[sf] = get_avg_awg_radiance(
        awg_channels_avg_transmi_arrs_dict[sf], solar_radiances)
photon_energy = 10**-34 * 6.626 * 299792458 / (1560 * 10**-9) #joules
#resolveable distance in meters 

st.write("Thus the pixel size is " + str(np.round(
    side_length_pixel_perp_to_orbit_dir)) + " meters by " +
    str(np.round(distance_swept_on_Earth, 1)) + " meters.")

def w_per_m2_to_counts_per_second(w_per_m2_arr):
    """
    Convert a list of average irradance values from W/m2 to a list of 
    average counts per second for each spectral channel of the AWG 
    photonic chip.

    Parameters
    ----------
    w_per_m2_arr : list
        Average irradance values in W/m2 for each spectral channel of the AWG 
        photonic chip.

    Returns
    -------
    list
        Average counts per second for each spectral channel of the AWG 
        photonic chip.

    Need to account for reflectivity of the ground.
    """
    return [x * (fov_radius_meters**2) * 
    solid_angle_of_telescope_viewed_by_Earth / (2 * np.pi * photon_energy)
    for x in w_per_m2_arr]
counts_per_second_arrs_dict = {}
for sf in scale_factors:
    counts_per_second_arrs_dict[sf] = w_per_m2_to_counts_per_second(
        awg_channels_avg_irr_arrs_dict[sf])
all_scales_counts_per_second = [x for x in counts_per_second_arrs_dict.values()]

def is_detectable(index):
    st.write("So the total difference in the number of counts per second between background methane and " 
        + str((np.round((scale_factors[index]-1)*100, 3))) + "% more methane is" )
    count_difference = time_exposure * (
        np.sum(counts_per_second_arrs_dict[scale_factors[0]]) - 
        np.sum(counts_per_second_arrs_dict[scale_factors[index]]))
    photon_noise_background_measurement = np.sqrt(time_exposure * 
        np.sum(counts_per_second_arrs_dict[scale_factors[0]]))
    photon_noise_not_from_background = np.sqrt(time_exposure *
        np.sum(counts_per_second_arrs_dict[scale_factors[index]]))
    photon_noise = np.sqrt(photon_noise_background_measurement**2 + 
        photon_noise_not_from_background**2)
    #error propagation since these would be independent measurements
    st.write("Difference in counts: ", count_difference)
    st.write("Photon Noise in counts: ", photon_noise)
    sensitive = (count_difference > photon_noise)
    if sensitive:
        st.write("We are sensitive to that difference in counts!")
    else:
        st.write("The difference in counts is less than the photon noise.")
    output_dict = {
        "count_difference": count_difference,
        "photon_noise": photon_noise,
        "detectable": sensitive
    }
    return output_dict
start_index = 1
for i in range(start_index, len(scale_factors)): 
    stats_dict = is_detectable(i)
    if stats_dict["detectable"]:
        break
background_methane_moles_per_square_meter = 0.7
excess_methane_moles_per_square_meter = ((scale_factors[i] - 1) * 
background_methane_moles_per_square_meter)
avg_wind_speed_meters_per_second = 3.5
time_to_move_methane = (side_length_pixel_perp_to_orbit_dir / 
    avg_wind_speed_meters_per_second)
max_approx_flow_rate_leak_moles_per_second = (excess_methane_moles_per_square_meter / 
    time_to_move_methane)
min_leak_detected_kg_per_hour = max_approx_flow_rate_leak_moles_per_second * 0.016 * 3600
#Difference in photon counts plus or minus photon noise
uncertainty_kg_per_hour = (min_leak_detected_kg_per_hour * stats_dict["photon_noise"] 
                           / stats_dict["count_difference"])
st.write("That would be a leak with a flow rate less than " + 
str(min_leak_detected_kg_per_hour) + "$ \\pm $ " + 
str(uncertainty_kg_per_hour) + " kg/h")