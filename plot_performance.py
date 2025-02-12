import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import astropy.constants as const
import astropy.units as unit

st.title("Interactive Web Application to Plot the Peformance of Drones and Cubesats Using AWG Photonic Chips")
a_surface = st.number_input("Please enter the altitude of the Cubesat Orbit in kilometers") * unit.km
st.write("Entered cubesat orbit altitude: ", a_surface)
cubesat_aperture = st.number_input("Please enter the aperture of the primary optic of the cubesat in cm") * unit.cm
st.write("Entered cubesat orbit altitude: ", cubesat_aperture)
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
p = 'spectralcalc_transmission/extraterrestrial_solar_spectrum_wavelengths.txt'
solar_microns, all_floats = txt_to_float_arr(p)
p = 'spectralcalc_transmission/extraterrestrial_solar_irradiance.txt'
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
    with open(p, 'r') as file:
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
smallest_neighbors = [get_smallest_neighbor(c, scale_1_wavelen_arr) for c in channel_edges]
""" 
Now we will make dictionaries. The first dictionary is callend wavelen_arrs_dict 
and each key is a scale_factor and each value is be the corresponding wavelengths used in the 
output of the spectralcalc simulation. The second dictionary is callend transmi_arrs_dict and each 
key is a scale_factor and each value is the corresponding transmission as a dimensionless number 
between 0 and 1 used in the output of the spectralcalc simulation. 
"""
wavelen_arrs_dict, transmi_arrs_dict, equal_lengths = {}, {}, {}
for i in range(len(scale_factors)):
    sf = scale_factors[i]
    wavelen_arrs_dict[sf], transmi_arrs_dict[sf] = read_spectralcalc_output(file_paths[i])
    equal_lengths[sf] = len(wavelen_arrs_dict[sf]) == len(transmi_arrs_dict[sf])
def get_channel_avg_transmi(transmi_arr):
    """Returns the average transmittance in a spectral channel"""
    return [np.average(transmi_arr[smallest_neighbors[i][0]:smallest_neighbors[i+1][0]]) 
    for i in range(20)]
awg_channels_avg_transmi_arrs_dict = {}
for sf in scale_factors:
    awg_channels_avg_transmi_arrs_dict[sf] = get_channel_avg_transmi(
        transmi_arrs_dict[sf])
def get_avg_awg_radiance(awg_channels_avg_transmi, solar_radiances):
    return [ (0.05/a_surface) *
        awg_channels_avg_transmi[i] * solar_radiances[solar_smallest_neighbors[i][0]] 
        for i in range(len(awg_channels_avg_transmi)) 
    ]
awg_channels_avg_irr_arrs_dict = {}
for sf in scale_factors:
    awg_channels_avg_irr_arrs_dict[sf] = get_avg_awg_radiance(
        awg_channels_avg_transmi_arrs_dict[sf], solar_radiances)
M, G = const.M_earth, const.G
a_center_earth = (a_surface + const.R_earth).to(unit.m)
P = (np.sqrt(4 * np.pi**2 * a_center_earth**3 / (G*M))).to(unit.second)
v = (np.sqrt(G*M/a_center_earth)).to(unit.m/unit.second)
orbit_fraction = v * 0.1 * unit.second / (2 * np.pi * a_center_earth) 
distance_swept_on_Earth = orbit_fraction * const.R_earth * 2 * np.pi #meters
fov_radius_meters = (1.22 * (1.5*10**-6 / 5 *10**-2)) * a_surface.to(unit.m)  
radians_viewed = 1.22 * (1.5*10**-6 / 5 *10**-2) #radians!
solid_angle_viewed = radians_viewed**2 #steradians
photon_energy = 10**-34 * 6.626 * 299792458 / (1560 * 10**-9)
area = cubesat_aperture.to(unit.m)**2 * np.pi #square meters
solid_angle_of_telescope_viewed_by_Earth = float((area / a_surface**2).to(""))
def w_per_m2_per_sr_to_counts_per_second(w_per_m2_per_sr_arr):
    return [x * area * solid_angle_of_telescope_viewed_by_Earth / photon_energy
        for x in w_per_m2_per_sr_arr]
counts_per_second_arrs_dict = {}
for sf in scale_factors:
    counts_per_second_arrs_dict[sf] = w_per_m2_per_sr_to_counts_per_second(
        awg_channels_avg_irr_arrs_dict[sf])



