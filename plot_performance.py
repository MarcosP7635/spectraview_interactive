import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

st.title("Interactive Web Application to Plot the Peformance of Drones and Cubesats Using AWG Photonic Chips")
number = st.number_input("Please enter the altitude of the Cubesat Orbit in kilometers")
st.write("Entered cubesat orbit altitude in km: ", number)
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