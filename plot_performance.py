import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sympy
import json

st.title(
"Interactive Web Application to Plot the Peformance of Drones and Cubesats Using AWG Photonic Chips")
cost_annualized_direct = 1000*st.number_input("Please enter the annual budget in thousands of dollars (USD)", 
    value=604, placeholder="Type a number...", step=604, min_value=604)
st.write("Entered annual budget: " + str(int(np.floor(cost_annualized_direct/1000))) + "k USD")
satellites_available = int(np.floor(cost_annualized_direct / 604000))
st.write("This budget allows for " + str(satellites_available) + 
         " satellites, which with a polar orbit allows for a chosen site to "
         + "be viewed every " +  str(12 / satellites_available) + " hours.")
satellite_number = st.number_input("Alternatively, you can enter the number of satellites here ", 
                value=1, placeholder="Type a number...", step=1, min_value=1)
st.write("Which would cost " + str(604 * satellite_number) + "k USD.")
st.write("Using a polar orbit with a period of 90 minutes, a chosen site on Earth can be visited"
      +   " every " + str(np.round(43200 / (3600 * satellite_number), 2)) + " hours.")
a_surface_km = st.number_input(
"Please enter the altitude of the Cubesat Orbit in kilometers",
value=500, placeholder="Type a number...", step=25) 
st.write("Entered cubesat orbit altitude: ", str(a_surface_km), " km")
cubesat_aperture_cm = st.number_input(
"Please enter the aperture of the primary optic of the cubesat in cm",
value=5.0, placeholder="Type a number...", step=1.0)
cubesat_aperture_m = cubesat_aperture_cm * 10**-2 
st.write("Entered cubesat aperture: ", str(cubesat_aperture_cm), " cm")
time_exposure_ms = st.number_input(
"Please enter the exposure time of the Cubesat's Detector in milliseconds",
value=100, placeholder="Type a number...") 
st.write("Entered cubesat exposure time: ", str(time_exposure_ms), " milliseconds")
#now convert to kilograms, meters, and seconds SI.
cubesat_aperture, a_surface = cubesat_aperture_cm * 10**-2, a_surface_km * 10**3
time_exposure = time_exposure_ms * 10**-3
def import_dict_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)
unconverted = import_dict_json("awg_channels_avg_irr_arrs_dict.json")
converted = {}
for key, value in unconverted.items():
    float_key = float(key)
    converted[float_key] = [np.float64(x) for x in value] 
awg_channels_avg_irr_arrs_dict = converted

scale_factors = [2**n for n in range(10)]
for n in range(16):
    scale_factors.append(1.5**(n+1))
for x in [1.001, 1.005, 1.01, 1.05, 1.1, 1.25, 1.75]:
    scale_factors.append(x)
scale_factors.sort()
    
D, lam, r = sympy.symbols("D"), sympy.symbols("lambda"), sympy.symbols("r")
lam_symbol, lam_val = sympy.symbols("lambda"), 1.662 * 10**-6
lam_def = "the center wavelength, in meters, of the AWG photonic chip"
D_symbol, D_val = sympy.symbols("D"), cubesat_aperture_m
D_def = "the diameter of the telescope's primary optic in meters"
r_symbol, r_val = sympy.symbols("r"), a_surface
r_def = "the altitude of the orbit of the cubesat in meters"
vars_dict = {} 
vars_dict[lam], vars_dict[D] = {"Definition": lam_def, "Value": lam_val}, {"Definition": D_def, 
                                                                            "Value": D_val}
vars_dict[r] = {"Definition": r_def, "Value": r_val}
resolved_distance = r * 1.22 * lam / D
vars_dict[lam]["Symbol"], vars_dict[D]["Symbol"], vars_dict[r]["Symbol"] = lam_symbol, D_symbol, r_symbol
#Calculate the rayleigh criterion as a float
resolved_distance_val = np.float64(resolved_distance.subs([(x, vars_dict[x]["Value"])
    for x in resolved_distance.free_symbols], numpy=True))
resolved_distance_symbol = sympy.symbols("L")
vars_dict[resolved_distance] = {"Value": resolved_distance_val, 
                                "Symbol": resolved_distance_symbol,
"Definition": ("the angular resolution in radians of the telescope on the cubesat perpendicular to the " +
    " motion of the cubesat")}
#print the result
radians_viewed = float(resolved_distance_val) / a_surface
fov_radius_meters = resolved_distance_val / 2
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
side_length_pixel_perp_to_orbit_dir = np.float64(resolved_distance_val)
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
    vars = get_vars(resolved_distance)
    vars.append(resolved_distance)
    st.write("$$\\Large " + sympy.latex(vars_dict[resolved_distance]["Symbol"])
     + " = " + sympy.latex(resolved_distance) + "$$")
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

photon_energy = 10**-34 * 6.626 * 299792458 / (1560 * 10**-9) #joules
#resolveable distance in meters 
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
def calc_min_leak_kg_per_hour(SNR, distance, wind_speed):
    #$W$ is the size of the region, $Q$ is the flow rate. Similarly, we can rearrange and use the 
    # SNR $\sigma$ to find the minimum detectable flow rate
    #\begin{align}
    #Q_{\text{min}} = \frac{M_{\text{CH}_4}UWp\sigma}{gM_\text{a}}
    #\end{align}
    #Which agrees with Equation 14 from Jacobs et al. 2016.  https://doi.org/10.5194/acp-16-14371-2016
    #They further multiply some dimensionless $q$ where "q takes on values of 2 for detection and 5 for 
    # quantification."
    #Please enter values in units of kg, hour, m, and m/hour
    #Sigma is the uncertainty or 1/SNR
    #return excess_methane * 0.7 * 0.016 * 3600 / (distance / wind_speed)
    background_methane_fraction = 1921.72 * (10**-9) #moles of methane per mole of dry air
    p = 101325 #pascals
    mCH4 = 0.016 #kg / mol
    mAir = 0.029 #kg / mol
    q_detect = 2 #dimensionless
    g = 9.81 #meters / seconds^2
    #pascals are really just kg/(m * seconds^2) , so the seconds^2 cancels from g
    #Write down the dimensions by hand and then check using sympy
    return distance * wind_speed * mCH4 * q_detect * p * background_methane_fraction / (
        SNR * mAir * g)
SNR =  np.sqrt(time_exposure * np.sum(counts_per_second_arrs_dict[scale_factors[0]]))
avg_wind_speed_meters_per_second = 3.5
avg_wind_speed_meters_per_hour = avg_wind_speed_meters_per_second * 3600 #meters per hour
min_flow_rate = calc_min_leak_kg_per_hour(SNR, side_length_pixel_perp_to_orbit_dir, 
    avg_wind_speed_meters_per_hour)
st.write("Since the Signal-to-Noise Ratio (SNR) is " + str(np.round(SNR, 2)) + 
", the smallest methane leak we can detect has a flow rate of " + 
str(np.round(min_flow_rate, 2)) + " kg/hour.")
#Fixed W! Consider that we have to consider a specific instant, not the whole exposure!   
#return [side_length_pixel_perp_to_orbit_dir, distance_swept_on_Earth, 
#    calc_min_leak_kg_per_hour(SNR, side_length_pixel_perp_to_orbit_dir, 
#                                    avg_wind_speed_meters_per_hour)]

t, n_dot = sympy.symbols("t"), sympy.symbols("\\frac{dn}{dt}")
W, U, sigma = sympy.symbols("W"), sympy.symbols("U"), sympy.sqrt(t*n_dot)
t_symbol, t_val = sympy.symbols("t"), time_exposure
t_def = "the exposure time of the detector in seconds"
n_dot_symbol, n_dot_val = sympy.symbols("(dn/dt)"), np.float64(np.sum(counts_per_second_arrs_dict[scale_factors[0]]))
n_dot_def = "the average number of photons per second incident on the detector summed over every spectral channel"
W_symbol, W_val = sympy.symbols("D"), side_length_pixel_perp_to_orbit_dir
W_def = "the width of the instantaneous field of view seen by the satellite in meters "
U_symbol, U_val = sympy.symbols("lambda"), avg_wind_speed_meters_per_hour
U_def = "the average wind speed in meters per hour"

flow_vars_dict = {} 
flow_vars_dict[t], flow_vars_dict[n_dot] = {"Definition": t_def, "Value": t_val}, {"Definition": n_dot_def, 
    "Value": n_dot_val}
flow_vars_dict[W], flow_vars_dict[U] = {"Definition": W_def, "Value": W_val}, {"Definition": U_def, 
    "Value": U_val}
flow_vars_dict[t]["Symbol"], flow_vars_dict[n_dot]["Symbol"], flow_vars_dict[W]["Symbol"], flow_vars_dict[U]["Symbol"] = t_symbol, n_dot_symbol, W_symbol, U_symbol

pressure, massCH4, massAir = sympy.symbols("p"), sympy.symbols("m_{\\text{CH}_4}"), sympy.symbols("m_\\text{air}")
g, g_val = sympy.symbols("g"), 9.81
g_def = "the acceleration due to gravity in meters per second^2"
p_def = "the atmospheric pressure on Earth's surface in pascals"
massCH4_def = "the molar mass of methane in kilograms"
massAir_def = "the molar mass of air in kilograms per mole"
p_val, massCH4_val, massAir_val = 101325, 0.016, 0.029
f, f_symbol = sympy.symbols("f"), sympy.symbols("f")
f_val = 1921.72 * (10**-9) #moles of methane per mole of dry air
f_dev = ("fraction of moles of dry air that are methane on " + 
    "average on Earth in an atmospeheric column ")

flow_vars_dict[pressure], flow_vars_dict[massCH4], flow_vars_dict[massAir] = {"Definition": p_def, "Value": p_val}, {"Definition": massCH4_def, 
    "Value": massCH4_val}, {"Definition": massAir_def, "Value": massAir_val}
pressure_symbol, massCH4_symbol, massAir_symbol = sympy.symbols("p"), sympy.symbols("m_{\\text{CH}_4}"), sympy.symbols("m_\\text{air}")
flow_vars_dict[pressure]["Symbol"], flow_vars_dict[massCH4]["Symbol"], flow_vars_dict[massAir]["Symbol"] = pressure_symbol, massCH4_symbol, massAir_symbol
flow_vars_dict[g] = {"Definition": g_def, "Value": g_val}
flow_vars_dict[g]["Symbol"] = sympy.symbols("g")
flow_vars_dict[f] = {"Definition": f_dev, "Value": f_val}
flow_vars_dict[f]["Symbol"] = sympy.symbols("f")

Q = sympy.symbols("Q_\\text{min}")
q_detect = 2 #dimensionless
Q_expr = U * W * massCH4 * q_detect * pressure * f / (
(sympy.sqrt(t*n_dot)) * massAir * g)
#Q = distance * wind_speed * mCH4 * q_detect * p * background_methane_fraction / (
#SNR * mAir * g

#Calculate the rayleigh criterion as a float
min_flow_rate_detected_val = np.float64(Q_expr.subs([(x, flow_vars_dict[x]["Value"])
    for x in Q_expr.free_symbols], numpy=True))
resolved_distance_symbol = sympy.symbols("L")
flow_vars_dict[Q] = {"Value": min_flow_rate_detected_val, 
    "Symbol": sympy.symbols("Q"),
"Definition": ("The flow rate in kilograms per hour of the smallest methane " + 
    "leak that the satellite can detect ")}


if st.button("Show work for the threshold flow rate?"):
    vars = get_vars(Q_expr)
    vars.append(Q)
    st.write("$$\\Large " + sympy.latex(Q) + " = " + sympy.latex(Q_expr) + " $$")
    for x in vars:
        val_str = str(flow_vars_dict[x]["Value"])
        if "e" in str(flow_vars_dict[x]["Value"]):
            floats, oom = val_str.split("e")
            floats = float(floats)
            st.write("$ "+ sympy.latex(flow_vars_dict[x]["Symbol"]), " $ is ", flow_vars_dict[x]["Definition"], 
                 " and has a value of $", str(np.round(floats,4)), "\\times 10^{", 
                 str(int(oom)), "}$")
        else:
            val_str = str(np.round(flow_vars_dict[x]["Value"],2))
            st.write("$ "+ sympy.latex(x), "$ is ", flow_vars_dict[x]["Definition"], 
                 " and has a value of $", val_str, "$")
            val_str = str(flow_vars_dict[x]["Value"])

