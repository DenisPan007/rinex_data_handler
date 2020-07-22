import math

import georinex as gr
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from ephemeris import Ephemeris

GM = 3.986005 * 10e14
PI = math.pi
OMEGA_DOT_EARTH = 7.292115 * 10e-5


def get_ephemeris(navigation_entry, satellite_name):
    satellite = navigation_entry.sel(sv=satellite_name)
    ephem_list = []

    M0_array = satellite.get('M0').data
    deltaN_array = satellite.get('DeltaN').data
    sqrtA_array = satellite.get('sqrtA').data
    time_array = satellite.get('time').data
    e_array = satellite.get('Eccentricity').data
    omega_array = satellite.get('omega').data
    Crs_array = satellite.get('Crs').data
    Crc_array = satellite.get('Crc').data
    Cuc_array = satellite.get('Cuc').data
    Cus_array = satellite.get('Cus').data
    Cic_array = satellite.get('Cic').data
    Cis_array = satellite.get('Cis').data
    IDOT_array = satellite.get('IDOT').data
    Io_array = satellite.get('Io').data
    omega0_array = satellite.get('Omega0').data
    omegaDot_array = satellite.get('OmegaDot').data
    Toe_array = satellite.get('Toe').data

    for index, time in enumerate(time_array):
        M0 = M0_array[index]
        e = e_array[index]
        omega = omega_array[index]
        Crs = Crs_array[index]
        Crc = Crc_array[index]
        Cuc = Cuc_array[index]
        Cus = Cus_array[index]
        Cic = Cic_array[index]
        Cis = Cis_array[index]
        IDOT = IDOT_array[index]
        Io = Io_array[index]
        omega0 = omega0_array[index]
        omegaDot = omegaDot_array[index]
        deltaN = deltaN_array[index]
        sqrtA = sqrtA_array[index]
        Toe = Toe_array[index]
        if not np.isnan(np.array([e, M0, deltaN, sqrtA, omega, Crs, Crc, Cuc, Cus,
                                  Cic, Cis, IDOT, Io, omega0, omegaDot, Toe])).any():
            ephemeris = Ephemeris(time, e, M0, deltaN, sqrtA, omega, Crs, Crc, Cuc, Cus, Cic, Cis, IDOT, Io, omega0,
                                  omegaDot, Toe)
            ephem_list.append(ephemeris)
    return ephem_list


def get_ephemeris_from_nearest_reference_time(ephemeris_list, epoch):
    nearest_time = ephemeris_list[0].time
    nearest_index = 0
    for index, entry in enumerate(ephemeris_list):
        ref_time = entry.time
        if abs(epoch - ref_time) < abs(epoch - nearest_time) and epoch - ref_time >= np.timedelta64(seconds=0):
            nearest_time = ref_time
            nearest_index = index
    return ephemeris_list[nearest_index]


def find_tk_in_sec(epoch_ref, Toc):
    result = Toc - epoch_ref
    result_sec = result.astype('timedelta64[s]')
    if result_sec > np.timedelta64(seconds=302400):
        result_sec -= np.timedelta64(seconds=604800)
    if result_sec < np.timedelta64(seconds=-302400):
        result_sec += np.timedelta64(seconds=604800)
    return result_sec / np.timedelta64(1, 's')


def find_eccentric_anomaly(e, Mk):
    func = lambda Ek: Ek - e * math.sin(Ek) - Mk
    Ek_solution = fsolve(func, Mk)
    return Ek_solution[0]


def find_vk(e, Ek):
    vk = math.atan((math.sqrt(1 - math.pow(e, 2)) * math.sin(Ek)) / (math.cos(Ek) - e))
    return vk


def find_satellite_xyz(epoch, navigation_entry):
    ephemeris_list_g13 = get_ephemeris(navigation_entry, 'G13')
    ephemeris = get_ephemeris_from_nearest_reference_time(ephemeris_list_g13, epoch)

    tk = find_tk_in_sec(ephemeris.time, epoch)
    A = math.pow(ephemeris.sqrtA, 2)
    n0 = math.sqrt(GM / math.pow(A, 3))
    n = n0 + ephemeris.deltaN
    Mk = ephemeris.M0 + n * tk
    Ek = find_eccentric_anomaly(ephemeris.e, Mk)
    vk = find_vk(ephemeris.e, Ek)
    Fk = vk + ephemeris.omega
    delta_uk = ephemeris.Cus * math.sin(2 * Fk) + ephemeris.Cuc * math.cos(2 * Fk)
    delta_rk = ephemeris.Crs * math.sin(2 * Fk) + ephemeris.Crc * math.cos(2 * Fk)
    delta_ik = ephemeris.Cis * math.sin(2 * Fk) + ephemeris.Cic * math.cos(2 * Fk)
    uk = Fk + delta_uk
    rk = A * (1 - ephemeris.e * math.cos(Ek)) + delta_rk
    ik = ephemeris.Io + delta_ik + ephemeris.IDOT * tk
    xk_orbital_plane = rk * math.cos(uk)
    yk_orbital_plane = rk * math.sin(uk)
    omegak = ephemeris.omega0 + (ephemeris.omegaDot - OMEGA_DOT_EARTH) * tk - OMEGA_DOT_EARTH * ephemeris.Toe
    x = xk_orbital_plane + math.cos(omegak) - yk_orbital_plane * math.cos(ik) * math.sin(omegak)
    y = xk_orbital_plane + math.sin(omegak) + yk_orbital_plane * math.cos(ik) * math.cos(omegak)
    z = yk_orbital_plane * math.sin(ik)
    return x, y, z


def find_distance_3d(x1, y1, z1, x2, y2, z2):
    return math.sqrt(math.pow(x2 - x1, 2) +
                     math.pow(y2 - y1, 2) +
                     math.pow(z2 - z1, 2))



navigation_entry = gr.load('gps_data.n')
observation_entry = gr.load('obs_data.o', use='G', interval=10)
receiver_position = gr.rinexheader('obs_data.o').get('position')
times = observation_entry.get('time')
distance_array = []
time_array = []
for index, epoch in enumerate(times.data):
    xyz_satellite = find_satellite_xyz(epoch, navigation_entry)
    distance = find_distance_3d(xyz_satellite[0], xyz_satellite[1], xyz_satellite[2], receiver_position[0],
                                receiver_position[1], receiver_position[2])
    distance_array.append(distance)
    time_array.append(epoch)

plt.plot(time_array, distance_array)
plt.show()
