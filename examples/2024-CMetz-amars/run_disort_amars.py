#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

import canoe
from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.harp import radiation_band
from canoe.athena import ParameterInput, Mesh
from atm_profile_utils import read_atm_profile
from numpy import linspace, ones, exp, ndarray
from netCDF4 import Dataset
from pylab import *

def plot_opacity():
    config = load_configure("amarsw-rt.yaml")
    band = radiation_band("B1", config, load_opacity=True)
    aH2O = band.get_absorber_by_name("H2O")
    
    atm = read_atm_profile("amarsw.atm")
    ilyr = 0
    atm = [atm['TEM'][ilyr], atm['H2O'][ilyr], 0., 0., 0., atm['PRE'][ilyr]]
    print(atm)
    print(aH2O.get_attenuation(50.0, 50.0, atm))

if __name__ == "__main__":
    mesh, inp = canoe.start_with_input("amarsw.inp")

    mb = mesh.meshblock(0)
    rad = mb.get_rad()

    #plot_opacity()

    atm = read_atm_profile("amarsw.atm")
    mb.modify_atm(atm)

    rad.cal_flux(mb)
    #print(rad.fluxup)


    # aH2S = band.get_absorber_by_name("H2S")

    # print(aH2S.get_attenuation(50.0, 50.0, atm))

    # num_layers = 128
    # band.resize(num_layers)

    # atm = create_atmosphere(num_layers)
    # band.set_spectral_properties(atm)

    # dtau = band.get_dtau()
    # print(dtau.shape)
    # print(dtau)

    # band.cal_flux()
