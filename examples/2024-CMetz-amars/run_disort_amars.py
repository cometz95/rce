#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

import canoe
from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.harp import radiation_band
from canoe.athena import ParameterInput, Mesh
from numpy import linspace, ones, exp, ndarray
from netCDF4 import Dataset
from pylab import *

if __name__ == "__main__":
    canoe.start()

    pin = ParameterInput()
    pin.load_from_file("amarsw.inp")

    vapors = pin.get_string("species", "vapor").split(", ")
    clouds = pin.get_string("species", "cloud").split(", ")
    # tracers = pin.get_string("species", "tracer").split(", ")

    def_species(vapors=vapors, clouds=clouds)
    def_thermo(pin)

    config = load_configure("amars.yaml")

    mesh = Mesh(pin)
    mesh.initialize(pin)

    mb = mesh.meshblock(0)
    rad = mb.get_rad()

    rad.cal_flux(mb)
    #print(rad.fluxup)

    # band = radiation_band("B1", config, load_opacity=True)

    # aH2O = band.get_absorber_by_name("H2O")
    # aH2S = band.get_absorber_by_name("H2S")

    # atm = [temp, xH2O, xH2S, xSO2, xCO2, v1, v2, v3, pres]
    # print(aH2O.get_attenuation(50.0, 50.0, atm))
    # print(aH2S.get_attenuation(50.0, 50.0, atm))

    # num_layers = 128
    # band.resize(num_layers)

    # atm = create_atmosphere(num_layers)
    # band.set_spectral_properties(atm)

    # dtau = band.get_dtau()
    # print(dtau.shape)
    # print(dtau)

    # band.cal_flux()
