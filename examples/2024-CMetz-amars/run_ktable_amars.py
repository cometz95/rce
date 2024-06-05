#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure, find_resource
from canoe.harp import radiation_band
from atm_profile_utils import read_atm_profile, get_species_names
from netCDF4 import Dataset
from numpy import *
from rfmlib import *

if __name__ == "__main__":
    hitran_file = find_resource("HITRAN2020.par")
    species = ['H2O']
    wghts = [18.]
    atm = read_atm_profile("amarsw-main.nc", species, wghts)

    config = load_configure("amarsw.yaml")
    def_species(vapors=species)

    for i in range(8):
        band = radiation_band(str(config["bands"][i]), config)

        nspec = band.get_num_specgrids()
        wmin, wmax = band.get_range()
        wres = (wmax - wmin) / (nspec - 1)

        # atmospheric properties
        num_layers = 128
        wav_grid = (wmin, wmax, wres)
        tem_grid = (5, -50, 50)

        driver = create_rfm_driver(wav_grid, tem_grid, species, hitran_file)

        # write rfm atmosphere file to file
        write_rfm_atm(atm)
        write_rfm_drv(driver)

        # run rfm and write kcoeff file
        run_rfm()

        fname = band.get_absorber_by_name("H2O").get_opacity_file()
        base_name = os.path.basename(fname)
        fname, _ = os.path.splitext(base_name)
        write_ktable(fname + "_B" + str(i + 1), species, atm, wav_grid, tem_grid)
