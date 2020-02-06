############################
# Date: 24/12/2019
# Title: Script for plotting the positions of IRDIS sources
# Description: Use this script make a plot comparing the different values of position and errors from different sources
# VIP version: 0.9.11 (Rainot edit.)
# Python version: 3 ONLY
############################
# Some definitions

## Load libraries
import __init__
import sys
import matplotlib
import vip_hci
from hciplot import plot_frames, plot_cubes
from vip_hci.metrics.fakecomp import cube_inject_companions, frame_inject_companion, normalize_psf
import numpy as np
import scipy
import astropy.io.fits as fits
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as mh
import csv

# Find the files containing the values
## MCMC
mcmc_pos_K1_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK1_MCMC_pos.txt'
mcmc_PA_K1_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK1_MCMC_PA.txt'
mcmc_flux_K1_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK1_MCMC_flux.txt'
mcmc_pos_K2_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK2_MCMC_pos.txt'
mcmc_PA_K2_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK2_MCMC_PA.txt'
mcmc_flux_K2_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/MCMC/IRDISK2_MCMC_flux.txt'

## Simplex
simplex_K1_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K1.txt'
simplex_K2_file = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt'

## Simplex errors
simplex_err_K1_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/Simplex_errors/'
simplex_err_K2_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/Simplex_errors/'

## Julia PSF photometry
julia_K1_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/Julias_PSF_fitting/K1/psf_photometry_K1.csv'
julia_K2_path = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/Julias_PSF_fitting/K2/psf_photometry_K2.csv'

# Central star aperture photometry flux
#central_star_flux = np.array([1793119.1,1033893.75])
central_star_flux = np.array([2024466.5281759,1138772.21601975])

# Center coordinates of the image
centy = 512
centx = 512

# Plots
pos_K1 = False
pos_K2 = False
PA_K1 = False
PA_K2 = False
flux_K1 = False
flux_K2 = False

pos_K1_err = False
pos_K2_err = False
PA_K1_err = False
PA_K2_err = False
flux_K1_err = False
flux_K2_err = False
dmag_simplex_K1_mcmc = False
dmag_simplex_K2_mcmc = False

x_simplex_K1_julia = False
y_simplex_K1_julia = False
x_simplex_K2_julia = False
y_simplex_K2_julia = False
dmag_simplex_K1_julia = False
dmag_simplex_K2_julia = False

x_mcmc_K1_julia = False
y_mcmc_K1_julia = False
x_mcmc_K2_julia = False
y_mcmc_K2_julia = False
dmag_mcmc_K1_julia = False
dmag_mcmc_K2_julia = False
dmag_mcmc_K1_K2 = True
dmag_julia_K1_K2 = True

# Define arrays
mcmc_pos_K1 = np.zeros((15,2))
mcmc_x_K1 = np.zeros((15,2))
mcmc_y_K1 = np.zeros((15,2))
mcmc_PA_K1 = np.zeros_like(mcmc_pos_K1)
mcmc_flux_K1 = np.zeros_like(mcmc_pos_K1)
mcmc_dmag_K1 = np.zeros_like(mcmc_pos_K1)
mcmc_pos_K2 = np.zeros((15,2))
mcmc_x_K2 = np.zeros((15,2))
mcmc_y_K2 = np.zeros((15,2))
mcmc_PA_K2 = np.zeros_like(mcmc_pos_K2)
mcmc_flux_K2 = np.zeros_like(mcmc_pos_K2)
mcmc_dmag_K2 = np.zeros_like(mcmc_pos_K2)

simplex_pos_K1 = np.zeros((19,2))
simplex_x_K1 = np.zeros((19,2))
simplex_y_K1 = np.zeros((19,2))
simplex_PA_K1 = np.zeros_like(simplex_pos_K1)
simplex_flux_K1 = np.zeros_like(simplex_pos_K1)
simplex_dmag_K1 = np.zeros_like(simplex_pos_K1)
simplex_pos_K2 = np.zeros((19,2))
simplex_x_K2 = np.zeros((19,2))
simplex_y_K2 = np.zeros((19,2))
simplex_PA_K2 = np.zeros_like(simplex_pos_K2)
simplex_flux_K2 = np.zeros_like(simplex_pos_K2)
simplex_dmag_K2 = np.zeros_like(simplex_pos_K2)

## Logarithmic arrays
mcmc_pos_K1_log = np.zeros((15,2))
mcmc_PA_K1_log = np.zeros_like(mcmc_pos_K1_log)
mcmc_flux_K1_log = np.zeros_like(mcmc_pos_K1_log)
mcmc_pos_K2_log = np.zeros((15,2))
mcmc_PA_K2_log = np.zeros_like(mcmc_pos_K2_log)
mcmc_flux_K2_log = np.zeros_like(mcmc_pos_K2_log)

simplex_pos_K1_log = np.zeros((19,2))
simplex_PA_K1_log = np.zeros_like(simplex_pos_K1_log)
simplex_flux_K1_log = np.zeros_like(simplex_pos_K1_log)
simplex_pos_K2_log = np.zeros((19,2))
simplex_PA_K2_log = np.zeros_like(simplex_pos_K2_log)
simplex_flux_K2_log = np.zeros_like(simplex_pos_K2_log)

## Julia's photometry
julia_x_K1 = np.zeros((16,2))
julia_y_K1 = np.zeros_like(julia_x_K1)
julia_pos_K1 = np.zeros_like(julia_x_K1)
julia_pos_K1_arcsec = np.zeros_like(julia_pos_K1)
julia_PA_K1 = np.zeros_like(julia_pos_K1)
julia_flux_K1 = np.zeros_like(julia_x_K1)
julia_dmag_K1 = np.zeros_like(julia_x_K1)
julia_x_K2 = np.zeros((16,2))
julia_y_K2 = np.zeros_like(julia_x_K2)
julia_pos_K2 = np.zeros_like(julia_x_K2)
julia_pos_K2_arcsec = np.zeros_like(julia_x_K2)
julia_PA_K2 = np.zeros_like(julia_pos_K2)
julia_flux_K2 = np.zeros_like(julia_x_K2)
julia_dmag_K2 = np.zeros_like(julia_x_K2)

# Open the files and save the values to previously defined arrays
j = 0 # Create count variable
with open(mcmc_pos_K1_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_pos_K1[j,0] = float(columns[0])
        mcmc_pos_K1[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(mcmc_PA_K1_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_PA_K1[j,0] = float(columns[0])
        mcmc_PA_K1[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(mcmc_flux_K1_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_flux_K1[j,0] = float(columns[0])
        mcmc_flux_K1[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(mcmc_pos_K2_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_pos_K2[j,0] = float(columns[0])
        mcmc_pos_K2[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(mcmc_PA_K2_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_PA_K2[j,0] = float(columns[0])
        mcmc_PA_K2[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(mcmc_flux_K2_file) as f:
    for line in f:
        line = line.strip()
        columns = line.split()
        mcmc_flux_K2[j,0] = float(columns[0])
        mcmc_flux_K2[j,1] = float(columns[1])
        j+=1
f.close

j = 0
with open(simplex_K1_file) as f:
    header1 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        simplex_pos_K1[j,0] = float(columns[0])
        simplex_PA_K1[j,0] = float(columns[1])
        simplex_flux_K1[j,0] = float(columns[2])
        j+=1
f.close

j = 0
with open(simplex_K2_file) as f:
    header1 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        simplex_pos_K2[j,0] = float(columns[0])
        simplex_PA_K2[j,0] = float(columns[1])
        simplex_flux_K2[j,0] = float(columns[2])
        j+=1
f.close
simplex_flux_K1[15,0] = 1.923620270692435241e+01
# Simplex errors
for i in range(0,17):
    with open(simplex_err_K1_path + 'Simplex_K1_S{}.txt'.format(i),'rb') as f:
        header1 = f.readline()
        header2 = f.readline()
        value = np.zeros(3)
        j=0
        for line in f:
            line = line.strip()
            value[j] = line
            j+=1
        simplex_pos_K1[i,1] = value[0]
        simplex_PA_K1[i,1] = value[1]
        simplex_flux_K1[i,1] = value[2]
    f.close
    with open(simplex_err_K2_path + 'Simplex_K2_S{}.txt'.format(i),'rb') as f:
        header1 = f.readline()
        header2 = f.readline()
        value = np.zeros(3)
        j=0
        for line in f:
            line = line.strip()
            value[j] = line
            j+=1
        simplex_pos_K2[i,1] = value[0]
        simplex_PA_K2[i,1] = value[1]
        simplex_flux_K2[i,1] = value[2]
    f.close
    j = 0
simplex_flux_K1[:,1] = [13.99371532,  9.53036689,  2.15003922,  2.4995157 ,  3.81966204,
        1.18707697,  0.86711963,  1.80330878,  1.53391872,  2.67004237,
        2.4007351 ,  2.48527375,  2.20384686,  2.50508398,  3.26460592,
        0.58607642,  2.5617156,0,0]
simplex_flux_K2[:,1] = [8.32191222, 7.71549515, 2.08340082, 2.61796825, 1.28565406,
       2.97004946, 2.6005361 , 3.32050519, 1.84634837, 1.36315893,
       0.90722712, 1.06433667, 1.6302751 , 2.78925085, 1.0655311 ,
       2.94670661, 2.12218066,0,0]
mcmc_dmag_K1[1,0] = 3.7
mcmc_dmag_K2[1,0] = 3.8

# Convert r, PA to X, Y coordinates
for i in range(len(mcmc_pos_K1)):
    mcmc_x_K1[i,0] = mcmc_pos_K1[i,0] * np.cos(np.deg2rad(mcmc_PA_K1[i,0])) + centx
    mcmc_x_K1[i,1] = np.sqrt((mh.cos(np.deg2rad(mcmc_PA_K1[i,0]))*mcmc_pos_K1[i,1])**2 + (mcmc_pos_K1[i,0]*mh.sin(np.deg2rad(mcmc_PA_K1[i,0]))*np.deg2rad(mcmc_PA_K1[i,1]))**2)
    mcmc_y_K1[i,0] = mcmc_pos_K1[i,0] * np.sin(np.deg2rad(mcmc_PA_K1[i,0])) + centy
    mcmc_y_K1[i,1] = np.sqrt((mh.sin(np.deg2rad(mcmc_PA_K1[i,0]))*mcmc_pos_K1[i,1])**2 + (mcmc_pos_K1[i,0]*mh.cos(np.deg2rad(mcmc_PA_K1[i,0]))*np.deg2rad(mcmc_PA_K1[i,1]))**2)
    mcmc_x_K2[i,0] = mcmc_pos_K2[i,0] * np.cos(np.deg2rad(mcmc_PA_K2[i,0])) + centx
    mcmc_x_K2[i,1] = np.sqrt((mh.cos(np.deg2rad(mcmc_PA_K2[i,0]))*mcmc_pos_K2[i,1])**2 + (mcmc_pos_K2[i,0]*mh.sin(np.deg2rad(mcmc_PA_K2[i,0]))*np.deg2rad(mcmc_PA_K2[i,1]))**2)
    mcmc_y_K2[i,0] = mcmc_pos_K2[i,0] * np.sin(np.deg2rad(mcmc_PA_K2[i,0])) + centy
    mcmc_y_K2[i,1] = np.sqrt((mh.sin(np.deg2rad(mcmc_PA_K2[i,0]))*mcmc_pos_K2[i,1])**2 + (mcmc_pos_K2[i,0]*mh.cos(np.deg2rad(mcmc_PA_K2[i,0]))*np.deg2rad(mcmc_PA_K2[i,1]))**2)

for i in range(len(simplex_pos_K1)):
    simplex_x_K1[i,0] = simplex_pos_K1[i,0] * np.cos(np.deg2rad(simplex_PA_K1[i,0])) + centx
    simplex_x_K1[i,1] = np.sqrt((mh.cos(np.deg2rad(simplex_PA_K1[i,0]))*simplex_pos_K1[i,1])**2 + (simplex_pos_K1[i,0]*mh.sin(np.deg2rad(simplex_PA_K1[i,0]))*np.deg2rad(simplex_PA_K1[i,1]))**2)
    simplex_y_K1[i,0] = simplex_pos_K1[i,0] * np.sin(np.deg2rad(simplex_PA_K1[i,0])) + centy
    simplex_y_K1[i,1] = np.sqrt((mh.sin(np.deg2rad(simplex_PA_K1[i,0]))*simplex_pos_K1[i,1])**2 + (simplex_pos_K1[i,0]*mh.cos(np.deg2rad(simplex_PA_K1[i,0]))*np.deg2rad(simplex_PA_K1[i,1]))**2)
    #simplex_y_K1[i,1] = np.sqrt((simplex_pos_K1[i,1]/simplex_pos_K1[i,0])**2 + (np.deg2rad(simplex_PA_K1[i,1])/mh.tan(np.deg2rad(simplex_PA_K1[i,0])))**2)
    simplex_x_K2[i,0] = simplex_pos_K2[i,0] * np.cos(np.deg2rad(simplex_PA_K2[i,0])) + centx
    simplex_x_K2[i,1] = np.sqrt((mh.cos(np.deg2rad(simplex_PA_K2[i,0]))*simplex_pos_K2[i,1])**2 + (simplex_pos_K2[i,0]*mh.sin(np.deg2rad(simplex_PA_K2[i,0]))*np.deg2rad(simplex_PA_K2[i,1]))**2)
    simplex_y_K2[i,0] = simplex_pos_K2[i,0] * np.sin(np.deg2rad(simplex_PA_K2[i,0])) + centy
    simplex_y_K2[i,1] = np.sqrt((mh.sin(np.deg2rad(simplex_PA_K2[i,0]))*simplex_pos_K2[i,1])**2 + (simplex_pos_K2[i,0]*mh.cos(np.deg2rad(simplex_PA_K2[i,0]))*np.deg2rad(simplex_PA_K2[i,1]))**2)

# Get magnitude contrasts
for i in range(len(mcmc_flux_K1)):
    if i==1:
        continue
    mcmc_dmag_K1[i,0] = -2.5*mh.log10(mcmc_flux_K1[i,0]/central_star_flux[0])
    # mcmc_dmag_K1[i,1] = 2.5*mcmc_flux_K1[i,1]/(mcmc_flux_K1[i,0]*mh.log(10))
    mcmc_dmag_K1[i,1] = np.sqrt((2.5*mcmc_flux_K1[i,1]/(mcmc_flux_K1[i,0]*np.log(10)))**2 + (2.5*21453.893/(2021121.4*np.log(10)))**2)
    mcmc_dmag_K2[i,0] = -2.5*mh.log10(mcmc_flux_K2[i,0]/central_star_flux[1])
    # mcmc_dmag_K2[i,1] = 2.5*mcmc_flux_K2[i,1]/(mcmc_flux_K2[i,0]*mh.log(10))
    mcmc_dmag_K2[i,1] = np.sqrt((2.5*mcmc_flux_K2[i,1]/(mcmc_flux_K2[i,0]*np.log(10)))**2 + (2.5*10503.037/(1149148.6*np.log(10)))**2)

for i in range(len(simplex_flux_K1)-2):
    simplex_dmag_K1[i,0] = -2.5*mh.log10(simplex_flux_K1[i,0]/central_star_flux[0])
    # simplex_dmag_K1[i,1] = 2.5*simplex_flux_K1[i,1]/(simplex_flux_K1[i,0]*mh.log(10))
    simplex_dmag_K1[i,1] = np.sqrt((2.5*simplex_flux_K1[i,1]/(simplex_flux_K1[i,0]*np.log(10)))**2 + (2.5*21453.893/(2021121.4*np.log(10)))**2)
    simplex_dmag_K2[i,0] = -2.5*mh.log10(simplex_flux_K2[i,0]/central_star_flux[1])
    # simplex_dmag_K2[i,1] = 2.5*simplex_flux_K2[i,1]/(simplex_flux_K2[i,0]*mh.log(10))
    simplex_dmag_K2[i,1] = np.sqrt((2.5*simplex_flux_K2[i,1]/(simplex_flux_K2[i,0]*np.log(10)))**2 + (2.5*10503.037/(1149148.6*np.log(10)))**2)

# Log
for i in range(len(mcmc_pos_K1)):
    if i==1:
        continue
    mcmc_pos_K1_log[i,0] = mh.log10(mcmc_pos_K1[i,0])
    mcmc_pos_K1_log[i,1] = mcmc_pos_K1[i,1]/(mcmc_pos_K1[i,0]*mh.log(10))
    mcmc_PA_K1_log[i,0] = mh.log10(mcmc_PA_K1[i,0])
    mcmc_PA_K1_log[i,1] = mcmc_PA_K1[i,1]/(mcmc_PA_K1[i,0]*mh.log(10))
    mcmc_flux_K1_log[i,0] = mh.log10(mcmc_flux_K1[i,0])
    mcmc_flux_K1_log[i,1] = mcmc_flux_K1[i,1]/(mcmc_flux_K1[i,0]*mh.log(10))
    mcmc_pos_K2_log[i,0] = mh.log10(mcmc_pos_K2[i,0])
    mcmc_pos_K2_log[i,1] = mcmc_pos_K2[i,1]/(mcmc_pos_K2[i,0]*mh.log(10))
    mcmc_PA_K2_log[i,0] = mh.log10(mcmc_PA_K2[i,0])
    mcmc_PA_K2_log[i,1] = mcmc_PA_K2[i,1]/(mcmc_PA_K2[i,0]*mh.log(10))
    mcmc_flux_K2_log[i,0] = mh.log10(mcmc_flux_K2[i,0])
    mcmc_flux_K2_log[i,1] = mcmc_flux_K2[i,1]/(mcmc_flux_K2[i,0]*mh.log(10))

for i in range(len(simplex_pos_K1)-2):
    simplex_pos_K1_log[i,0] = mh.log10(simplex_pos_K1[i,0])
    simplex_pos_K1_log[i,1] = simplex_pos_K1[i,1]/(simplex_pos_K1[i,0]*mh.log(10))
    simplex_PA_K1_log[i,0] = mh.log10(simplex_PA_K1[i,0])
    simplex_PA_K1_log[i,1] = simplex_PA_K1[i,1]/(simplex_PA_K1[i,0]*mh.log(10))
    simplex_flux_K1_log[i,0] = mh.log10(simplex_flux_K1[i,0])
    simplex_flux_K1_log[i,1] = simplex_flux_K1[i,1]/(simplex_flux_K1[i,0]*mh.log(10))
    simplex_pos_K2_log[i,0] = mh.log10(simplex_pos_K2[i,0])
    simplex_pos_K2_log[i,1] = simplex_pos_K2[i,1]/(simplex_pos_K2[i,0]*mh.log(10))
    simplex_PA_K2_log[i,0] = mh.log10(simplex_PA_K2[i,0])
    simplex_PA_K2_log[i,1] = simplex_PA_K2[i,1]/(simplex_PA_K2[i,0]*mh.log(10))
    simplex_flux_K2_log[i,0] = mh.log10(simplex_flux_K2[i,0])
    simplex_flux_K2_log[i,1] = simplex_flux_K2[i,1]/(simplex_flux_K2[i,0]*mh.log(10))

# Julia PSF photometry
## K1
with open(julia_K1_path) as csvfile:
   reader = csv.reader(csvfile, delimiter=';')
   next(reader, None) #Skip header
   next(reader, None) #Skip header
   j = 0
   for row in reader:
       julia_x_K1[j,0] = float(row[1])
       julia_x_K1[j,1] = float(row[2])
       julia_y_K1[j,0] = float(row[3])
       julia_y_K1[j,1] = float(row[4])
       julia_pos_K1[j,0] = np.sqrt((julia_x_K1[j,0]-centx)**2+(julia_y_K1[j,0]-centy)**2)
       julia_pos_K1[j,1] = np.sqrt((julia_x_K1[j,0]/julia_pos_K1[j,0])**2*julia_x_K1[j,1]**2 + (julia_y_K1[j,0]/julia_pos_K1[j,0])**2*julia_y_K1[j,1]**2)
       julia_flux_K1[j,0] = float(row[6])
       julia_flux_K1[j,1] = float(row[7])#+(4671.604084//float(row[5]))**2)
       # julia_flux_K1[j,0] = float(row[5])/48508723.081159
       # julia_flux_K1[j,1] = np.sqrt((float(row[6])/48508723.081159)**2)#+(4671.604084//float(row[5]))**2)
       julia_dmag_K1[j,0] = -2.5*mh.log10(julia_flux_K1[j,0]/3748275.952813)
       julia_dmag_K1[j,1] = 2.5*julia_flux_K1[j,1]/(julia_flux_K1[j,0]*mh.log(10))
       julia_dmag_K1[j,1] = np.sqrt((2.5*julia_flux_K1[j,1]/(julia_flux_K1[j,0]*np.log(10)))**2 + (2.5*177.215684/(3748275.952813*np.log(10)))**2)
       j+=1

## K2
with open(julia_K2_path) as csvfile:
   reader = csv.reader(csvfile, delimiter=';')
   next(reader, None) #Skip header
   next(reader, None) #Skip header
   j = 0
   for row in reader:
       julia_x_K2[j,0] = float(row[1])
       julia_x_K2[j,1] = float(row[2])
       julia_y_K2[j,0] = float(row[3])
       julia_y_K2[j,1] = float(row[4])
       julia_pos_K2[j,0] = np.sqrt((julia_x_K2[j,0]-centx)**2+(julia_y_K2[j,0]-centy)**2)
       julia_pos_K2[j,1] = np.sqrt((julia_x_K2[j,0]/julia_pos_K2[j,0])**2*julia_x_K2[j,1]**2 + (julia_y_K2[j,0]/julia_pos_K2[j,0])**2*julia_y_K2[j,1]**2)
       julia_flux_K2[j,0] = np.abs(float(row[6]))
       julia_flux_K2[j,1] = float(row[7])
       # julia_flux_K2[j,0] = float(row[5])/27874094.331858
       # julia_flux_K2[j,1] = julia_flux_K2[j,0] * np.sqrt((float(row[6])/float(row[5]))**2+(3266.130664/27874094.331858)**2)
       julia_dmag_K2[j,0] = -2.5*mh.log10(julia_flux_K2[j,0]/2158544.286605)
       #julia_dmag_K2[j,1] = 2.5*julia_flux_K2[j,1]/(julia_flux_K2[j,0]*mh.log(10))
       julia_dmag_K2[j,1] = np.sqrt((2.5*julia_flux_K2[j,1]/(julia_flux_K2[j,0]*np.log(10)))**2 + (2.5*121.569173/(2158544.286605*np.log(10)))**2)
       j+=1


# Plotting
colors = cm.rainbow(np.linspace(0, 1, len(simplex_pos_K1)))

## Position K1
if pos_K1 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("r - K1")
    plt.errorbar(mcmc_pos_K1[0,0], simplex_pos_K1[0,0],xerr=mcmc_pos_K1[0,1],yerr=simplex_pos_K1[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_pos_K1[1,0], simplex_pos_K1[1,0],xerr=mcmc_pos_K1[1,1],yerr=simplex_pos_K1[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_pos_K1[2,0], simplex_pos_K1[2,0],xerr=mcmc_pos_K1[2,1],yerr=simplex_pos_K1[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_pos_K1[3,0], simplex_pos_K1[3,0],xerr=mcmc_pos_K1[3,1],yerr=simplex_pos_K1[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_pos_K1[4,0], simplex_pos_K1[4,0],xerr=mcmc_pos_K1[4,1],yerr=simplex_pos_K1[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_pos_K1[5,0], simplex_pos_K1[5,0],xerr=mcmc_pos_K1[5,1],yerr=simplex_pos_K1[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_pos_K1[6,0], simplex_pos_K1[6,0],xerr=mcmc_pos_K1[6,1],yerr=simplex_pos_K1[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_pos_K1[7,0], simplex_pos_K1[7,0],xerr=mcmc_pos_K1[7,1],yerr=simplex_pos_K1[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_pos_K1[8,0], simplex_pos_K1[8,0],xerr=mcmc_pos_K1[8,1],yerr=simplex_pos_K1[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_pos_K1[9,0], simplex_pos_K1[9,0],xerr=mcmc_pos_K1[9,1],yerr=simplex_pos_K1[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_pos_K1[10,0], simplex_pos_K1[10,0],xerr=mcmc_pos_K1[10,1],yerr=simplex_pos_K1[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_pos_K1[11,0], simplex_pos_K1[11,0],xerr=mcmc_pos_K1[11,1],yerr=simplex_pos_K1[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_pos_K1[12,0], simplex_pos_K1[12,0],xerr=mcmc_pos_K1[12,1],yerr=simplex_pos_K1[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_pos_K1[13,0], simplex_pos_K1[13,0],xerr=mcmc_pos_K1[13,1],yerr=simplex_pos_K1[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_pos_K1[14,0], simplex_pos_K1[14,0],xerr=mcmc_pos_K1[14,1],yerr=simplex_pos_K1[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_pos_K1[15,0],yerr=simplex_pos_K1[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_pos_K1[16,0],yerr=simplex_pos_K1[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_pos_K1[17,0],yerr=simplex_pos_K1[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_pos_K1[18,0],yerr=simplex_pos_K1[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,3,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Position K2
if pos_K2 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("r - K2")
    plt.errorbar(mcmc_pos_K2[0,0], simplex_pos_K2[0,0],xerr=mcmc_pos_K2[0,1],yerr=simplex_pos_K2[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_pos_K2[1,0], simplex_pos_K2[1,0],xerr=mcmc_pos_K2[1,1],yerr=simplex_pos_K2[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_pos_K2[2,0], simplex_pos_K2[2,0],xerr=mcmc_pos_K2[2,1],yerr=simplex_pos_K2[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_pos_K2[3,0], simplex_pos_K2[3,0],xerr=mcmc_pos_K2[3,1],yerr=simplex_pos_K2[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_pos_K2[4,0], simplex_pos_K2[4,0],xerr=mcmc_pos_K2[4,1],yerr=simplex_pos_K2[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_pos_K2[5,0], simplex_pos_K2[5,0],xerr=mcmc_pos_K2[5,1],yerr=simplex_pos_K2[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_pos_K2[6,0], simplex_pos_K2[6,0],xerr=mcmc_pos_K2[6,1],yerr=simplex_pos_K2[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_pos_K2[7,0], simplex_pos_K2[7,0],xerr=mcmc_pos_K2[7,1],yerr=simplex_pos_K2[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_pos_K2[8,0], simplex_pos_K2[8,0],xerr=mcmc_pos_K2[8,1],yerr=simplex_pos_K2[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_pos_K2[9,0], simplex_pos_K2[9,0],xerr=mcmc_pos_K2[9,1],yerr=simplex_pos_K2[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_pos_K2[10,0], simplex_pos_K2[10,0],xerr=mcmc_pos_K2[10,1],yerr=simplex_pos_K2[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_pos_K2[11,0], simplex_pos_K2[11,0],xerr=mcmc_pos_K2[11,1],yerr=simplex_pos_K2[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_pos_K2[12,0], simplex_pos_K2[12,0],xerr=mcmc_pos_K2[12,1],yerr=simplex_pos_K2[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_pos_K2[13,0], simplex_pos_K2[13,0],xerr=mcmc_pos_K2[13,1],yerr=simplex_pos_K2[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_pos_K2[14,0], simplex_pos_K2[14,0],xerr=mcmc_pos_K2[14,1],yerr=simplex_pos_K2[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_pos_K2[15,0],yerr=simplex_pos_K2[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_pos_K2[16,0],yerr=simplex_pos_K2[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_pos_K2[17,0],yerr=simplex_pos_K2[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_pos_K2[18,0],yerr=simplex_pos_K2[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,3,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Position-Angle K1
if PA_K1 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("PA - K1")
    plt.errorbar(mcmc_PA_K1[0,0], simplex_PA_K1[0,0],xerr=mcmc_PA_K1[0,1],yerr=simplex_PA_K1[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_PA_K1[1,0], simplex_PA_K1[1,0],xerr=mcmc_PA_K1[1,1],yerr=simplex_PA_K1[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_PA_K1[2,0], simplex_PA_K1[2,0],xerr=mcmc_PA_K1[2,1],yerr=simplex_PA_K1[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_PA_K1[3,0], simplex_PA_K1[3,0],xerr=mcmc_PA_K1[3,1],yerr=simplex_PA_K1[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_PA_K1[4,0], simplex_PA_K1[4,0],xerr=mcmc_PA_K1[4,1],yerr=simplex_PA_K1[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_PA_K1[5,0], simplex_PA_K1[5,0],xerr=mcmc_PA_K1[5,1],yerr=simplex_PA_K1[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_PA_K1[6,0], simplex_PA_K1[6,0],xerr=mcmc_PA_K1[6,1],yerr=simplex_PA_K1[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_PA_K1[7,0], simplex_PA_K1[7,0],xerr=mcmc_PA_K1[7,1],yerr=simplex_PA_K1[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_PA_K1[8,0], simplex_PA_K1[8,0],xerr=mcmc_PA_K1[8,1],yerr=simplex_PA_K1[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_PA_K1[9,0], simplex_PA_K1[9,0],xerr=mcmc_PA_K1[9,1],yerr=simplex_PA_K1[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_PA_K1[10,0], simplex_PA_K1[10,0],xerr=mcmc_PA_K1[10,1],yerr=simplex_PA_K1[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_PA_K1[11,0], simplex_PA_K1[11,0],xerr=mcmc_PA_K1[11,1],yerr=simplex_PA_K1[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_PA_K1[12,0], simplex_PA_K1[12,0],xerr=mcmc_PA_K1[12,1],yerr=simplex_PA_K1[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_PA_K1[13,0], simplex_PA_K1[13,0],xerr=mcmc_PA_K1[13,1],yerr=simplex_PA_K1[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_PA_K1[14,0], simplex_PA_K1[14,0],xerr=mcmc_PA_K1[14,1],yerr=simplex_PA_K1[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_PA_K1[15,0],yerr=simplex_PA_K1[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_PA_K1[16,0],yerr=simplex_PA_K1[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_PA_K1[17,0],yerr=simplex_PA_K1[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_PA_K1[18,0],yerr=simplex_PA_K1[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,10,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Position-Angle K2
if PA_K2 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("PA - K2")
    plt.errorbar(mcmc_PA_K2[0,0], simplex_PA_K2[0,0],xerr=mcmc_PA_K2[0,1],yerr=simplex_PA_K2[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_PA_K2[1,0], simplex_PA_K2[1,0],xerr=mcmc_PA_K2[1,1],yerr=simplex_PA_K2[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_PA_K2[2,0], simplex_PA_K2[2,0],xerr=mcmc_PA_K2[2,1],yerr=simplex_PA_K2[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_PA_K2[3,0], simplex_PA_K2[3,0],xerr=mcmc_PA_K2[3,1],yerr=simplex_PA_K2[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_PA_K2[4,0], simplex_PA_K2[4,0],xerr=mcmc_PA_K2[4,1],yerr=simplex_PA_K2[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_PA_K2[5,0], simplex_PA_K2[5,0],xerr=mcmc_PA_K2[5,1],yerr=simplex_PA_K2[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_PA_K2[6,0], simplex_PA_K2[6,0],xerr=mcmc_PA_K2[6,1],yerr=simplex_PA_K2[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_PA_K2[7,0], simplex_PA_K2[7,0],xerr=mcmc_PA_K2[7,1],yerr=simplex_PA_K2[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_PA_K2[8,0], simplex_PA_K2[8,0],xerr=mcmc_PA_K2[8,1],yerr=simplex_PA_K2[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_PA_K2[9,0], simplex_PA_K2[9,0],xerr=mcmc_PA_K2[9,1],yerr=simplex_PA_K2[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_PA_K2[10,0], simplex_PA_K2[10,0],xerr=mcmc_PA_K2[10,1],yerr=simplex_PA_K2[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_PA_K2[11,0], simplex_PA_K2[11,0],xerr=mcmc_PA_K2[11,1],yerr=simplex_PA_K2[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_PA_K2[12,0], simplex_PA_K2[12,0],xerr=mcmc_PA_K2[12,1],yerr=simplex_PA_K2[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_PA_K2[13,0], simplex_PA_K2[13,0],xerr=mcmc_PA_K2[13,1],yerr=simplex_PA_K2[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_PA_K2[14,0], simplex_PA_K2[14,0],xerr=mcmc_PA_K2[14,1],yerr=simplex_PA_K2[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_PA_K2[15,0],yerr=simplex_PA_K2[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_PA_K2[16,0],yerr=simplex_PA_K2[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_PA_K2[17,0],yerr=simplex_PA_K2[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_PA_K2[18,0],yerr=simplex_PA_K2[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,10,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Flux K1
if flux_K1 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("flux - K1")
    plt.errorbar(mcmc_flux_K1[0,0], simplex_flux_K1[0,0],xerr=mcmc_flux_K1[0,1],yerr=simplex_flux_K1[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_flux_K1[1,0], simplex_flux_K1[1,0],xerr=mcmc_flux_K1[1,1],yerr=simplex_flux_K1[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_flux_K1[2,0], simplex_flux_K1[2,0],xerr=mcmc_flux_K1[2,1],yerr=simplex_flux_K1[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_flux_K1[3,0], simplex_flux_K1[3,0],xerr=mcmc_flux_K1[3,1],yerr=simplex_flux_K1[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_flux_K1[4,0], simplex_flux_K1[4,0],xerr=mcmc_flux_K1[4,1],yerr=simplex_flux_K1[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_flux_K1[5,0], simplex_flux_K1[5,0],xerr=mcmc_flux_K1[5,1],yerr=simplex_flux_K1[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_flux_K1[6,0], simplex_flux_K1[6,0],xerr=mcmc_flux_K1[6,1],yerr=simplex_flux_K1[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_flux_K1[7,0], simplex_flux_K1[7,0],xerr=mcmc_flux_K1[7,1],yerr=simplex_flux_K1[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_flux_K1[8,0], simplex_flux_K1[8,0],xerr=mcmc_flux_K1[8,1],yerr=simplex_flux_K1[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_flux_K1[9,0], simplex_flux_K1[9,0],xerr=mcmc_flux_K1[9,1],yerr=simplex_flux_K1[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_flux_K1[10,0], simplex_flux_K1[10,0],xerr=mcmc_flux_K1[10,1],yerr=simplex_flux_K1[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_flux_K1[11,0], simplex_flux_K1[11,0],xerr=mcmc_flux_K1[11,1],yerr=simplex_flux_K1[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_flux_K1[12,0], simplex_flux_K1[12,0],xerr=mcmc_flux_K1[12,1],yerr=simplex_flux_K1[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_flux_K1[13,0], simplex_flux_K1[13,0],xerr=mcmc_flux_K1[13,1],yerr=simplex_flux_K1[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_flux_K1[14,0], simplex_flux_K1[14,0],xerr=mcmc_flux_K1[14,1],yerr=simplex_flux_K1[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_flux_K1[15,0],yerr=simplex_flux_K1[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_flux_K1[16,0],yerr=simplex_flux_K1[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_flux_K1[17,0],yerr=simplex_flux_K1[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_flux_K1[18,0],yerr=simplex_flux_K1[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,10,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Flux K2
if flux_K2 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Simplex", fontsize=20)
    plt.xlabel('MCMC', fontsize=20)
    plt.title("flux - K2")
    plt.errorbar(mcmc_flux_K2[0,0], simplex_flux_K2[0,0],xerr=mcmc_flux_K2[0,1],yerr=simplex_flux_K2[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_flux_K2[1,0], simplex_flux_K2[1,0],xerr=mcmc_flux_K2[1,1],yerr=simplex_flux_K2[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_flux_K2[2,0], simplex_flux_K2[2,0],xerr=mcmc_flux_K2[2,1],yerr=simplex_flux_K2[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_flux_K2[3,0], simplex_flux_K2[3,0],xerr=mcmc_flux_K2[3,1],yerr=simplex_flux_K2[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_flux_K2[4,0], simplex_flux_K2[4,0],xerr=mcmc_flux_K2[4,1],yerr=simplex_flux_K2[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_flux_K2[5,0], simplex_flux_K2[5,0],xerr=mcmc_flux_K2[5,1],yerr=simplex_flux_K2[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_flux_K2[6,0], simplex_flux_K2[6,0],xerr=mcmc_flux_K2[6,1],yerr=simplex_flux_K2[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_flux_K2[7,0], simplex_flux_K2[7,0],xerr=mcmc_flux_K2[7,1],yerr=simplex_flux_K2[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_flux_K2[8,0], simplex_flux_K2[8,0],xerr=mcmc_flux_K2[8,1],yerr=simplex_flux_K2[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_flux_K2[9,0], simplex_flux_K2[9,0],xerr=mcmc_flux_K2[9,1],yerr=simplex_flux_K2[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_flux_K2[10,0], simplex_flux_K2[10,0],xerr=mcmc_flux_K2[10,1],yerr=simplex_flux_K2[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_flux_K2[11,0], simplex_flux_K2[11,0],xerr=mcmc_flux_K2[11,1],yerr=simplex_flux_K2[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_flux_K2[12,0], simplex_flux_K2[12,0],xerr=mcmc_flux_K2[12,1],yerr=simplex_flux_K2[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_flux_K2[13,0], simplex_flux_K2[13,0],xerr=mcmc_flux_K2[13,1],yerr=simplex_flux_K2[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_flux_K2[14,0], simplex_flux_K2[14,0],xerr=mcmc_flux_K2[14,1],yerr=simplex_flux_K2[14,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(0, simplex_flux_K2[15,0],yerr=simplex_flux_K2[15,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(0, simplex_flux_K2[16,0],yerr=simplex_flux_K2[16,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(0, simplex_flux_K2[17,0],yerr=simplex_flux_K2[17,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(0, simplex_flux_K2[18,0],yerr=simplex_flux_K2[18,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(0,10,num=2)
    plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    plt.legend(prop={'size': 12})
    plt.show()

## Position K1
if pos_K1_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex}$ [pixels]", fontsize=20)
    plt.xlabel('r$_{MCMC}$ - r$_{Simplex}$ [pixels]', fontsize=20)
    plt.title("r - K1")
    plt.scatter(mcmc_pos_K1[0,0] - simplex_pos_K1[0,0],mcmc_pos_K1[0,1]-simplex_pos_K1[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_pos_K1[1,0] - simplex_pos_K1[1,0],mcmc_pos_K1[1,1]-simplex_pos_K1[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_pos_K1[2,0]-simplex_pos_K1[2,0],mcmc_pos_K1[2,1]-simplex_pos_K1[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_pos_K1[3,0]- simplex_pos_K1[3,0],mcmc_pos_K1[3,1]-simplex_pos_K1[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_pos_K1[4,0]- simplex_pos_K1[4,0],mcmc_pos_K1[4,1]-simplex_pos_K1[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_pos_K1[5,0]- simplex_pos_K1[5,0],mcmc_pos_K1[5,1]-simplex_pos_K1[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_pos_K1[6,0]- simplex_pos_K1[6,0],mcmc_pos_K1[6,1]-simplex_pos_K1[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_pos_K1[7,0]- simplex_pos_K1[7,0],mcmc_pos_K1[7,1]-simplex_pos_K1[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_pos_K1[8,0]- simplex_pos_K1[8,0],mcmc_pos_K1[8,1]-simplex_pos_K1[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_pos_K1[9,0]- simplex_pos_K1[9,0],mcmc_pos_K1[9,1]-simplex_pos_K1[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_pos_K1[10,0]- simplex_pos_K1[10,0],mcmc_pos_K1[10,1]-simplex_pos_K1[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_pos_K1[11,0]- simplex_pos_K1[11,0],mcmc_pos_K1[11,1]-simplex_pos_K1[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_pos_K1[12,0]- simplex_pos_K1[12,0],mcmc_pos_K1[12,1]-simplex_pos_K1[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_pos_K1[13,0]- simplex_pos_K1[13,0],mcmc_pos_K1[13,1]-simplex_pos_K1[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_pos_K1[14,0]- simplex_pos_K1[14,0],mcmc_pos_K1[14,1]-simplex_pos_K1[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_pos_K1[15,0],simplex_pos_K1[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0- simplex_pos_K1[16,0],simplex_pos_K1[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0- simplex_pos_K1[17,0],simplex_pos_K1[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0- simplex_pos_K1[18,0],simplex_pos_K1[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-1.5,2),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-2,50),lw=0.8,c="black")
    plt.xlim(-1.5,2)
    plt.ylim(-2,50)
    plt.legend(prop={'size': 12})
    plt.show()

## Position K2
if pos_K2_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex} [pixels]$", fontsize=20)
    plt.xlabel('r$_{MCMC}$ - r$_{Simplex}$ [pixels]', fontsize=20)
    plt.title("r - K2")
    plt.scatter(mcmc_pos_K2[0,0]-simplex_pos_K2[0,0],mcmc_pos_K2[0,1]-simplex_pos_K2[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_pos_K2[1,0]-simplex_pos_K2[1,0],mcmc_pos_K2[1,1]-simplex_pos_K2[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_pos_K2[2,0]-simplex_pos_K2[2,0],mcmc_pos_K2[2,1]-simplex_pos_K2[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_pos_K2[3,0]-simplex_pos_K2[3,0],mcmc_pos_K2[3,1]-simplex_pos_K2[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_pos_K2[4,0]-simplex_pos_K2[4,0],mcmc_pos_K2[4,1]-simplex_pos_K2[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_pos_K2[5,0]-simplex_pos_K2[5,0],mcmc_pos_K2[5,1]-simplex_pos_K2[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_pos_K2[6,0]-simplex_pos_K2[6,0],mcmc_pos_K2[6,1]-simplex_pos_K2[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_pos_K2[7,0]-simplex_pos_K2[7,0],mcmc_pos_K2[7,1]-simplex_pos_K2[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_pos_K2[8,0]-simplex_pos_K2[8,0],mcmc_pos_K2[8,1]-simplex_pos_K2[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_pos_K2[9,0]-simplex_pos_K2[9,0],mcmc_pos_K2[9,1]-simplex_pos_K2[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_pos_K2[10,0]-simplex_pos_K2[10,0],mcmc_pos_K2[10,1]-simplex_pos_K2[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_pos_K2[11,0]-simplex_pos_K2[11,0],mcmc_pos_K2[11,1]-simplex_pos_K2[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_pos_K2[12,0]-simplex_pos_K2[12,0],mcmc_pos_K2[12,1]-simplex_pos_K2[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_pos_K2[13,0]-simplex_pos_K2[13,0],mcmc_pos_K2[13,1]-simplex_pos_K2[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_pos_K2[14,0]-simplex_pos_K2[14,0],mcmc_pos_K2[14,1]-simplex_pos_K2[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_pos_K2[15,0],simplex_pos_K2[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0-simplex_pos_K2[16,0],simplex_pos_K2[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0-simplex_pos_K2[17,0],simplex_pos_K2[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0-simplex_pos_K2[18,0],simplex_pos_K2[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-4,6),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-2,60),lw=0.8,c="black")
    plt.xlim(-4,6)
    plt.ylim(-2,60)
    plt.legend(prop={'size': 12})
    plt.show()

## Position-Angle K1
if PA_K1_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex}$ [degrees]", fontsize=20)
    plt.xlabel('$PA_{MCMC} - PA_{Simplex}$ [degrees]', fontsize=20)
    plt.title("PA - K1")
    plt.scatter(mcmc_PA_K1[0,0]-simplex_PA_K1[0,0],mcmc_PA_K1[0,1]-simplex_PA_K1[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_PA_K1[1,0]-simplex_PA_K1[1,0],mcmc_PA_K1[1,1]-simplex_PA_K1[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_PA_K1[2,0]-simplex_PA_K1[2,0],mcmc_PA_K1[2,1]-simplex_PA_K1[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_PA_K1[3,0]-simplex_PA_K1[3,0],mcmc_PA_K1[3,1]-simplex_PA_K1[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_PA_K1[4,0]-simplex_PA_K1[4,0],mcmc_PA_K1[4,1]-simplex_PA_K1[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_PA_K1[5,0]-simplex_PA_K1[5,0],mcmc_PA_K1[5,1]-simplex_PA_K1[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_PA_K1[6,0]-simplex_PA_K1[6,0],mcmc_PA_K1[6,1]-simplex_PA_K1[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_PA_K1[7,0]-simplex_PA_K1[7,0],mcmc_PA_K1[7,1]-simplex_PA_K1[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_PA_K1[8,0]-simplex_PA_K1[8,0],mcmc_PA_K1[8,1]-simplex_PA_K1[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_PA_K1[9,0]-simplex_PA_K1[9,0],mcmc_PA_K1[9,1]-simplex_PA_K1[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_PA_K1[10,0]-simplex_PA_K1[10,0],mcmc_PA_K1[10,1]-simplex_PA_K1[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_PA_K1[11,0]-simplex_PA_K1[11,0],mcmc_PA_K1[11,1]-simplex_PA_K1[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_PA_K1[12,0]-simplex_PA_K1[12,0],mcmc_PA_K1[12,1]-simplex_PA_K1[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_PA_K1[13,0]-simplex_PA_K1[13,0],mcmc_PA_K1[13,1]-simplex_PA_K1[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_PA_K1[14,0]-simplex_PA_K1[14,0],mcmc_PA_K1[14,1]-simplex_PA_K1[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_PA_K1[15,0],simplex_PA_K1[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0-simplex_PA_K1[16,0],simplex_PA_K1[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0-simplex_PA_K1[17,0],simplex_PA_K1[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0-simplex_PA_K1[18,0],simplex_PA_K1[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-130,35),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-2,40),lw=0.8,c="black")
    plt.xlim(-130,35)
    plt.ylim(-2,40)
    plt.legend(prop={'size': 12})
    plt.show()

## Position-Angle K2
if PA_K2_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex}$ [degrees]", fontsize=20)
    plt.xlabel('$PA_{MCMC} - PA_{Simplex}$ [degrees]', fontsize=20)
    plt.title("PA - K2")
    plt.scatter(mcmc_PA_K2[0,0]-simplex_PA_K2[0,0],mcmc_PA_K2[0,1]-simplex_PA_K2[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_PA_K2[1,0]-simplex_PA_K2[1,0],mcmc_PA_K2[1,1]-simplex_PA_K2[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_PA_K2[2,0]-simplex_PA_K2[2,0],mcmc_PA_K2[2,1]-simplex_PA_K2[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_PA_K2[3,0]-simplex_PA_K2[3,0],mcmc_PA_K2[3,1]-simplex_PA_K2[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_PA_K2[4,0]-simplex_PA_K2[4,0],mcmc_PA_K2[4,1]-simplex_PA_K2[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_PA_K2[5,0]-simplex_PA_K2[5,0],mcmc_PA_K2[5,1]-simplex_PA_K2[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_PA_K2[6,0]-simplex_PA_K2[6,0],mcmc_PA_K2[6,1]-simplex_PA_K2[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_PA_K2[7,0]-simplex_PA_K2[7,0],mcmc_PA_K2[7,1]-simplex_PA_K2[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_PA_K2[8,0]-simplex_PA_K2[8,0],mcmc_PA_K2[8,1]-simplex_PA_K2[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_PA_K2[9,0]-simplex_PA_K2[9,0],mcmc_PA_K2[9,1]-simplex_PA_K2[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_PA_K2[10,0]-simplex_PA_K2[10,0],mcmc_PA_K2[10,1]-simplex_PA_K2[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_PA_K2[11,0]-simplex_PA_K2[11,0],mcmc_PA_K2[11,1]-simplex_PA_K2[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_PA_K2[12,0]-simplex_PA_K2[12,0],mcmc_PA_K2[12,1]-simplex_PA_K2[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_PA_K2[13,0]-simplex_PA_K2[13,0],mcmc_PA_K2[13,1]-simplex_PA_K2[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_PA_K2[14,0]-simplex_PA_K2[14,0],mcmc_PA_K2[14,1]-simplex_PA_K2[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_PA_K2[15,0],simplex_PA_K2[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0-simplex_PA_K2[16,0],simplex_PA_K2[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0-simplex_PA_K2[17,0],simplex_PA_K2[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0-simplex_PA_K2[18,0],simplex_PA_K2[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-0.6,1),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-2,55),lw=0.8,c="black")
    plt.xlim(-0.6,1)
    plt.ylim(-2,55)
    plt.legend(prop={'size': 12})
    plt.show()

## Flux K1
if flux_K1_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex}$ [ADU/s]", fontsize=20)
    plt.xlabel('f$_{MCMC}$ - f$_{Simplex}$ [ADU/s]', fontsize=20)
    plt.title("flux - K1")
    plt.scatter(mcmc_flux_K1[0,0]-simplex_flux_K1[0,0],mcmc_flux_K1[0,1]-simplex_flux_K1[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_flux_K1[1,0]-simplex_flux_K1[1,0],mcmc_flux_K1[1,1]-simplex_flux_K1[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_flux_K1[2,0]-simplex_flux_K1[2,0],mcmc_flux_K1[2,1]-simplex_flux_K1[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_flux_K1[3,0]-simplex_flux_K1[3,0],mcmc_flux_K1[3,1]-simplex_flux_K1[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_flux_K1[4,0]-simplex_flux_K1[4,0],mcmc_flux_K1[4,1]-simplex_flux_K1[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_flux_K1[5,0]-simplex_flux_K1[5,0],mcmc_flux_K1[5,1]-simplex_flux_K1[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_flux_K1[6,0]-simplex_flux_K1[6,0],mcmc_flux_K1[6,1]-simplex_flux_K1[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_flux_K1[7,0]-simplex_flux_K1[7,0],mcmc_flux_K1[7,1]-simplex_flux_K1[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_flux_K1[8,0]-simplex_flux_K1[8,0],mcmc_flux_K1[8,1]-simplex_flux_K1[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_flux_K1[9,0]-simplex_flux_K1[9,0],mcmc_flux_K1[9,1]-simplex_flux_K1[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_flux_K1[10,0]-simplex_flux_K1[10,0],mcmc_flux_K1[10,1]-simplex_flux_K1[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_flux_K1[11,0]-simplex_flux_K1[11,0],mcmc_flux_K1[11,1]-simplex_flux_K1[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_flux_K1[12,0]-simplex_flux_K1[12,0],mcmc_flux_K1[12,1]-simplex_flux_K1[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_flux_K1[13,0]-simplex_flux_K1[13,0],mcmc_flux_K1[13,1]-simplex_flux_K1[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_flux_K1[14,0]-simplex_flux_K1[14,0],mcmc_flux_K1[14,1]-simplex_flux_K1[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_flux_K1[15,0],simplex_flux_K1[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0-simplex_flux_K1[16,0],simplex_flux_K1[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0-simplex_flux_K1[17,0],simplex_flux_K1[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0-simplex_flux_K1[18,0],simplex_flux_K1[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-550,50),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-3,25),lw=0.8,c="black")
    plt.xlim(-550,50)
    plt.ylim(-3,25)
    plt.legend(prop={'size': 12})
    plt.show()

## Flux K2
if flux_K2_err == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{MCMC} - \sigma_{Simplex}$ [ADU/s]", fontsize=20)
    plt.xlabel('f$_{MCMC}$ - f$_{Simplex}$ [ADU/s]', fontsize=20)
    plt.title("flux - K2")
    plt.scatter(mcmc_flux_K2[0,0]-simplex_flux_K2[0,0],mcmc_flux_K2[0,1]-simplex_flux_K2[0,1],s=50,c=colors[0],label= "Ad")
    # plt.scatter(mcmc_flux_K2[1,0]-simplex_flux_K2[1,0],mcmc_flux_K2[1,1]-simplex_flux_K2[1,1],s=50,c=colors[1],label= "Ab")
    plt.scatter(mcmc_flux_K2[2,0]-simplex_flux_K2[2,0],mcmc_flux_K2[2,1]-simplex_flux_K2[2,1],s=50,c=colors[2],label= "E")
    plt.scatter(mcmc_flux_K2[3,0]-simplex_flux_K2[3,0],mcmc_flux_K2[3,1]-simplex_flux_K2[3,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_flux_K2[4,0]-simplex_flux_K2[4,0],mcmc_flux_K2[4,1]-simplex_flux_K2[4,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_flux_K2[5,0]-simplex_flux_K2[5,0],mcmc_flux_K2[5,1]-simplex_flux_K2[5,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_flux_K2[6,0]-simplex_flux_K2[6,0],mcmc_flux_K2[6,1]-simplex_flux_K2[6,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_flux_K2[7,0]-simplex_flux_K2[7,0],mcmc_flux_K2[7,1]-simplex_flux_K2[7,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_flux_K2[8,0]-simplex_flux_K2[8,0],mcmc_flux_K2[8,1]-simplex_flux_K2[8,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_flux_K2[9,0]-simplex_flux_K2[9,0],mcmc_flux_K2[9,1]-simplex_flux_K2[9,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_flux_K2[10,0]-simplex_flux_K2[10,0],mcmc_flux_K2[10,1]-simplex_flux_K2[10,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_flux_K2[11,0]-simplex_flux_K2[11,0],mcmc_flux_K2[11,1]-simplex_flux_K2[11,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_flux_K2[12,0]-simplex_flux_K2[12,0],mcmc_flux_K2[12,1]-simplex_flux_K2[12,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_flux_K2[13,0]-simplex_flux_K2[13,0],mcmc_flux_K2[13,1]-simplex_flux_K2[13,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_flux_K2[14,0]-simplex_flux_K2[14,0],mcmc_flux_K2[14,1]-simplex_flux_K2[14,1],s=50,c=colors[14],label= "S12")
    # plt.scatter(0-simplex_flux_K2[15,0],simplex_flux_K2[15,1],s=50,c=colors[15],label= "S13")
    # plt.scatter(0-simplex_flux_K2[16,0],simplex_flux_K2[16,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(0-simplex_flux_K2[17,0],simplex_flux_K2[17,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(0-simplex_flux_K2[18,0],simplex_flux_K2[18,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-400,30),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-1000,20),lw=0.8,c="black")
    plt.xlim(-400,30)
    plt.ylim(-1000,20)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia X
if x_simplex_K1_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('x$_{simplex}$ - x$_{Julia}$', fontsize=20)
    plt.title("x - K1")
    plt.scatter(simplex_x_K1[3,0]- julia_x_K1[0,0],simplex_x_K1[3,1]-julia_x_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_x_K1[4,0]- julia_x_K1[1,0],simplex_x_K1[4,1]-julia_x_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_x_K1[5,0]- julia_x_K1[2,0],simplex_x_K1[5,1]-julia_x_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_x_K1[6,0]- julia_x_K1[3,0],simplex_x_K1[6,1]-julia_x_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_x_K1[7,0]- julia_x_K1[4,0],simplex_x_K1[7,1]-julia_x_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_x_K1[8,0]- julia_x_K1[5,0],simplex_x_K1[8,1]-julia_x_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_x_K1[9,0]- julia_x_K1[6,0],simplex_x_K1[9,1]-julia_x_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_x_K1[10,0]- julia_x_K1[7,0],simplex_x_K1[10,1]-julia_x_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_x_K1[11,0]- julia_x_K1[8,0],simplex_x_K1[11,1]-julia_x_K1[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_x_K1[12,0]- julia_x_K1[9,0],simplex_x_K1[12,1]-julia_x_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_x_K1[13,0]- julia_x_K1[10,0],simplex_x_K1[13,1]-julia_x_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_x_K1[14,0]- julia_x_K1[11,0],simplex_x_K1[14,1]-julia_x_K1[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_x_K1[15,0]- julia_x_K1[12,0],simplex_x_K1[15,1]-julia_x_K1[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_x_K1[16,0]- julia_x_K1[13,0],simplex_x_K1[16,1]-julia_x_K1[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_x_K1[17,0]- julia_x_K1[14,0],simplex_x_K1[17,1]-julia_x_K1[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_x_K1[18,0]- julia_x_K1[15,0],simplex_x_K1[18,1]-julia_x_K1[15,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-0.6,1),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-0.03,0.19),lw=0.8,c="black")
    plt.xlim(-0.6,1)
    plt.ylim(-0.03,0.185)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Y
if y_simplex_K1_julia == True:
    #
    # plt.grid(True)
    # plt.legend()
    # plt.figure(figsize=(12, 9))
    # ax = plt.subplot(111)
    # ax = plt.subplot(111)
    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    # plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    # plt.xticks(fontsize=18)
    # plt.yticks(fontsize=18)
    # plt.ylabel("$x_{Julia}$ ", fontsize=20)
    # plt.xlabel('$x_{Simplex}$', fontsize=20)
    # plt.title("$x$ - K1")
    # plt.errorbar(simplex_x_K1[3,0],simplex_x_K1[3,0]-julia_x_K1[0,0],xerr=simplex_x_K1[3,1],yerr=np.sqrt(simplex_x_K1[3,1]**2+julia_x_K1[0,1]**2),fmt='--o',c=colors[3],label= "S1")
    # plt.errorbar(simplex_x_K1[4,0],simplex_x_K1[4,0]-julia_x_K1[1,0],xerr=simplex_x_K1[4,1],yerr=np.sqrt(simplex_x_K1[4,1]**2+julia_x_K1[1,1]**2),fmt='--o',c=colors[4],label= "S2")
    # plt.errorbar(simplex_x_K1[5,0],simplex_x_K1[5,0]-julia_x_K1[2,0],xerr=simplex_x_K1[5,1],yerr=np.sqrt(simplex_x_K1[5,1]**2+julia_x_K1[2,1]**2),fmt='--o',c=colors[5],label= "S3")
    # plt.errorbar(simplex_x_K1[6,0],simplex_x_K1[6,0]-julia_x_K1[3,0],xerr=simplex_x_K1[6,1],yerr=np.sqrt(simplex_x_K1[6,1]**2+julia_x_K1[3,1]**2),fmt='--o',c=colors[6],label= "S4")
    # plt.errorbar(simplex_x_K1[7,0],simplex_x_K1[7,0]-julia_x_K1[4,0],xerr=simplex_x_K1[7,1],yerr=np.sqrt(simplex_x_K1[7,1]**2+julia_x_K1[4,1]**2),fmt='--o',c=colors[7],label= "S5")
    # plt.errorbar(simplex_x_K1[8,0],simplex_x_K1[8,0]-julia_x_K1[5,0],xerr=simplex_x_K1[8,1],yerr=np.sqrt(simplex_x_K1[8,1]**2+julia_x_K1[5,1]**2),fmt='--o',c=colors[8],label= "S6")
    # plt.errorbar(simplex_x_K1[9,0],simplex_x_K1[9,0]-julia_x_K1[6,0],xerr=simplex_x_K1[9,1],yerr=np.sqrt(simplex_x_K1[9,1]**2+julia_x_K1[6,1]**2),fmt='--o',c=colors[9],label= "S7")
    # plt.errorbar(simplex_x_K1[10,0],simplex_x_K1[10,0]-julia_x_K1[7,0],xerr=simplex_x_K1[10,1],yerr=np.sqrt(simplex_x_K1[10,1]**2+julia_x_K1[7,1]**2),fmt='--o',c=colors[10],label= "S8")
    # plt.errorbar(simplex_x_K1[11,0],simplex_x_K1[11,0]-julia_x_K1[8,0],xerr=simplex_x_K1[11,1],yerr=np.sqrt(simplex_x_K1[11,1]**2+julia_x_K1[8,1]**2),fmt='--o',c=colors[11],label= "S9")
    # plt.errorbar(simplex_x_K1[12,0],simplex_x_K1[12,0]-julia_x_K1[9,0],xerr=simplex_x_K1[12,1],yerr=np.sqrt(simplex_x_K1[12,1]**2+julia_x_K1[9,1]**2),fmt='--o',c=colors[12],label= "S10")
    # plt.errorbar(simplex_x_K1[13,0],simplex_x_K1[13,0]-julia_x_K1[10,0],xerr=simplex_x_K1[13,1],yerr=np.sqrt(simplex_x_K1[13,1]**2+julia_x_K1[10,1]**2),fmt='--o',c=colors[13],label= "S11")
    # plt.errorbar(simplex_x_K1[14,0],simplex_x_K1[14,0]-julia_x_K1[11,0],xerr=simplex_x_K1[14,1],yerr=np.sqrt(simplex_x_K1[14,1]**2+julia_x_K1[11,1]**2),fmt='--o',c=colors[14],label= "S12")
    # plt.errorbar(simplex_x_K1[15,0],simplex_x_K1[15,0]-julia_x_K1[12,0],xerr=simplex_x_K1[15,1],yerr=np.sqrt(simplex_x_K1[15,1]**2+julia_x_K1[12,1]**2),fmt='--o',c=colors[15],label= "S13")
    # plt.errorbar(simplex_x_K1[16,0],simplex_x_K1[16,0]-julia_x_K1[13,0],xerr=simplex_x_K1[16,1],yerr=np.sqrt(simplex_x_K1[16,1]**2+julia_x_K1[13,1]**2),fmt='--o',c=colors[16],label= "S14")
    # # plt.errorbar(simplex_x_K1[17,0],simplex_x_K1[17,0]-julia_x_K1[14,0],xerr=0,yerr=np.sqrt(simplex_x_K1[17,1]**2+julia_x_K1[14,1]**2),fmt='--o',c=colors[17],label= "S15")
    # # plt.errorbar(simplex_x_K1[18,0],simplex_x_K1[18,0]-julia_x_K1[15,0],xerr=0,yerr=np.sqrt(simplex_x_K1[18,1]**2+julia_x_K1[15,1]**2),fmt='--o',c=colors[18],label= "S16")
    # # x = np.linspace(0,900,num=2)
    # # plt.plot(x,x,lw=2.8,c="royalblue",alpha=0.5)
    # plt.legend(prop={'size': 12})
    # plt.show()

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('y$_{simplex}$ - y$_{Julia}$', fontsize=20)
    plt.title("y - K1")
    plt.scatter(simplex_y_K1[3,0]- julia_y_K1[0,0],simplex_y_K1[3,1]-julia_y_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_y_K1[4,0]- julia_y_K1[1,0],simplex_y_K1[4,1]-julia_y_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_y_K1[5,0]- julia_y_K1[2,0],simplex_y_K1[5,1]-julia_y_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_y_K1[6,0]- julia_y_K1[3,0],simplex_y_K1[6,1]-julia_y_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_y_K1[7,0]- julia_y_K1[4,0],simplex_y_K1[7,1]-julia_y_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_y_K1[8,0]- julia_y_K1[5,0],simplex_y_K1[8,1]-julia_y_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_y_K1[9,0]- julia_y_K1[6,0],simplex_y_K1[9,1]-julia_y_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_y_K1[10,0]- julia_y_K1[7,0],simplex_y_K1[10,1]-julia_y_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_y_K1[11,0]- julia_y_K1[8,0],simplex_y_K1[11,1]-julia_y_K1[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_y_K1[12,0]- julia_y_K1[9,0],simplex_y_K1[12,1]-julia_y_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_y_K1[13,0]- julia_y_K1[10,0],simplex_y_K1[13,1]-julia_y_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_y_K1[14,0]- julia_y_K1[11,0],simplex_y_K1[14,1]-julia_y_K1[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_y_K1[15,0]- julia_y_K1[12,0],simplex_y_K1[15,1]-julia_y_K1[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_y_K1[16,0]- julia_y_K1[13,0],simplex_y_K1[16,1]-julia_y_K1[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_y_K1[17,0]- julia_y_K1[14,0],simplex_y_K1[17,1]-julia_y_K1[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_y_K1[18,0]- julia_y_K1[15,0],simplex_y_K1[18,1]-julia_y_K1[15,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-0.15,0.55),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-0.1,0.25),lw=0.8,c="black")
    plt.xlim(-0.15,0.55)
    plt.ylim(-0.1,0.25)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia X
if x_simplex_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('x$_{simplex}$ - x$_{Julia}$', fontsize=20)
    plt.title("x - K2")
    plt.scatter(simplex_x_K2[3,0]- julia_x_K2[0,0],simplex_x_K2[3,1]-julia_x_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_x_K2[4,0]- julia_x_K2[1,0],simplex_x_K2[4,1]-julia_x_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_x_K2[5,0]- julia_x_K2[2,0],simplex_x_K2[5,1]-julia_x_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_x_K2[6,0]- julia_x_K2[3,0],simplex_x_K2[6,1]-julia_x_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_x_K2[7,0]- julia_x_K2[4,0],simplex_x_K2[7,1]-julia_x_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_x_K2[8,0]- julia_x_K2[5,0],simplex_x_K2[8,1]-julia_x_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_x_K2[9,0]- julia_x_K2[6,0],simplex_x_K2[9,1]-julia_x_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_x_K2[10,0]- julia_x_K2[7,0],simplex_x_K2[10,1]-julia_x_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_x_K2[11,0]- julia_x_K2[8,0],simplex_x_K2[11,1]-julia_x_K2[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_x_K2[12,0]- julia_x_K2[9,0],simplex_x_K2[12,1]-julia_x_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_x_K2[13,0]- julia_x_K2[10,0],simplex_x_K2[13,1]-julia_x_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_x_K2[14,0]- julia_x_K2[11,0],simplex_x_K2[14,1]-julia_x_K2[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_x_K2[15,0]- julia_x_K2[12,0],simplex_x_K2[15,1]-julia_x_K2[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_x_K2[16,0]- julia_x_K2[13,0],simplex_x_K2[16,1]-julia_x_K2[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_x_K2[17,0]- julia_x_K2[14,0],simplex_x_K2[17,1]-julia_x_K2[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_x_K2[18,0]- julia_x_K2[15,0],simplex_x_K2[18,1]-julia_x_K2[15,1],s=50,c=colors[18],label= "S16")
    # plt.plot(np.linspace(-1.5,1.5),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-1,13),lw=0.8,c="black")
    # plt.xlim(-1.5,1.5)
    # plt.ylim(-1,13)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Y
if y_simplex_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('y$_{simplex}$ - y$_{Julia}$', fontsize=20)
    plt.title("y - K2")
    plt.scatter(simplex_y_K2[3,0]- julia_y_K2[0,0],simplex_y_K2[3,1]-julia_y_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_y_K2[4,0]- julia_y_K2[1,0],simplex_y_K2[4,1]-julia_y_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_y_K2[5,0]- julia_y_K2[2,0],simplex_y_K2[5,1]-julia_y_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_y_K2[6,0]- julia_y_K2[3,0],simplex_y_K2[6,1]-julia_y_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_y_K2[7,0]- julia_y_K2[4,0],simplex_y_K2[7,1]-julia_y_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_y_K2[8,0]- julia_y_K2[5,0],simplex_y_K2[8,1]-julia_y_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_y_K2[9,0]- julia_y_K2[6,0],simplex_y_K2[9,1]-julia_y_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_y_K2[10,0]- julia_y_K2[7,0],simplex_y_K2[10,1]-julia_y_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_y_K2[11,0]- julia_y_K2[8,0],simplex_y_K2[11,1]-julia_y_K2[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_y_K2[12,0]- julia_y_K2[9,0],simplex_y_K2[12,1]-julia_y_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_y_K2[13,0]- julia_y_K2[10,0],simplex_y_K2[13,1]-julia_y_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_y_K2[14,0]- julia_y_K2[11,0],simplex_y_K2[14,1]-julia_y_K2[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_y_K2[15,0]- julia_y_K2[12,0],simplex_y_K2[15,1]-julia_y_K2[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_y_K2[16,0]- julia_y_K2[13,0],simplex_y_K2[16,1]-julia_y_K2[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_y_K2[17,0]- julia_y_K2[14,0],simplex_y_K2[17,1]-julia_y_K2[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_y_K2[18,0]- julia_y_K2[15,0],simplex_y_K2[18,1]-julia_y_K2[15,1],s=50,c=colors[18],label= "S16")
    # plt.plot(np.linspace(-1.2,3.2),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-5,100),lw=0.8,c="black")
    # plt.xlim(-1.2,3.2)
    # plt.ylim(-5,100)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K1
if dmag_simplex_K1_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{simplex}$ - $\Delta m_{Julia}$', fontsize=20)
    plt.title("$\Delta$m - K1")
    plt.scatter(simplex_dmag_K1[3,0]- julia_dmag_K1[0,0],simplex_dmag_K1[3,1]-julia_dmag_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_dmag_K1[4,0]- julia_dmag_K1[1,0],simplex_dmag_K1[4,1]-julia_dmag_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_dmag_K1[5,0]- julia_dmag_K1[2,0],simplex_dmag_K1[5,1]-julia_dmag_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_dmag_K1[6,0]- julia_dmag_K1[3,0],simplex_dmag_K1[6,1]-julia_dmag_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_dmag_K1[7,0]- julia_dmag_K1[4,0],simplex_dmag_K1[7,1]-julia_dmag_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_dmag_K1[8,0]- julia_dmag_K1[5,0],simplex_dmag_K1[8,1]-julia_dmag_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_dmag_K1[9,0]- julia_dmag_K1[6,0],simplex_dmag_K1[9,1]-julia_dmag_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_dmag_K1[10,0]- julia_dmag_K1[7,0],simplex_dmag_K1[10,1]-julia_dmag_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_dmag_K1[11,0]- julia_dmag_K1[8,0],simplex_dmag_K1[11,1]-julia_dmag_K1[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_dmag_K1[12,0]- julia_dmag_K1[9,0],simplex_dmag_K1[12,1]-julia_dmag_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_dmag_K1[13,0]- julia_dmag_K1[10,0],simplex_dmag_K1[13,1]-julia_dmag_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_dmag_K1[14,0]- julia_dmag_K1[11,0],simplex_dmag_K1[14,1]-julia_dmag_K1[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_dmag_K1[15,0]- julia_dmag_K1[12,0],simplex_dmag_K1[15,1]-julia_dmag_K1[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_dmag_K1[16,0]- julia_dmag_K1[13,0],simplex_dmag_K1[16,1]-julia_dmag_K1[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_dmag_K1[17,0]- julia_dmag_K1[14,0],simplex_dmag_K1[17,1]-julia_dmag_K1[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_dmag_K1[18,0]- julia_dmag_K1[15,0],simplex_dmag_K1[18,1]-julia_dmag_K1[15,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-0.8,0.2),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-0.13,0.04),lw=0.8,c="black")
    plt.xlim(-0.8,0.2)
    plt.ylim(-0.13,0.04)
    plt.legend(prop={'size': 12})
    plt.show()

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{simplex}$', fontsize=20)
    plt.title("$\Delta m$ - K1")
    plt.errorbar(simplex_dmag_K1[3,0],julia_dmag_K1[0,0],xerr=simplex_dmag_K1[3,1],yerr=julia_dmag_K1[0,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(simplex_dmag_K1[4,0],julia_dmag_K1[1,0],xerr=simplex_dmag_K1[4,1],yerr=julia_dmag_K1[1,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(simplex_dmag_K1[5,0],julia_dmag_K1[2,0],xerr=simplex_dmag_K1[5,1],yerr=julia_dmag_K1[2,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(simplex_dmag_K1[6,0],julia_dmag_K1[3,0],xerr=simplex_dmag_K1[6,1],yerr=julia_dmag_K1[3,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(simplex_dmag_K1[7,0],julia_dmag_K1[4,0],xerr=simplex_dmag_K1[7,1],yerr=julia_dmag_K1[4,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(simplex_dmag_K1[8,0],julia_dmag_K1[5,0],xerr=simplex_dmag_K1[8,1],yerr=julia_dmag_K1[5,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(simplex_dmag_K1[9,0],julia_dmag_K1[6,0],xerr=simplex_dmag_K1[9,1],yerr=julia_dmag_K1[6,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(simplex_dmag_K1[10,0],julia_dmag_K1[7,0],xerr=simplex_dmag_K1[10,1],yerr=julia_dmag_K1[7,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(simplex_dmag_K1[11,0],julia_dmag_K1[8,0],xerr=simplex_dmag_K1[11,1],yerr=julia_dmag_K1[8,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(simplex_dmag_K1[12,0],julia_dmag_K1[9,0],xerr=simplex_dmag_K1[12,1],yerr=julia_dmag_K1[9,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(simplex_dmag_K1[13,0],julia_dmag_K1[10,0],xerr=simplex_dmag_K1[13,1],yerr=julia_dmag_K1[10,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(simplex_dmag_K1[14,0],julia_dmag_K1[11,0],xerr=simplex_dmag_K1[14,1],yerr=julia_dmag_K1[11,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(simplex_dmag_K1[15,0],julia_dmag_K1[12,0],xerr=simplex_dmag_K1[15,1],yerr=julia_dmag_K1[12,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(simplex_dmag_K1[16,0],julia_dmag_K1[13,0],xerr=simplex_dmag_K1[16,1],yerr=julia_dmag_K1[13,1],fmt='--o',c=colors[16],label= "S14")
    # plt.errorbar(simplex_dmag_K1[17,0],simplex_dmag_K1[17,0]-julia_dmag_K1[14,0],xerr=0,yerr=np.sqrt(simplex_dmag_K1[17,1]**2+julia_dmag_K1[14,1]**2),fmt='--o',c=colors[17],label= "S15")
    # plt.errorbar(simplex_dmag_K1[18,0],simplex_dmag_K1[18,0]-julia_dmag_K1[15,0],xerr=0,yerr=np.sqrt(simplex_dmag_K1[18,1]**2+julia_dmag_K1[15,1]**2),fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K2
if dmag_simplex_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{simplex} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{simplex}$ - $\Delta m_{Julia}$', fontsize=20)
    plt.title("$\Delta$m - K2")
    plt.scatter(simplex_dmag_K2[3,0]- julia_dmag_K2[0,0],simplex_dmag_K2[3,1]-julia_dmag_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(simplex_dmag_K2[4,0]- julia_dmag_K2[1,0],simplex_dmag_K2[4,1]-julia_dmag_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(simplex_dmag_K2[5,0]- julia_dmag_K2[2,0],simplex_dmag_K2[5,1]-julia_dmag_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(simplex_dmag_K2[6,0]- julia_dmag_K2[3,0],simplex_dmag_K2[6,1]-julia_dmag_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(simplex_dmag_K2[7,0]- julia_dmag_K2[4,0],simplex_dmag_K2[7,1]-julia_dmag_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(simplex_dmag_K2[8,0]- julia_dmag_K2[5,0],simplex_dmag_K2[8,1]-julia_dmag_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(simplex_dmag_K2[9,0]- julia_dmag_K2[6,0],simplex_dmag_K2[9,1]-julia_dmag_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(simplex_dmag_K2[10,0]- julia_dmag_K2[7,0],simplex_dmag_K2[10,1]-julia_dmag_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(simplex_dmag_K2[11,0]- julia_dmag_K2[8,0],simplex_dmag_K2[11,1]-julia_dmag_K2[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(simplex_dmag_K2[12,0]- julia_dmag_K2[9,0],simplex_dmag_K2[12,1]-julia_dmag_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(simplex_dmag_K2[13,0]- julia_dmag_K2[10,0],simplex_dmag_K2[13,1]-julia_dmag_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(simplex_dmag_K2[14,0]- julia_dmag_K2[11,0],simplex_dmag_K2[14,1]-julia_dmag_K2[11,1],s=50,c=colors[14],label= "S12")
    plt.scatter(simplex_dmag_K2[15,0]- julia_dmag_K2[12,0],simplex_dmag_K2[15,1]-julia_dmag_K2[12,1],s=50,c=colors[15],label= "S13")
    plt.scatter(simplex_dmag_K2[16,0]- julia_dmag_K2[13,0],simplex_dmag_K2[16,1]-julia_dmag_K2[13,1],s=50,c=colors[16],label= "S14")
    # plt.scatter(simplex_dmag_K2[17,0]- julia_dmag_K2[14,0],simplex_dmag_K2[17,1]-julia_dmag_K2[14,1],s=50,c=colors[17],label= "S15")
    # plt.scatter(simplex_dmag_K2[18,0]- julia_dmag_K2[15,0],simplex_dmag_K2[18,1]-julia_dmag_K2[15,1],s=50,c=colors[18],label= "S16")
    plt.plot(np.linspace(-0.9,0),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-1.8,0.),lw=0.8,c="black")
    plt.xlim(-0.9,0)
    plt.ylim(-1.8,0)
    plt.legend(prop={'size': 12})
    plt.show()

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{simplex}$', fontsize=20)
    plt.title("$\Delta m$ - K2")
    plt.errorbar(simplex_dmag_K2[3,0],julia_dmag_K2[0,0],xerr=simplex_dmag_K2[3,1],yerr=julia_dmag_K2[0,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(simplex_dmag_K2[4,0],julia_dmag_K2[1,0],xerr=simplex_dmag_K2[4,1],yerr=julia_dmag_K2[1,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(simplex_dmag_K2[5,0],julia_dmag_K2[2,0],xerr=simplex_dmag_K2[5,1],yerr=julia_dmag_K2[2,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(simplex_dmag_K2[6,0],julia_dmag_K2[3,0],xerr=simplex_dmag_K2[6,1],yerr=julia_dmag_K2[3,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(simplex_dmag_K2[7,0],julia_dmag_K2[4,0],xerr=simplex_dmag_K2[7,1],yerr=julia_dmag_K2[4,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(simplex_dmag_K2[8,0],julia_dmag_K2[5,0],xerr=simplex_dmag_K2[8,1],yerr=julia_dmag_K2[5,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(simplex_dmag_K2[9,0],julia_dmag_K2[6,0],xerr=simplex_dmag_K2[9,1],yerr=julia_dmag_K2[6,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(simplex_dmag_K2[10,0],julia_dmag_K2[7,0],xerr=simplex_dmag_K2[10,1],yerr=julia_dmag_K2[7,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(simplex_dmag_K2[11,0],julia_dmag_K2[8,0],xerr=simplex_dmag_K2[11,1],yerr=julia_dmag_K2[8,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(simplex_dmag_K2[12,0],julia_dmag_K2[9,0],xerr=simplex_dmag_K2[12,1],yerr=julia_dmag_K2[9,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(simplex_dmag_K2[13,0],julia_dmag_K2[10,0],xerr=simplex_dmag_K2[13,1],yerr=julia_dmag_K2[10,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(simplex_dmag_K2[14,0],julia_dmag_K2[11,0],xerr=simplex_dmag_K2[14,1],yerr=julia_dmag_K2[11,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(simplex_dmag_K2[15,0],julia_dmag_K2[12,0],xerr=simplex_dmag_K2[15,1],yerr=julia_dmag_K2[12,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(simplex_dmag_K2[16,0],julia_dmag_K2[13,0],xerr=simplex_dmag_K2[16,1],yerr=julia_dmag_K2[13,1],fmt='--o',c=colors[16],label= "S14")
    # plt.errorbar(simplex_dmag_K1[17,0],simplex_dmag_K1[17,0]-julia_dmag_K1[14,0],xerr=0,yerr=np.sqrt(simplex_dmag_K1[17,1]**2+julia_dmag_K1[14,1]**2),fmt='--o',c=colors[17],label= "S15")
    # plt.errorbar(simplex_dmag_K1[18,0],simplex_dmag_K1[18,0]-julia_dmag_K1[15,0],xerr=0,yerr=np.sqrt(simplex_dmag_K1[18,1]**2+julia_dmag_K1[15,1]**2),fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia X
if x_mcmc_K1_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('x$_{mcmc}$ - x$_{Julia}$', fontsize=20)
    plt.title("x - K1")
    plt.scatter(mcmc_x_K1[3,0]- julia_x_K1[0,0],mcmc_x_K1[3,1]-julia_x_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_x_K1[4,0]- julia_x_K1[1,0],mcmc_x_K1[4,1]-julia_x_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_x_K1[5,0]- julia_x_K1[2,0],mcmc_x_K1[5,1]-julia_x_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_x_K1[6,0]- julia_x_K1[3,0],mcmc_x_K1[6,1]-julia_x_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_x_K1[7,0]- julia_x_K1[4,0],mcmc_x_K1[7,1]-julia_x_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_x_K1[8,0]- julia_x_K1[5,0],mcmc_x_K1[8,1]-julia_x_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_x_K1[9,0]- julia_x_K1[6,0],mcmc_x_K1[9,1]-julia_x_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_x_K1[10,0]- julia_x_K1[7,0],mcmc_x_K1[10,1]-julia_x_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_x_K1[11,0]- julia_x_K1[8,0],mcmc_x_K1[11,1]-julia_x_K1[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_x_K1[12,0]- julia_x_K1[9,0],mcmc_x_K1[12,1]-julia_x_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_x_K1[13,0]- julia_x_K1[10,0],mcmc_x_K1[13,1]-julia_x_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_x_K1[14,0]- julia_x_K1[11,0],mcmc_x_K1[14,1]-julia_x_K1[11,1],s=50,c=colors[14],label= "S12")
    #
    # plt.plot(np.linspace(-400,70),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-5,110),lw=0.8,c="black")
    # plt.xlim(-400,70)
    # plt.ylim(-5,110)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Y
if y_mcmc_K1_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('y$_{mcmc}$ - y$_{Julia}$', fontsize=20)
    plt.title("y - K1")
    plt.scatter(mcmc_y_K1[3,0]- julia_y_K1[0,0],mcmc_y_K1[3,1]-julia_y_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_y_K1[4,0]- julia_y_K1[1,0],mcmc_y_K1[4,1]-julia_y_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_y_K1[5,0]- julia_y_K1[2,0],mcmc_y_K1[5,1]-julia_y_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_y_K1[6,0]- julia_y_K1[3,0],mcmc_y_K1[6,1]-julia_y_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_y_K1[7,0]- julia_y_K1[4,0],mcmc_y_K1[7,1]-julia_y_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_y_K1[8,0]- julia_y_K1[5,0],mcmc_y_K1[8,1]-julia_y_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_y_K1[9,0]- julia_y_K1[6,0],mcmc_y_K1[9,1]-julia_y_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_y_K1[10,0]- julia_y_K1[7,0],mcmc_y_K1[10,1]-julia_y_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_y_K1[11,0]- julia_y_K1[8,0],mcmc_y_K1[11,1]-julia_y_K1[8,1],s=50,c=colors[11],label= "S9")
    # plt.scatter(mcmc_y_K1[12,0]- julia_y_K1[9,0],mcmc_y_K1[12,1]-julia_y_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_y_K1[13,0]- julia_y_K1[10,0],mcmc_y_K1[13,1]-julia_y_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_y_K1[14,0]- julia_y_K1[11,0],mcmc_y_K1[14,1]-julia_y_K1[11,1],s=50,c=colors[14],label= "S12")
    # plt.plot(np.linspace(-105,5),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-2,15),lw=0.8,c="black")
    # plt.xlim(-105,5)
    # plt.ylim(-2,15)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia X
if x_mcmc_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('x$_{mcmc}$ - x$_{Julia}$', fontsize=20)
    plt.title("x - K2")
    plt.scatter(mcmc_x_K2[3,0]- julia_x_K2[0,0],mcmc_x_K2[3,1]-julia_x_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_x_K2[4,0]- julia_x_K2[1,0],mcmc_x_K2[4,1]-julia_x_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_x_K2[5,0]- julia_x_K2[2,0],mcmc_x_K2[5,1]-julia_x_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_x_K2[6,0]- julia_x_K2[3,0],mcmc_x_K2[6,1]-julia_x_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_x_K2[7,0]- julia_x_K2[4,0],mcmc_x_K2[7,1]-julia_x_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_x_K2[8,0]- julia_x_K2[5,0],mcmc_x_K2[8,1]-julia_x_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_x_K2[9,0]- julia_x_K2[6,0],mcmc_x_K2[9,1]-julia_x_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_x_K2[10,0]- julia_x_K2[7,0],mcmc_x_K2[10,1]-julia_x_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_x_K2[11,0]- julia_x_K2[8,0],mcmc_x_K2[11,1]-julia_x_K2[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_x_K2[12,0]- julia_x_K2[9,0],mcmc_x_K2[12,1]-julia_x_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_x_K2[13,0]- julia_x_K2[10,0],mcmc_x_K2[13,1]-julia_x_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_x_K2[14,0]- julia_x_K2[11,0],mcmc_x_K2[14,1]-julia_x_K2[11,1],s=50,c=colors[14],label= "S12")
    # plt.plot(np.linspace(-1,7),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-1,130),lw=0.8,c="black")
    # plt.xlim(-1,7)
    # plt.ylim(-1,130)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Y
if y_mcmc_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('y$_{mcmc}$ - y$_{Julia}$', fontsize=20)
    plt.title("y - K2")
    plt.scatter(mcmc_y_K2[3,0]- julia_y_K2[0,0],mcmc_y_K2[3,1]-julia_y_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_y_K2[4,0]- julia_y_K2[1,0],mcmc_y_K2[4,1]-julia_y_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_y_K2[5,0]- julia_y_K2[2,0],mcmc_y_K2[5,1]-julia_y_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_y_K2[6,0]- julia_y_K2[3,0],mcmc_y_K2[6,1]-julia_y_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_y_K2[7,0]- julia_y_K2[4,0],mcmc_y_K2[7,1]-julia_y_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_y_K2[8,0]- julia_y_K2[5,0],mcmc_y_K2[8,1]-julia_y_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_y_K2[9,0]- julia_y_K2[6,0],mcmc_y_K2[9,1]-julia_y_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_y_K2[10,0]- julia_y_K2[7,0],mcmc_y_K2[10,1]-julia_y_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_y_K2[11,0]- julia_y_K2[8,0],mcmc_y_K2[11,1]-julia_y_K2[8,1],s=50,c=colors[11],label= "S9")
    # plt.scatter(mcmc_y_K2[12,0]- julia_y_K2[9,0],mcmc_y_K2[12,1]-julia_y_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_y_K2[13,0]- julia_y_K2[10,0],mcmc_y_K2[13,1]-julia_y_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_y_K2[14,0]- julia_y_K2[11,0],mcmc_y_K2[14,1]-julia_y_K2[11,1],s=50,c=colors[14],label= "S12")
    # plt.plot(np.linspace(-0.5,1.2),np.zeros(50),lw=0.8,c="black")
    # plt.plot(np.zeros(50),np.linspace(-5,55),lw=0.8,c="black")
    # plt.xlim(-0.5,1.2)
    # plt.ylim(-5,55)
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K1
if dmag_mcmc_K1_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{mcmc}$ - $\Delta m_{Julia}$', fontsize=20)
    plt.title("$\Delta$m - K1")
    plt.scatter(mcmc_dmag_K1[3,0]- julia_dmag_K1[0,0],mcmc_dmag_K1[3,1]-julia_dmag_K1[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_dmag_K1[4,0]- julia_dmag_K1[1,0],mcmc_dmag_K1[4,1]-julia_dmag_K1[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_dmag_K1[5,0]- julia_dmag_K1[2,0],mcmc_dmag_K1[5,1]-julia_dmag_K1[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_dmag_K1[6,0]- julia_dmag_K1[3,0],mcmc_dmag_K1[6,1]-julia_dmag_K1[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_dmag_K1[7,0]- julia_dmag_K1[4,0],mcmc_dmag_K1[7,1]-julia_dmag_K1[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_dmag_K1[8,0]- julia_dmag_K1[5,0],mcmc_dmag_K1[8,1]-julia_dmag_K1[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_dmag_K1[9,0]- julia_dmag_K1[6,0],mcmc_dmag_K1[9,1]-julia_dmag_K1[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_dmag_K1[10,0]- julia_dmag_K1[7,0],mcmc_dmag_K1[10,1]-julia_dmag_K1[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_dmag_K1[11,0]- julia_dmag_K1[8,0],mcmc_dmag_K1[11,1]-julia_dmag_K1[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_dmag_K1[12,0]- julia_dmag_K1[9,0],mcmc_dmag_K1[12,1]-julia_dmag_K1[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_dmag_K1[13,0]- julia_dmag_K1[10,0],mcmc_dmag_K1[13,1]-julia_dmag_K1[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_dmag_K1[14,0]- julia_dmag_K1[11,0],mcmc_dmag_K1[14,1]-julia_dmag_K1[11,1],s=50,c=colors[14],label= "S12")
    plt.plot(np.linspace(-0.7,0.2),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-0.10,0.),lw=0.8,c="black")
    plt.xlim(-0.7,0.2)
    plt.ylim(-0.10,0.)
    plt.legend(prop={'size': 12})
    plt.show()

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{MCMC}$', fontsize=20)
    plt.title("$\Delta m$ - K1")
    plt.errorbar(mcmc_dmag_K1[3,0],julia_dmag_K1[0,0],xerr=mcmc_dmag_K1[3,1],yerr=julia_dmag_K1[0,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_dmag_K1[4,0],julia_dmag_K1[1,0],xerr=mcmc_dmag_K1[4,1],yerr=julia_dmag_K1[1,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_dmag_K1[5,0],julia_dmag_K1[2,0],xerr=mcmc_dmag_K1[5,1],yerr=julia_dmag_K1[2,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_dmag_K1[6,0],julia_dmag_K1[3,0],xerr=mcmc_dmag_K1[6,1],yerr=julia_dmag_K1[3,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_dmag_K1[7,0],julia_dmag_K1[4,0],xerr=mcmc_dmag_K1[7,1],yerr=julia_dmag_K1[4,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_dmag_K1[8,0],julia_dmag_K1[5,0],xerr=mcmc_dmag_K1[8,1],yerr=julia_dmag_K1[5,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_dmag_K1[9,0],julia_dmag_K1[6,0],xerr=mcmc_dmag_K1[9,1],yerr=julia_dmag_K1[6,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_dmag_K1[10,0],julia_dmag_K1[7,0],xerr=mcmc_dmag_K1[10,1],yerr=julia_dmag_K1[7,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_dmag_K1[11,0],julia_dmag_K1[8,0],xerr=mcmc_dmag_K1[11,1],yerr=julia_dmag_K1[8,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_dmag_K1[12,0],julia_dmag_K1[9,0],xerr=mcmc_dmag_K1[12,1],yerr=julia_dmag_K1[9,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_dmag_K1[13,0],julia_dmag_K1[10,0],xerr=mcmc_dmag_K1[13,1],yerr=julia_dmag_K1[10,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_dmag_K1[14,0],julia_dmag_K1[11,0],xerr=mcmc_dmag_K1[14,1],yerr=julia_dmag_K1[11,1],fmt='--o',c=colors[14],label= "S12")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K2
if dmag_mcmc_K2_julia == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\sigma_{mcmc} - \sigma_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{mcmc}$ - $\Delta m_{Julia}$', fontsize=20)
    plt.title("$\Delta$m - K2")
    plt.scatter(mcmc_dmag_K2[3,0]- julia_dmag_K2[0,0],mcmc_dmag_K2[3,1]-julia_dmag_K2[0,1],s=50,c=colors[3],label= "S1")
    plt.scatter(mcmc_dmag_K2[4,0]- julia_dmag_K2[1,0],mcmc_dmag_K2[4,1]-julia_dmag_K2[1,1],s=50,c=colors[4],label= "S2")
    plt.scatter(mcmc_dmag_K2[5,0]- julia_dmag_K2[2,0],mcmc_dmag_K2[5,1]-julia_dmag_K2[2,1],s=50,c=colors[5],label= "S3")
    plt.scatter(mcmc_dmag_K2[6,0]- julia_dmag_K2[3,0],mcmc_dmag_K2[6,1]-julia_dmag_K2[3,1],s=50,c=colors[6],label= "S4")
    plt.scatter(mcmc_dmag_K2[7,0]- julia_dmag_K2[4,0],mcmc_dmag_K2[7,1]-julia_dmag_K2[4,1],s=50,c=colors[7],label= "S5")
    plt.scatter(mcmc_dmag_K2[8,0]- julia_dmag_K2[5,0],mcmc_dmag_K2[8,1]-julia_dmag_K2[5,1],s=50,c=colors[8],label= "S6")
    plt.scatter(mcmc_dmag_K2[9,0]- julia_dmag_K2[6,0],mcmc_dmag_K2[9,1]-julia_dmag_K2[6,1],s=50,c=colors[9],label= "S7")
    plt.scatter(mcmc_dmag_K2[10,0]- julia_dmag_K2[7,0],mcmc_dmag_K2[10,1]-julia_dmag_K2[7,1],s=50,c=colors[10],label= "S8")
    plt.scatter(mcmc_dmag_K2[11,0]- julia_dmag_K2[8,0],mcmc_dmag_K2[11,1]-julia_dmag_K2[8,1],s=50,c=colors[11],label= "S9")
    plt.scatter(mcmc_dmag_K2[12,0]- julia_dmag_K2[9,0],mcmc_dmag_K2[12,1]-julia_dmag_K2[9,1],s=50,c=colors[12],label= "S10")
    plt.scatter(mcmc_dmag_K2[13,0]- julia_dmag_K2[10,0],mcmc_dmag_K2[13,1]-julia_dmag_K2[10,1],s=50,c=colors[13],label= "S11")
    plt.scatter(mcmc_dmag_K2[14,0]- julia_dmag_K2[11,0],mcmc_dmag_K2[14,1]-julia_dmag_K2[11,1],s=50,c=colors[14],label= "S12")
    plt.plot(np.linspace(-1.0,0.65),np.zeros(50),lw=0.8,c="black")
    plt.plot(np.zeros(50),np.linspace(-0.11,0.045),lw=0.8,c="black")
    plt.xlim(-1.0,0.65)
    plt.ylim(-0.11,0.045)
    plt.legend(prop={'size': 12})
    plt.show()

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{Julia}$", fontsize=20)
    plt.xlabel('$\Delta m_{MCMC}$', fontsize=20)
    plt.title("$\Delta m$ - K2")
    plt.errorbar(mcmc_dmag_K2[3,0],julia_dmag_K2[0,0],xerr=mcmc_dmag_K2[3,1],yerr=julia_dmag_K2[0,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_dmag_K2[4,0],julia_dmag_K2[1,0],xerr=mcmc_dmag_K2[4,1],yerr=julia_dmag_K2[1,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_dmag_K2[5,0],julia_dmag_K2[2,0],xerr=mcmc_dmag_K2[5,1],yerr=julia_dmag_K2[2,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_dmag_K2[6,0],julia_dmag_K2[3,0],xerr=mcmc_dmag_K2[6,1],yerr=julia_dmag_K2[3,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_dmag_K2[7,0],julia_dmag_K2[4,0],xerr=mcmc_dmag_K2[7,1],yerr=julia_dmag_K2[4,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_dmag_K2[8,0],julia_dmag_K2[5,0],xerr=mcmc_dmag_K2[8,1],yerr=julia_dmag_K2[5,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_dmag_K2[9,0],julia_dmag_K2[6,0],xerr=mcmc_dmag_K2[9,1],yerr=julia_dmag_K2[6,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_dmag_K2[10,0],julia_dmag_K2[7,0],xerr=mcmc_dmag_K2[10,1],yerr=julia_dmag_K2[7,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_dmag_K2[11,0],julia_dmag_K2[8,0],xerr=mcmc_dmag_K2[11,1],yerr=julia_dmag_K2[8,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_dmag_K2[12,0],julia_dmag_K2[9,0],xerr=mcmc_dmag_K2[12,1],yerr=julia_dmag_K2[9,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_dmag_K2[13,0],julia_dmag_K2[10,0],xerr=mcmc_dmag_K2[13,1],yerr=julia_dmag_K2[10,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_dmag_K2[14,0],julia_dmag_K2[11,0],xerr=mcmc_dmag_K2[14,1],yerr=julia_dmag_K2[11,1],fmt='--o',c=colors[14],label= "S12")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K1
if dmag_mcmc_K1_K2 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{MCMC}K1$", fontsize=20)
    plt.xlabel('$\Delta m_{MCMC}K1 - \Delta m_{MCMC}K2$', fontsize=20)
    plt.title("K1 vs. K1-K2")
    plt.errorbar(mcmc_dmag_K1[0,0]-mcmc_dmag_K2[0,0],mcmc_dmag_K1[0,0],xerr=mcmc_dmag_K1[0,1],yerr=mcmc_dmag_K1[0,1],fmt='--o',c=colors[0],label= "Ad")
    plt.errorbar(mcmc_dmag_K1[1,0]-mcmc_dmag_K2[1,0],mcmc_dmag_K1[1,0],xerr=mcmc_dmag_K1[1,1],yerr=mcmc_dmag_K1[1,1],fmt='--o',c=colors[1],label= "Ab")
    plt.errorbar(mcmc_dmag_K1[2,0]-mcmc_dmag_K2[2,0],mcmc_dmag_K1[2,0],xerr=mcmc_dmag_K1[2,1],yerr=mcmc_dmag_K1[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_dmag_K1[3,0]-mcmc_dmag_K2[3,0],mcmc_dmag_K1[3,0],xerr=mcmc_dmag_K1[3,1],yerr=mcmc_dmag_K1[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_dmag_K1[4,0]-mcmc_dmag_K2[4,0],mcmc_dmag_K1[4,0],xerr=mcmc_dmag_K1[4,1],yerr=mcmc_dmag_K1[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_dmag_K1[5,0]-mcmc_dmag_K2[5,0],mcmc_dmag_K1[5,0],xerr=mcmc_dmag_K1[5,1],yerr=mcmc_dmag_K1[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_dmag_K1[6,0]-mcmc_dmag_K2[6,0],mcmc_dmag_K1[6,0],xerr=mcmc_dmag_K1[6,1],yerr=mcmc_dmag_K1[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_dmag_K1[7,0]-mcmc_dmag_K2[7,0],mcmc_dmag_K1[7,0],xerr=mcmc_dmag_K1[7,1],yerr=mcmc_dmag_K1[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_dmag_K1[8,0]-mcmc_dmag_K2[8,0],mcmc_dmag_K1[8,0],xerr=mcmc_dmag_K1[8,1],yerr=mcmc_dmag_K1[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_dmag_K1[9,0]-mcmc_dmag_K2[9,0],mcmc_dmag_K1[9,0],xerr=mcmc_dmag_K1[9,1],yerr=mcmc_dmag_K1[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_dmag_K1[10,0]-mcmc_dmag_K2[10,0],mcmc_dmag_K1[10,0],xerr=mcmc_dmag_K1[10,1],yerr=mcmc_dmag_K1[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_dmag_K1[11,0]-mcmc_dmag_K2[11,0],mcmc_dmag_K1[11,0],xerr=mcmc_dmag_K1[11,1],yerr=mcmc_dmag_K1[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_dmag_K1[12,0]-mcmc_dmag_K2[12,0],mcmc_dmag_K1[12,0],xerr=mcmc_dmag_K1[12,1],yerr=mcmc_dmag_K1[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_dmag_K1[13,0]-mcmc_dmag_K2[13,0],mcmc_dmag_K1[13,0],xerr=mcmc_dmag_K1[13,1],yerr=mcmc_dmag_K1[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_dmag_K1[14,0]-mcmc_dmag_K2[14,0],mcmc_dmag_K1[14,0],xerr=mcmc_dmag_K1[14,1],yerr=mcmc_dmag_K1[14,1],fmt='--o',c=colors[14],label= "S12")
    x = np.linspace(11,13.5,num=2)
    # plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

if dmag_julia_K1_K2 == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{julia}K1$", fontsize=20)
    plt.xlabel('$\Delta m_{julia}K1 - \Delta m_{julia}K2$', fontsize=20)
    plt.title("K1 vs. K1-K2")
    plt.errorbar(julia_dmag_K1[0,0]-julia_dmag_K2[0,0],julia_dmag_K1[0,0],xerr=julia_dmag_K1[0,1],yerr=julia_dmag_K1[0,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(julia_dmag_K1[1,0]-julia_dmag_K2[1,0],julia_dmag_K1[1,0],xerr=julia_dmag_K1[1,1],yerr=julia_dmag_K1[1,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(julia_dmag_K1[2,0]-julia_dmag_K2[2,0],julia_dmag_K1[2,0],xerr=julia_dmag_K1[2,1],yerr=julia_dmag_K1[2,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(julia_dmag_K1[3,0]-julia_dmag_K2[3,0],julia_dmag_K1[3,0],xerr=julia_dmag_K1[3,1],yerr=julia_dmag_K1[3,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(julia_dmag_K1[4,0]-julia_dmag_K2[4,0],julia_dmag_K1[4,0],xerr=julia_dmag_K1[4,1],yerr=julia_dmag_K1[4,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(julia_dmag_K1[5,0]-julia_dmag_K2[5,0],julia_dmag_K1[5,0],xerr=julia_dmag_K1[5,1],yerr=julia_dmag_K1[5,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(julia_dmag_K1[6,0]-julia_dmag_K2[6,0],julia_dmag_K1[6,0],xerr=julia_dmag_K1[6,1],yerr=julia_dmag_K1[6,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(julia_dmag_K1[7,0]-julia_dmag_K2[7,0],julia_dmag_K1[7,0],xerr=julia_dmag_K1[7,1],yerr=julia_dmag_K1[7,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(julia_dmag_K1[8,0]-julia_dmag_K2[8,0],julia_dmag_K1[8,0],xerr=julia_dmag_K1[8,1],yerr=julia_dmag_K1[8,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(julia_dmag_K1[9,0]-julia_dmag_K2[9,0],julia_dmag_K1[9,0],xerr=julia_dmag_K1[9,1],yerr=julia_dmag_K1[9,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(julia_dmag_K1[10,0]-julia_dmag_K2[10,0],julia_dmag_K1[10,0],xerr=julia_dmag_K1[10,1],yerr=julia_dmag_K1[10,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(julia_dmag_K1[11,0]-julia_dmag_K2[11,0],julia_dmag_K1[11,0],xerr=julia_dmag_K1[11,1],yerr=julia_dmag_K1[11,1],fmt='--o',c=colors[14],label= "S12")
    plt.errorbar(julia_dmag_K1[12,0]-julia_dmag_K2[12,0],julia_dmag_K1[12,0],xerr=julia_dmag_K1[12,1],yerr=julia_dmag_K1[12,1],fmt='--o',c=colors[15],label= "S13")
    plt.errorbar(julia_dmag_K1[13,0]-julia_dmag_K2[13,0],julia_dmag_K1[13,0],xerr=julia_dmag_K1[13,1],yerr=julia_dmag_K1[13,1],fmt='--o',c=colors[16],label= "S14")
    plt.errorbar(julia_dmag_K1[14,0]-julia_dmag_K2[14,0],julia_dmag_K1[14,0],xerr=julia_dmag_K1[14,1],yerr=julia_dmag_K1[14,1],fmt='--o',c=colors[17],label= "S15")
    plt.errorbar(julia_dmag_K1[15,0]-julia_dmag_K2[15,0],julia_dmag_K1[15,0],xerr=julia_dmag_K1[15,1],yerr=julia_dmag_K1[15,1],fmt='--o',c=colors[18],label= "S16")
    x = np.linspace(11,13.5,num=2)
    # plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with Julia Dmag K1
if dmag_simplex_K1_mcmc == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{simplex}$", fontsize=20)
    plt.xlabel('$\Delta m_{MCMC}$', fontsize=20)
    plt.title("$\Delta m$ - K1")
    # plt.errorbar(mcmc_dmag_K1[0,0],simplex_dmag_K1[0,0],xerr=mcmc_dmag_K1[0,1],yerr=simplex_dmag_K1[0,1],fmt='--o',c=colors[0],label= "Ad")
    # plt.errorbar(mcmc_dmag_K1[2,0],simplex_dmag_K1[2,0],xerr=mcmc_dmag_K1[2,1],yerr=simplex_dmag_K1[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_dmag_K1[3,0],simplex_dmag_K1[3,0],xerr=mcmc_dmag_K1[3,1],yerr=simplex_dmag_K1[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_dmag_K1[4,0],simplex_dmag_K1[4,0],xerr=mcmc_dmag_K1[4,1],yerr=simplex_dmag_K1[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_dmag_K1[5,0],simplex_dmag_K1[5,0],xerr=mcmc_dmag_K1[5,1],yerr=simplex_dmag_K1[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_dmag_K1[6,0],simplex_dmag_K1[6,0],xerr=mcmc_dmag_K1[6,1],yerr=simplex_dmag_K1[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_dmag_K1[7,0],simplex_dmag_K1[7,0],xerr=mcmc_dmag_K1[7,1],yerr=simplex_dmag_K1[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_dmag_K1[8,0],simplex_dmag_K1[8,0],xerr=mcmc_dmag_K1[8,1],yerr=simplex_dmag_K1[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_dmag_K1[9,0],simplex_dmag_K1[9,0],xerr=mcmc_dmag_K1[9,1],yerr=simplex_dmag_K1[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_dmag_K1[10,0],simplex_dmag_K1[10,0],xerr=mcmc_dmag_K1[10,1],yerr=simplex_dmag_K1[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_dmag_K1[11,0],simplex_dmag_K1[11,0],xerr=mcmc_dmag_K1[11,1],yerr=simplex_dmag_K1[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_dmag_K1[12,0],simplex_dmag_K1[12,0],xerr=mcmc_dmag_K1[12,1],yerr=simplex_dmag_K1[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_dmag_K1[13,0],simplex_dmag_K1[13,0],xerr=mcmc_dmag_K1[13,1],yerr=simplex_dmag_K1[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_dmag_K1[14,0],simplex_dmag_K1[14,0],xerr=mcmc_dmag_K1[14,1],yerr=simplex_dmag_K1[14,1],fmt='--o',c=colors[14],label= "S12")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()

## Comparison with simplex Dmag K2
if dmag_simplex_K2_mcmc == True:

    plt.grid(True)
    plt.legend()
    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)
    ax = plt.subplot(111)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("$\Delta m_{simplex}$", fontsize=20)
    plt.xlabel('$\Delta m_{MCMC}$', fontsize=20)
    plt.title("$\Delta m$ - K2")
    # plt.errorbar(mcmc_dmag_K2[0,0],simplex_dmag_K2[0,0],xerr=mcmc_dmag_K2[0,1],yerr=simplex_dmag_K2[0,1],fmt='--o',c=colors[0],label= "Ad")
    # plt.errorbar(mcmc_dmag_K2[2,0],simplex_dmag_K2[2,0],xerr=mcmc_dmag_K2[2,1],yerr=simplex_dmag_K2[2,1],fmt='--o',c=colors[2],label= "E")
    plt.errorbar(mcmc_dmag_K2[3,0],simplex_dmag_K2[3,0],xerr=mcmc_dmag_K2[3,1],yerr=simplex_dmag_K2[3,1],fmt='--o',c=colors[3],label= "S1")
    plt.errorbar(mcmc_dmag_K2[4,0],simplex_dmag_K2[4,0],xerr=mcmc_dmag_K2[4,1],yerr=simplex_dmag_K2[4,1],fmt='--o',c=colors[4],label= "S2")
    plt.errorbar(mcmc_dmag_K2[5,0],simplex_dmag_K2[5,0],xerr=mcmc_dmag_K2[5,1],yerr=simplex_dmag_K2[5,1],fmt='--o',c=colors[5],label= "S3")
    plt.errorbar(mcmc_dmag_K2[6,0],simplex_dmag_K2[6,0],xerr=mcmc_dmag_K2[6,1],yerr=simplex_dmag_K2[6,1],fmt='--o',c=colors[6],label= "S4")
    plt.errorbar(mcmc_dmag_K2[7,0],simplex_dmag_K2[7,0],xerr=mcmc_dmag_K2[7,1],yerr=simplex_dmag_K2[7,1],fmt='--o',c=colors[7],label= "S5")
    plt.errorbar(mcmc_dmag_K2[8,0],simplex_dmag_K2[8,0],xerr=mcmc_dmag_K2[8,1],yerr=simplex_dmag_K2[8,1],fmt='--o',c=colors[8],label= "S6")
    plt.errorbar(mcmc_dmag_K2[9,0],simplex_dmag_K2[9,0],xerr=mcmc_dmag_K2[9,1],yerr=simplex_dmag_K2[9,1],fmt='--o',c=colors[9],label= "S7")
    plt.errorbar(mcmc_dmag_K2[10,0],simplex_dmag_K2[10,0],xerr=mcmc_dmag_K2[10,1],yerr=simplex_dmag_K2[10,1],fmt='--o',c=colors[10],label= "S8")
    plt.errorbar(mcmc_dmag_K2[11,0],simplex_dmag_K2[11,0],xerr=mcmc_dmag_K2[11,1],yerr=simplex_dmag_K2[11,1],fmt='--o',c=colors[11],label= "S9")
    plt.errorbar(mcmc_dmag_K2[12,0],simplex_dmag_K2[12,0],xerr=mcmc_dmag_K2[12,1],yerr=simplex_dmag_K2[12,1],fmt='--o',c=colors[12],label= "S10")
    plt.errorbar(mcmc_dmag_K2[13,0],simplex_dmag_K2[13,0],xerr=mcmc_dmag_K2[13,1],yerr=simplex_dmag_K2[13,1],fmt='--o',c=colors[13],label= "S11")
    plt.errorbar(mcmc_dmag_K2[14,0],simplex_dmag_K2[14,0],xerr=mcmc_dmag_K2[14,1],yerr=simplex_dmag_K2[14,1],fmt='--o',c=colors[14],label= "S12")
    x = np.linspace(11,13.5,num=2)
    plt.plot(x,x,lw=2.8,c="black")
    plt.legend(prop={'size': 12})
    plt.show()
