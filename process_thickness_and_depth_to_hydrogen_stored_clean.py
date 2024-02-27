#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 16:27:41 2023

@author: lesarmstrong
"""

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog

import time
import os

from cavern_physical_model import sphere
from cavern_physical_model import horizontal_cylinder
from cavern_physical_model import vertical_cylinder

##############################################################################################################



# Example class representing objects
class Grid_Cell:
    def __init__(self, depth, thickness, LAT=0, LON=0, lat=0, lon=0):
        self.lat = lat # Coordinates CRS: EPSG:4326 CRS (also known as WGS84) # Degrees
        self.lon = lon # Coordinates CRS: EPSG:4326 CRS (also known as WGS84) # Degrees
        self.LAT = LAT # Coordinates CRS: EPSG:26989 CRS (also known as "NAD83 / Michigan Central") [meters]
        self.LON = LON # Coordinates CRS: EPSG:26989 CRS (also known as "NAD83 / Michigan Central") [meters]
        self.depth = depth # [m]
        self.thickness = thickness # [m]


_start_time = time.time()

def tic():
    global _start_time 
    _start_time = time.time()

def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec,60)
    (t_hour,t_min) = divmod(t_min,60) 
    print('Processing time: {}hour:{}min:{}sec'.format(t_hour,t_min,t_sec))

# Size of grid cell is 500m
a = 500
# Area of a cell
A = a * a
# Buffer of cavern to top and bottom of salt layer
buffer = 15
discount_rate = 0.055
cavern_lifetime = 30

# From Salina in table
n = 4.1  # Unitless
T_star = 8715  # Kelvin
A_ref = 2.7752 * (10**5)  # [MPa^n s]


def beddedsalts_capacity_calculation_from_csv(path='', shape='sphere',  a=500):
    """
    Calculate the capacity of bedded salts from a CSV file.

    This function reads a CSV file containing data about bedded salts and calculates their capacity.
    The CSV file should have specific columns (please specify these columns here).

    Parameters:
    filename (str): The path to the CSV file to read.

    Returns:
    float: The calculated capacity of the bedded salts.
    """

    # Track time to run script
    tic()
    # Grid Cell Area size (for number of caverns that fit in a cell later)
    A = a*a
    
    # Select .csv file with depth and thickness of bedded salt layer
    if path == '':    
        print('Please select .csv file with depth and thickness of bedded salt layer: ')
        print('')
            
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename()
        print("Path to selected salt layer is: ")
        print(path)
        print('')
    
    # Read in .csv file
    df = pd.read_csv(path)
    # Initialize list of Grid_Cell objects
    grid_cell_list = []

    # Loop through each row of .csv file
    for i in range(len(df)):
        
        # Create Grid_Cell object
        grid_cell = Grid_Cell(df['depth'][i], df['thickness'][i], LAT=df['LAT'][i], LON=df['LON'][i])
        
        # Buffer in salt that cavern needs [m]
        buffer = 15
        
        if shape == 'sphere':
            
            ## For Spheres ##
            radius = (df['thickness'][i] - buffer*2) / 2
            
            cavern = sphere(radius, df['depth'][i], n=n, T_star=T_star, A_ref=A_ref)
        
            Ncavernsthatfitinasinglecell_sphere = A / (np.pi * ((cavern.radius*4)**2))
            grid_cell.N_caverns_spheres = Ncavernsthatfitinasinglecell_sphere
            
        
            # Working gas
            masskg = cavern.working_gas_capacity #[kg]
            massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
            
            total_working_gas_sphere_capacity =  massmetric_ton * Ncavernsthatfitinasinglecell_sphere
            grid_cell.sphere_working_gas_cell = total_working_gas_sphere_capacity 
            
            
            # Cushion gas
            cushion_gas_masskg = cavern.cushion_gas
            cushion_gas_mass_metricton =cushion_gas_masskg / 1000
            
            # To get cushion gas of single cavern divide total cushion gas / N
            total_cushion_gas_sphere =  cushion_gas_mass_metricton * Ncavernsthatfitinasinglecell_sphere
            grid_cell.sphere_cushion_gas_cell = total_cushion_gas_sphere 
            
                            
            # Energy Density
            energy_sphere =  33.33 * total_working_gas_sphere_capacity * 1000 # metric ton --> kWh
            grid_cell.energy_density_sphere = energy_sphere / (cavern.volume * Ncavernsthatfitinasinglecell_sphere) # [kWh/m3]
        
        if shape == 'horizontal_cylinder':
            ## For horizontal cylindrical caverns ## 
            radius = (df['thickness'][i] - buffer*2) / 2
            lenght = a
            cavern = horizontal_cylinder(radius, df['depth'][i], lenght)
                                        
            Ncavernsthatfitinasinglecell_horizantalcylinder = A / (cavern.lenght * (cavern.radius*4))                             
            grid_cell.N_caverns_horizontal_cylinder = Ncavernsthatfitinasinglecell_horizantalcylinder
            
            # Working gas
            masskg = cavern.working_gas_capacity #[kg]
            massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
            
            total_working_gas_horizontal_cylinder_capacity =  massmetric_ton * Ncavernsthatfitinasinglecell_horizantalcylinder
            grid_cell.horizontal_cylinder_working_gas_cell = total_working_gas_horizontal_cylinder_capacity
            
            
            # Cushion gas
            cushion_gas_masskg = cavern.cushion_gas
            cushion_gas_mass_metricton =cushion_gas_masskg / 1000
            
            total_cushion_gas_horizontal_cylinder =  cushion_gas_mass_metricton * Ncavernsthatfitinasinglecell_horizantalcylinder
            grid_cell.horizontal_cylinder_working_gas_cell = total_cushion_gas_horizontal_cylinder

            
            # Energy Density
            energy_horizontal_cylinder =  33.33 * total_working_gas_horizontal_cylinder_capacity * 1000 # metric ton --> kWh
            grid_cell.energy_density_horizontal_cylinder = energy_horizontal_cylinder / (cavern.volume * Ncavernsthatfitinasinglecell_horizantalcylinder) # [kWh/m3]
            grid_cell_list.append(grid_cell)
        
        print('')    
        print('Object List Created')
        print('')
            
        # Create a list of attribute names
        attribute_names = list(vars(grid_cell_list[0]).keys())
        
        # Create a dictionary to store the attributes of the objects
        data = {attr: [getattr(obj, attr) for obj in grid_cell_list] for attr in attribute_names}
        
        # Create a DataFrame from the dictionary
        df = pd.DataFrame(data)
        
        # Display the resulting DataFrame
        print(df)
            
        processed_path = path[:-4] + '_processed.csv'
        df.to_csv(processed_path)
    
        print('')
        print('Script complete')
        tac()



def saltdome_capacity_calculation_from_csv(path, save_objects=False):
    """
    Calculate the capacity of salt domes from a CSV file.

    This function reads a CSV file containing data about salt domes and calculates their capacity.
    The CSV file should have specific columns (please specify these columns here).

    Parameters:
    path (str): The path to the CSV file to read.
    save_objects (bool, optional): If True, the function will save intermediate objects to disk. Defaults to False.

    Returns:
    float: The calculated capacity of the salt domes.

    """
    # Select .csv file with depth and thickness of bedded salt layer
    if path == '':    
        print('Please select .csv file with depth and thickness of bedded salt layer: ')
        print('')
            
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename()
        print("Path to selected salt layer is: ")
        print(path)
        print('')
        
    df = pd.read_csv(path)

    # List of working gas capacity of single caverns
    working_gas_sphere_single_cavern_list = []
    working_gas_sphere_list = []

    cushion_gas_sphere_single_cavern_list = []
    cushion_gas_sphere_list = []
    
    sphere_object_list = [] 
    
    working_gas_vertical_cylinder_single_cavern_list = []
    working_gas_vertical_cylinder_list = []
    
    cushion_gas_vertical_cylinder_single_cavern_list = []
    cushion_gas_vertical_cylinder_list = []    
    vertical_cylinder_object_list = []


    for i in range(len(df)):
        
        #Radius of workable salt dome
        dome_radius = df['radius'][i]
        
        # Buffer in salt that cavern needs [m] (ref: Michael Susan et al 2019)
        ## REVIEW BUFFER ## 
        buffer_between_caverns = 210
        buffer_from_walls = 150
        
        # Why max radius?? Pulled this value from Michael Susan thesis
        ## REVIEW MAX RADIUS ##
        radius_cavern = 40 #[m]
        
        
        
        # Creates array of possible depths of construction that are deeper than when the salt domes begins. 7000m is used as an exagerated number
        depth_range = np.array(df['depth'][i]) + 70000 
                 
        # Choose actual available depth that is closest to the optimal depth
        # For now optimal depth is the calculated 1500m. May optimze later.
        ## REVIEW OPTIMAL DEPTH ##
        
        optimal_depth = 1100                                                      
        idx = (np.abs(depth_range - optimal_depth)).argmin()   
        depth_cavern = depth_range[idx]  

        # Circle packaging (ref: https://www.quora.com/How-can-I-work-out-how-many-small-circles-I-can-fit-into-a-big-circle-See-additional-information)                
          
        # Number of caverns that fit in dome. Take away wall buffers from salt dome radius and add buffer between caverns to cavern radius to account for appropriate buffering      
        N_caverns_in_saltdome = (np.pi * ((dome_radius - buffer_from_walls)**2) ) / ( ((radius_cavern + buffer_between_caverns)**2) * np.sqrt(12))


        ## For vertical caverns ##
        
        # Assume height
       
        ## HOW DO WE KNOW THE HEIGHT OF VERTICAL CAVERN??  THIS IS A PLACEHOLDER BASED ON HEURISTICS
        ## REVIEW HEIGHT OF CAVERN ##
        height_cavern = radius_cavern * 4       
          
        
        cavern = vertical_cylinder(radius_cavern, depth_cavern, height_cavern)
           
        # Working gas
        masskg = cavern.working_gas_capacity #[kg]
        massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
        working_gas_sphere_single_cavern_list.append(massmetric_ton)
        
        total_working_gas_sphere_capacity =  massmetric_ton * N_caverns_in_saltdome
        working_gas_sphere_list.append(total_working_gas_sphere_capacity)
        
        # Cushion gas
        cushion_gas_masskg = cavern.cushion_gas
        cushion_gas_mass_metricton =cushion_gas_masskg / 1000
        cushion_gas_sphere_single_cavern_list.append(cushion_gas_mass_metricton)
        
        total_cushion_gas_sphere =  cushion_gas_mass_metricton * N_caverns_in_saltdome
        cushion_gas_sphere_list.append(total_cushion_gas_sphere)
        
        # Cavern object data
        if save_objects == True:
            sphere_object_list.append(cavern)
            
        ## For Spheres ##        
        
        cavern = sphere(radius_cavern, depth_cavern)
       
        # Working gas
        masskg = cavern.working_gas_capacity #[kg]
        massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
        working_gas_sphere_single_cavern_list.append(massmetric_ton)
        
        total_working_gas_sphere_capacity =  massmetric_ton * N_caverns_in_saltdome
        working_gas_sphere_list.append(total_working_gas_sphere_capacity)
        
        # Cushion gas
        cushion_gas_masskg = cavern.cushion_gas
        cushion_gas_mass_metricton =cushion_gas_masskg / 1000
        cushion_gas_sphere_single_cavern_list.append(cushion_gas_mass_metricton)
        
        total_cushion_gas_sphere =  cushion_gas_mass_metricton * N_caverns_in_saltdome
        cushion_gas_sphere_list.append(total_cushion_gas_sphere)
            
        # Cavern object data
        if save_objects == True:
            sphere_object_list.append(cavern)
        

        df['vertical_cylinder_working_gas_total']  = working_gas_vertical_cylinder_list 
        df['vertical_cylinder_working_gas_single_cavern']  = working_gas_vertical_cylinder_single_cavern_list  
        df['vertical_cylinder_number_of_caverns'] = df['vertical_cylinder_working_gas_total'] /  df['vertical_cylinder_working_gas_single_cavern']     
        df['vertical_cylinder_cushion_gas_total'] = cushion_gas_vertical_cylinder_list
        df['vertical_cylinder_cushion_gas_single_cavern'] = cushion_gas_vertical_cylinder_single_cavern_list

        processed_path = path[:-4] + '_processed.csv'
        
        df.to_csv('processed_data/' + processed_path)
        
        # Cavern object data
        if save_objects == True:
            df['cavern_sphere_object'] = sphere_object_list
            df['caern_vertical_cylinder_object'] = vertical_cylinder_object_list
            
            df.to_pickle('processed_data_objects_pickle/' + processed_path)



