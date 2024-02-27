#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:25:41 2023

@author: lesarmstrong
"""

import numpy as np
from eos import create_eos
import preos
import time
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP



# from cavern_cost_calculator import compressor_stage_number
# from cavern_cost_calculator import compressor_calculator
# from cavern_cost_calculator import cavern_cost_calculation


##############################################################################################################

_start_time = time.time()


def tic():
    global _start_time
    _start_time = time.time()


def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Processing time: {}hour:{}min:{}sec'.format(t_hour, t_min, t_sec))



def moles_to_pressure(n, T, V, fluid='hydrogen'):
    # Set the state using moles, volume, and temperature
    CP.set_reference_state(fluid, 'DEF')
    molar_mass = CP.PropsSI('M', fluid)  # kg/mol
    mass = n * molar_mass  # total mass in kg
    P = CP.PropsSI('P', 'T', T, 'D', mass/V, fluid) # Pa
    P = P / 1000000
    
    return P

def pressure_to_moles(p, V, T, fluid='hydrogen'):
    
    p = p * 1000000 # MPa --> Pa
    # Set the state using pressure, volume, and temperature
    CP.set_reference_state(fluid, 'DEF')
    molar_mass = CP.PropsSI('M', fluid)  # kg/mol
    density = CP.PropsSI('D', 'T', T, 'P', p, fluid)
    mass = density * V  # total mass in kg
    n = mass / molar_mass
    return n




# From sheet
# Currently using median values from table that Bradford prepared
#n = 4.73  # Unitless
#T_star = 6495  # Kelvin
#A_ref = 452.31  # [MPa^n s]

class Cavern:
    '''The `Cavern` class represents a cavern with attributes for depth, shape, height, latitude, longitude, and temperature.
    **Attributes:**
    - `depth` (float): The depth to top of the cavern.
    - `shape` (str): The shape of the cavern. Can be 'sphere' or 'cylinder'.
    - `height` (float): The height of the cavern.
    - `n` (float, optional): The convergence number. Defaults to 4.73.
    - `T_star` (float, optional): The reference temperature. Defaults to 6495.
    - `A_ref` (float, optional): The reference constant of the power flow law. Defaults to 452.31.
    - `lat` (str, optional): The latitude of the cavern. Defaults to an empty string.
    - `lon` (str, optional): The longitude of the cavern. Defaults to an empty string.
    - `temperature` (float): The temperature of the cavern.'''


    def __init__(self, depth, shape, height, n=4.73, T_star=6495, A_ref=452.31, lon='', lat=''):

       
        self.D = 6 * 10 ** (-7)  # MPa^-n year^-1

        # Convergence number. Has to do with salt properties
        self.n = n
        self.T_star = T_star
        self.A_ref = A_ref

        self.depth = depth
        self.height = height
        self.shape = shape

        self.lat = lat
        self.lon = lon

        self.temperature = 290 + 0.023 * (self.depth + (self.height) / 2)

    def __str__(self):
        return (f'r={self.radius}, d={self.depth} h={self.height}, shape={self.shape}')



# The legal maximum diameter for caverns in the Netherlands is 90 m


class sphere(Cavern):
    ''' The `sphere` class represents a spherical cavern and inherits from the `Cavern` class. It calculates various attributes related to the cavern based on the given parameters.

    **Attributes:**
    - `radius` (float): The radius of the cavern.
    - `depth` (float): The depth of the cavern.
    - `shape` (str, optional): The shape of the cavern. Defaults to 'sphere'.
    - `n` (float, optional): The convergence number. Defaults to 4.73.
    - `T_star` (float, optional): The reference temperature. Defaults to 6495.
    - `A_ref` (float, optional): The reference constant of the power flow law.
    - `lat` (str, optional): The latitude of the cavern. Defaults to an empty string.
    - `lon` (str, optional): The longitude of the cavern. Defaults to an empty string.'''

    def __init__(self, radius, depth, shape='sphere', n=4.73, T_star=6495, A_ref=452.31, lat='', lon=''):

        Cavern.__init__(self, depth, shape, height=radius * 2, n=n, T_star=T_star, A_ref=A_ref, lat=lat, lon=lon)

        self.radius = radius

        self.volume = 4 * (np.pi) * (self.radius ** 3) / 3

        # Varies between 2200 and 2300
        row_rock = 2300

        # Gas constant (m^3 Pa / K mol) or (J / K mol)
        R = 8.3144

        # Gravity (m/s^2)
        g = 9.81

        # Overburden pressures [MPa] / Ref: Caglayan
        self.P_over_top = (row_rock * g * self.depth) / 1000000 # [MPa]
        self.P_over_bot = (row_rock * g * (self.depth + self.radius * 2)) / 1000000 # [MPa]

        # Overburden pressure at the bottom is the geostatic pressure
        p_geostatic = self.P_over_bot

        # Max pressure [MPa] / Can change safety constant 0.8 to something else
        self.P_max = self.P_over_top * 0.8
        
        # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
        p_min1 = self.P_over_bot * 0.24
        self.p_min_stress = p_min1

        # Algebra eq (29) from Ma et al paper where p_internal is Pc. / Ref: Xuqiang Ma et al., “Creep Deformation Analysis of Gas Storage in Salt Caverns,” International Journal of Rock Mechanics and Mining Sciences 139 (March 1, 2021): 104635, https://doi.org/10.1016/j.ijrmms.2021.104635.
        # p_min2 = p_geostatic - ( np.log(1 - delta_v) / ( t * (-3/2) * self.D * (3/(2*self.n))**self.n )  )**(1/self.n)

        ## FROM BRADFORD's NOTES ##
        # For 30% collapse over 30 years, D' (the rate of collapse) is 0.012/year
        # This means that sigma_star = 1.n * [0.012/Dref * alpha^(n+1)]^n
        # Where Dref = Aref * exp(-T_star/Tref)
        # Where Aref, T_star, and n come from SaltRheology.xlsx sheet. Values calculated as median from values in table compiled by Bradford Hager

        # For sphere #
        alpha = 3 / 2

        # Calculated from above values
        T_ref = self.temperature
        D_ref = self.A_ref * np.exp(-self.T_star / T_ref)

        sigma_star = self.n * (0.012 / (alpha ** (self.n + 1) * D_ref)) ** (1 / self.n)

        # Critical pressure when simga_star = P_lith - P_cavern
        p_min2 = p_geostatic - sigma_star

        self.sigma_star = sigma_star

        self.p_min_convergence = p_min2

        # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
        if p_min2 > p_min1:

            self.P_min = p_min2

        else:
            self.P_min = p_min1

        # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
        # Does not allow for negative mass, instead sets to 0
        if self.P_min >= self.P_max:
            self.P_min = self.P_max

        # Calculates max and min possible mass given pressure contraints
        # Goes onto calculating working gas capacity and cushion gas
        self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000  # [kg]
        self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000  # [kg]

        # Max and min mass give us the paramount Working Gas Capacity
        self.working_gas_capacity = self.Mass_max - self.Mass_min  # [kg]

        # Cushion gas is the amount of hydrogen that always needs to be present in the well
        self.cushion_gas = self.Mass_min

        ## ASSUMING THAT THE INITIAL PRESSURE IS PMIN
        self.pressure = self.P_min

        # Assigning variables to cavern after assuming that the initial pressure is pmi
        self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
        self.mass = self.moles * 2.016 / 1000  # [Kg]
        self.energy = self.working_gas_capacity * 33.33 / 1000  # [MW/h]

        ##### Number of intervals of moles calculated between max and min pressure for plotting
        n_intervals = 10
        n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
        n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
        n_array = np.linspace(n_min, n_max, n_intervals)

        # Calculating max/min/current compressibility factor for H2 as a function of
        self.Zmin = (self.P_min * self.volume) / (n_min * R * self.temperature)  # PV/nRT
        self.Zmax = (self.P_max * self.volume) / (n_max * R * self.temperature)
        self.zcompressibility = (self.pressure * self.volume) / (self.moles * R * self.temperature)

        self.n_array = n_array

        pressures = np.array([])
        for i in n_array:
            pressures = np.append(moles_to_pressure(n=i, T=self.temperature, V=self.volume), pressures)
        pressures = np.flip(pressures)
        self.pressures = pressures

        self.m_array = n_array * 2.016 / 1000000  # [metric tons]

        self.e_array = self.m_array * 33.33 / 1000  # [MW/h]



class horizontal_cylinder(Cavern):
    
    def __init__(self, radius, depth, lenght, n=4.73, T_star=6495, A_ref=452.31, lat='', lon='', shape='horizontal_cylinder'):

        Cavern.__init__(self, depth, shape, height=radius * 2, n=n, T_star=T_star, A_ref=A_ref, lat=lat, lon=lon)

        self.radius = radius
        self.lenght = lenght

        T = self.temperature
        self.volume = (np.pi) * (self.radius ** 2) * self.lenght

        # Steady state creep (Lankof 2022)
        n = 4.089

        # Varies between 2200 and 2300
        row_rock = 2300


        # Gas constant (m^3 Pa / K mol) or (J / K mol)
        R = 8.3144

        # Gravity (m/s^2)
        g = 9.81

        # Overburden pressures [MPa] / Ref: Caglayan
        self.P_over_top = (row_rock * g * self.depth) / 1000000
        self.P_over_bot = (row_rock * g * (self.depth + self.radius * 2)) / 1000000

        # Max pressure [MPa] / Can change safety constant 0.8 to something else
        self.P_max = self.P_over_top * 0.8

        # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
        p_min1 = self.P_over_bot * 0.24
        self.p_min_stress = p_min1

        # pmin_2 is related to creep.

        # Overburden pressure at the bottom is the geostatic pressure
        p_geostatic = self.P_over_bot

        # Algebra eq (29) from Ma et al paper where p_internal is Pc. / Ref: Xuqiang Ma et al., “Creep Deformation Analysis of Gas Storage in Salt Caverns,” International Journal of Rock Mechanics and Mining Sciences 139 (March 1, 2021): 104635, https://doi.org/10.1016/j.ijrmms.2021.104635.
        # p_min2 = p_geostatic - ( np.log(1 - delta_v) / ( t * (-3/2) * self.D * (3/(2*self.n))**self.n )  )**(1/self.n)

        # For 30% collapse over 30 years, D' (the rate of collapse) is 0.012/year
        # This means that sigma_star = 1.n * [0.012/Dref * alpha^(n+1)]^n
        # Where Dref = Aref * exp(-T_star/Tref)
        # Where Aref, T_star, and n come from SaltRheology.xlsx sheet. Values calculated as median from values in table compiled by Bradford Hager

        # For Cylinder #
        alpha = np.sqrt(3)

        # Calculated from above values
        T_ref = self.temperature
        D_ref = self.A_ref * np.exp(-self.T_star / T_ref)

        sigma_star = n * (0.012 / (alpha ** (self.n + 1) * D_ref)) ** (1 / self.n)

        # Critical pressure when simga_star = P_lith - P_cavern
        p_min2 = p_geostatic - sigma_star

        self.p_min_convergence = p_min2

        # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
        if p_min2 > p_min1:
            self.P_min = p_min2

        else:
            self.P_min = p_min1

        # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
        if self.P_min >= self.P_max:
            self.P_min = self.P_max

        # Calculates max and min possible mass given pressure contraints
        # Goes onto calculating working gas capacity and cushion gas
        self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000  # [kg]
        self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000  # [kg]

        # Max and min mass give us the paramount Working Gas Capacity
        self.working_gas_capacity = self.Mass_max - self.Mass_min  # [kg]

        # Cushion gas is the amount of hydrogen that always needs to be present in the well
        self.cushion_gas = self.Mass_min

        ## ASSUMING THAT THE INITIAL PRESSURE IS PMIN
        self.pressure = self.P_min

        # Assigning variables to cavern after assuming that the initial pressure is pmi
        self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
        self.mass = self.moles * 2.016 / 1000  # [Kg]
        self.energy = self.working_gas_capacity * 33.33 / 1000  # [MW/h]

        # For potential plotting
        ##### Number of intervals of moles calculated between max and min pressure
        n_intervals = 10
        n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
        n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
        n_array = np.linspace(n_min, n_max, n_intervals)

        # Calculating max/min/current compressibility factor for H2 as a function of
        self.Zmin = (self.P_min * self.volume) / (n_min * R * self.temperature)  # PV/nRT
        self.Zmax = (self.P_max * self.volume) / (n_max * R * self.temperature)
        self.zcompressibility = (self.pressure * self.volume) / (self.moles * R * self.temperature)

        self.n_array = n_array

        self.n_array = n_array

        pressures = np.array([])
        for i in n_array:
            pressures = np.append(moles_to_pressure(n=i, T=self.temperature, V=self.volume), pressures)
        pressures = np.flip(pressures)
        self.pressures = pressures

        self.m_array = n_array * 2.016 / 1000000  # [metric tons]

        self.e_array = self.m_array * 33.33 / 1000  # [MW/h]





class vertical_cylinder(Cavern):

    def __init__(self, radius, depth, height, n=4.73, T_star=6495, A_ref=452.31, lat='', lon='', shape='vertical cylinder'):

        Cavern.__init__(self, depth, shape, height, n=n, T_star=T_star, A_ref=A_ref, lat=lat, lon=lon)

        self.radius = radius
        self.volume = (np.pi) * (self.radius ** 2) * self.height

        # Steady state creep (Lankof 2022)
        n = 4.089

        # Varies between 2200 and 2300
        row_rock = 2300


        # Gas constant (m^3 Pa / K mol) or (J / K mol)
        R = 8.3144

        # Gravity (m/s^2)
        g = 9.81

        # Overburden pressures [MPa] / Ref: Caglayan
        self.P_over_top = (row_rock * g * self.depth) / 1000000
        self.P_over_bot = (row_rock * g * (self.depth + self.height)) / 1000000

        # Max pressure [MPa] / Can change safety constant 0.8 to something else
        self.P_max = self.P_over_top * 0.8

        ## ESPECIALLY MIN PRESSURE
        # Min pressure [MPa]. Two kinds of pressures that compete with each other: Stress & Creep. We take the largest value of the two as the base of the working gas pressure.

        # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
        p_min1 = self.P_over_bot * 0.24
        self.p_min_stress = p_min1

        # pmin_2 is related to creep.

        # Overburden pressure at the bottom is the geostatic pressure
        p_geostatic = self.P_over_bot
        
        # For Cylinder #
        alpha = np.sqrt(3)

        # Calculated from above values
        T_ref = self.temperature
        D_ref = A_ref * np.exp(-T_star / T_ref)

        sigma_star = n * (0.012 / (alpha ** (n + 1) * D_ref)) ** (1 / n)

        # Critical pressure when simga_star = P_lith - P_cavern
        p_min2 = p_geostatic - sigma_star

        self.p_min_convergence = p_min2

        # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
        if p_min2 > p_min1:

            self.P_min = p_min2

        else:
            self.P_min = p_min1

        # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
        if self.P_min >= self.P_max:
            self.P_min = self.P_max

        # Calculates max and min possible mass given pressure contraints
        # Goes onto calculating working gas capacity and cushion gas
        self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000  # [kg]
        self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000  # [kg]

        # Max and min mass give us the paramount Working Gas Capacity
        self.working_gas_capacity = self.Mass_max - self.Mass_min  # [kg]

        ### ????? ###
        # self.cushion_gas =  self.working_gas_capacity / self.Mass_min

        # Cushion gas is the amount of hydrogen that always needs to be present in the well
        self.cushion_gas = self.Mass_min
        # Assuming the initial pressure is p_min
        self.pressure = self.P_min

        # Assigning variables to cavern after assuming that the initial pressure is pmi
        self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
        self.mass = self.moles * 2.016 / 1000  # [Kg]
        self.energy = self.working_gas_capacity * 33.33 / 1000  # [MW/h]


        ##### Number of intervals of moles calculated between max and min pressure
        n_intervals = 10
        n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
        n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
        n_array = np.linspace(n_min, n_max, n_intervals)

        # Calculating max/min/current compressibility factor for H2 as a function of
        self.Zmin = (self.P_min * self.volume) / (n_min * R * self.temperature)  # PV/nRT
        self.Zmax = (self.P_max * self.volume) / (n_max * R * self.temperature)
        self.zcompressibility = (self.pressure * self.volume) / (self.moles * R * self.temperature)

        self.n_array = n_array

        self.n_array = n_array

        pressures = np.array([])
        for i in n_array:
            pressures = np.append(moles_to_pressure(n=i, T=self.temperature, V=self.volume), pressures)
        pressures = np.flip(pressures)
        self.pressures = pressures

        self.m_array = n_array * 2.016 / 1000000  # [metric tons]

        self.e_array = self.m_array * 33.33 / 1000  # [MW/h]



# From Salina in table
n = 4.1  # Unitless
T_star = 8715  # Kelvin
A_ref = 2.7752 * (10**5)  # [MPa^n s]

# From median in table
n = 4.73  # Unitless
T_star = 6495  # Kelvin
A_ref = 425.31  # [MPa^n s]


def plot_delta_pressure_chart(n=4.1, T_star=8715, A_ref=2.7752 * (10**5)):
    '''
    Plots a chart showing the pressure ranges for a spherical cavern of radius 40 meters at different depths.

    Args:
    - n (float, optional): The convergence number. Defaults to 4.1.
    - T_star (float, optional): The reference temperature. Defaults to 8715.
    - A_ref (float, optional): The reference constant of the power flow law. Defaults to 2.7752 * (10**5).

    Returns:
    None

    Examples:
    plot_delta_pressure_chart()  # Plots the pressure chart using default parameters.
    plot_delta_pressure_chart(n=4.5, T_star=8000)  # Plots the pressure chart with custom parameters.
    '''

    radius = 40
            
    # Create an array of depth values from 150m to 3000m
    depth_values = np.arange(150, 2501, 1)
    
    # Initialize lists to store attribute values
    P_min_values = []
    P_max_values = []
    p_min_convergence_values = []
    p_min_stress_values = []
    
    # Calculate attributes for each depth value
    for depth in depth_values:
        cavern = sphere(radius=radius, depth=depth, n=n, T_star=T_star, A_ref=A_ref)
        
        P_min_values.append(cavern.P_min)
        P_max_values.append(cavern.P_max)
        p_min_convergence_values.append(cavern.p_min_convergence)
        p_min_stress_values.append(cavern.p_min_stress)
    
    # Find the depth where p_min_convergence overtakes p_min_stress
    convergence_overtake_depth = None
    for i, (p_convergence, p_stress) in enumerate(zip(p_min_convergence_values, p_min_stress_values)):
        if p_convergence > p_stress:
            convergence_overtake_depth = depth_values[i]
            print(convergence_overtake_depth)
            break
    
    # Set a higher DPI for better plot quality
    plt.figure(figsize=(8, 5), dpi=1200)

    
    # Dotted lines for p_min_convergence and p_min_stress with custom colors
    plt.plot(depth_values, p_min_convergence_values, label='Minimum Pressure by Convergence', linestyle='dashed', color='purple')
    plt.plot(depth_values, p_min_stress_values, label='Minimum Pressure by Stress', linestyle='dashed', color='green')
    
    # Solid lines for P_min and P_max using default Matplotlib colors
    plt.plot(depth_values, P_min_values, label='Minimum Pressure', linestyle='solid')
    plt.plot(depth_values, P_max_values, label='Maximum Pressure', linestyle='solid')
    
    # Red dotted line where p_min_convergence overtakes p_min_stress
    if convergence_overtake_depth is not None:
        plt.plot([convergence_overtake_depth, convergence_overtake_depth], [0, max(p_min_convergence_values)], color='red', linestyle='dashed', label='Depth Sweet-Spot')
    
    plt.ylim(0, None) 
    
    # Add labels and legend
    plt.xlabel('Depth to Top of Cavern [m]')
    plt.ylabel('Pressure [MPa]')
    
    # Reorder the legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2, 3, 0, 1, 4]  # Reorder the legend entries as desired
    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
    
    # Remove the grid
    plt.grid(False)
    
    # Title for the plot
    plt.title('Pressure ranges for a spherical cavern of radius 40 meters')

    
    # Save the plot as a high-quality PNG image
    plt.savefig('pressure_chart.png', dpi=300, bbox_inches='tight')
    
    # Display the plot
    plt.show()

# Call the function to generate the plot and save it as a high-quality image
plot_delta_pressure_chart(n=n, T_star=T_star, A_ref=A_ref)








