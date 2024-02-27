#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 20:42:05 2023

@author: lesarmstrong
"""
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import openpyxl

import geopandas as gpd
from shapely.geometry import Point, Polygon
import geopandas as gpd
import pandas as pd
import numpy as np


# Area of a cell
A = 500 * 500
# Buffer of cavern to top and bottom of salt layer
buffer = 15
discount_rate = 0.055
cavern_lifetime = 30
working_gas = True

#df_all = pd.read_excel('/Users/lesarmstrong/Documents/GitHub/h2geostorage/cavern_optimization_scripts/salt_basin_data/ALL_salina_rheology_excluded_analysis.xlsx', sheet_name='Original_Data')
#df_all = pd.read_csv('/Users/lesarmstrong/Documents/GitHub/h2geostorage/cavern_optimization_scripts/salt_basin_data/MI_and_APP_salina_data_processed_python_excluded.csv')
df_all = pd.read_csv('/Users/lesarmstrong/Documents/GitHub/h2geostorage/cavern_optimization_scripts/salt_basin_data/MI_and_APP_salina_data_processed_python_excluded.csv')


# Create a GeoDataFrame


gdf_all = gpd.GeoDataFrame(df_all, geometry=[Point(xy) for xy in zip(df_all.LON, df_all.LAT)])


print('CREATED GEODATAFRAMES')


# Compute total bounds for each GeoDataFrame
minx1, miny1, maxx1, maxy1 = gdf_all .total_bounds

# Determine the overall minimum and maximum bounds
minx = minx1
miny = miny1
maxx = maxx1
maxy = maxy1

print(f"Overall bounds: {minx, miny, maxx, maxy}")


# for longitude, one degree is approximately 77.5 km in Michigan
dx = 500 / 77500
# 500 meters in degrees
dy = 500 / 111100  # for latitude, one degree is approximately 111.1 km


nx = int((maxx-minx)/dx)
ny = int((maxy-miny)/dy)

# Create polygon grid
grid_polys = []
for i in range(nx):
    for j in range(ny):
        minx2 = minx + dx * i
        miny2 = miny + dy * j
        maxx2 = minx + dx * (i + 1)
        maxy2 = miny + dy * (j + 1)
        grid_polys.append(Polygon([(minx2, miny2), (minx2, maxy2), (maxx2, maxy2), (maxx2, miny2)]))

print('POLYGRID CREATED')

# Create a geodataframe from the polygons.
#grid = gpd.GeoDataFrame(grid_polys, columns=['geometry'], crs=gdf_A2.crs)
grid = gpd.GeoDataFrame(grid_polys, columns=['geometry'])


def process_grid_data(grid=grid, gdf_all=gdf_all):
    def get_max_working_gas_and_row(row, gdfs):
        max_working_gas = -np.inf
        max_row = pd.Series(dtype='float')  # Create an empty Series to store the data
        
        for gdf in gdfs:
            intersecting_polygons = gdf[gdf.geometry.intersects(row.geometry)]
            if not intersecting_polygons.empty:
                current_max = intersecting_polygons['sphere_working_gas_cell'].max()
                if current_max > max_working_gas:
                    max_working_gas = current_max
                    max_row = intersecting_polygons.loc[intersecting_polygons['sphere_working_gas_cell'].idxmax()]
        return max_row

    max_rows_list = grid.apply(get_max_working_gas_and_row, gdfs=[gdf_all], axis=1)
    max_rows_df = pd.concat(max_rows_list, axis=1).transpose()

    grid['geometry'] = grid['geometry'].centroid
    grid['LON'] = grid['geometry'].x
    grid['LAT'] = grid['geometry'].y
    if 'index' in grid.columns:
        grid.drop('index', axis=1, inplace=True)
    grid = grid.replace([0, -np.inf], np.nan)
    grid.dropna(subset=['max_working_gas'], inplace=True)

    duplicates = gdf_all.duplicated(keep=False)
    num_duplicates = duplicates.sum()
    if num_duplicates > 0:
        duplicate_indices = gdf_all.index[duplicates]
        raise ValueError(f"There are {num_duplicates} duplicated 'sphere_working_gas_cell' values in the salt layer dataframe files at indices: {duplicate_indices.tolist()}")
    
    merged_all = pd.merge(grid[['LAT', 'LON', 'sphere_working_gas_cell']], gdf_all, on=['LAT', 'LON', 'sphere_working_gas_cell'], how='inner')
    return merged_all

# Example usage of the function:
# merged_result = process_grid_data(grid, gdf_all)
# merged_result.to_csv('clustering/APP_overlayed_salt_layers.csv')


# Example usage of the function:
# merged_result = process_grid_data(grid, gdf_all)
# merged_result.to_csv('clustering/APP_overlayed_salt_layers.csv')


grid = gpd.GeoDataFrame(grid_polys, columns=['geometry'])

def get_max_working_gas_and_row(row, gdfs):
    max_working_gas = -np.inf
    max_row = pd.Series(dtype='float')
    
    for gdf in gdfs:
        intersecting_polygons = gdf[gdf.geometry.intersects(row.geometry)]
        
        if not intersecting_polygons.empty:

            max_gas = intersecting_polygons['sphere_working_gas_cell'].max()
            #print(max_gas)
            if max_gas > max_working_gas:
                max_working_gas = max_gas
                #print(max_working_gas)
                # Capture the entire row, not just the max value
                max_row = intersecting_polygons.loc[intersecting_polygons['sphere_working_gas_cell'].idxmax()]
                print(max_row)
    return max_row



# Apply the function to each row in the grid DataFrame
# This will result in a Series for each row, which you can then concatenate into a DataFrame
max_rows_list = grid.apply(get_max_working_gas_and_row, gdfs=[gdf_all], axis=1)

# Original number of rows before dropping
original_row_count = max_rows_list.shape[0]

# Dropping rows where all elements are NaN or None
max_rows_list_dropped = max_rows_list.dropna(how='all')

# Number of rows after dropping
remaining_row_count = max_rows_list_dropped.shape[0]

# Number of rows dropped
dropped_row_count = original_row_count - remaining_row_count

print("Number of rows dropped:", dropped_row_count)
print("Number of rows remaining:", remaining_row_count)

# Concatenate all the Series into a DataFrame
#max_rows_df = pd.concat(max_rows_list, axis=1).transpose()

df_overlayed_excluded = max_rows_list_dropped.drop(columns=['Unnamed: 0', 'geometry'])

df_overlayed_excluded.to_csv('/Users/lesarmstrong/Documents/GitHub/h2geostorage/cavern_optimization_scripts/salt_basin_data/MI_and_APP_salina_data_processed_overlayed_and_pythonexcluded.csv')




# Save the result to a CSV file
#result.to_csv('clustering/APP_overlayed_salt_layers_max_working_gas.csv')
