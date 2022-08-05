#Created by Sarah Weischer
#20220805

########################################################################################################################
# Imports
import pandas as pd
from picasso import io
import seaborn as sns
import matplotlib.pyplot as plt
import os.path as _ospath
import numpy as np
from picasso import render
import matplotlib.pyplot as plt
from picasso import postprocess

# Settings
sns.set_palette("Set1", 8, .75)
#sns.set_theme(style="whitegrid")
sns.set_style("ticks")
%matplotlib inline

# Filter
def mean_frame_filter(df):
    min_mean_frame = df.frame.mean() - df.frame.std()
    max_mean_frame = df.frame.mean() + df.frame.std()
    print("Mean frame:", df.frame.mean().round(), "Lower bound:", min_mean_frame.round(), "Upper bound:",
          max_mean_frame.round())
    df_filter = df.groupby(df.group).filter(lambda x: x["frame"].mean() > min_mean_frame)
    df_filter = df_filter.groupby(df_filter.group).filter(lambda x: x["frame"].mean() > min_mean_frame)
    print("Locs before filtering:", df.shape[0], "; Locs after filtering:", df_filter.shape[0])

    return df_filter


def standard_dev_filter(df, min_std=2500, max_std=12500):
    df_filter = df.groupby(df.group).filter(lambda x: x["frame"].std() > min_std)
    df_filter = df_filter.groupby(df.group).filter(lambda x: x["frame"].std() < max_std)

    print("Locs before filtering:", df.shape[0], "; Locs after filtering:", df_filter.shape[0])

    return df_filter


######################################################################################################################
# Paramters/User input
# To Do: Make config file ?

# Variables
pixel_size = 160 # in nm
em = 580 # emission wavelength in nm
author_date = "SW_20220805"

# Filter params
sx_center = 1.0
sy_center = 1.0
radius = 0.2

loc_prec_max = 0.2

# not yet implemented:
#photon_min = 500
#photon_max = 7000
#bg_max = 300

# DBSCAN params (radius, min density points)
DBSCAN_1_params = (0.08, 10)
DBSCAN_2_params = (0.08, 10)
DBSCAN_3_params = (0.05, 20)

loc_path = r"D:\PROJECTS\MiN_Data\Workgroups\Sarah\Project_DNA-PAINT\Talin\20220727\Processed\Cell1\Filtering_Parameters\Cell1_1_Zoom.hdf5"

########################################################################################################################
# MAIN

# Load localizations
locs, info = io.load_locs(loc_path)
print('Loaded {} locs.'.format(len(locs)))

# Calculate median localization precision and NeNA
med_loc_prec = np.median(postprocess.localization_precision(locs.photons, locs.sx, locs.bg, em))
NeNA = postprocess.nena(locs, info)[1]
print("Median localization precision:", med_loc_prec, "pixel, " , med_loc_prec * pixelsize,"nm.")
print("NeNA precision:", NeNA, "pixel, " , NeNA * pixelsize,"nm.")

# Filter on sx and sy
df_locs = pd.DataFrame.from_records(locs)
to_keep = (df_locs.sx-sx_center)**2 + (df_locs.sy-sy_center)**2 < radius**2
filtered_locs = df_locs[to_keep]
print('Length of locs before filtering {}, after filtering {}.'.format(len(df_locs),len(filtered_locs)))

# Filter on localization precision
len_before = len(filtered_locs)
to_keep = filtered_locs.lpx < loc_prec_max
filtered_locs = filtered_locs[to_keep]
to_keep = filtered_locs.lpy < loc_prec_max
filtered_locs = filtered_locs[to_keep]
print('Length of locs before filtering {}, after filtering {}.'.format(len_before ,len(filtered_locs)))

med_loc_prec = np.median(postprocess.localization_precision(filtered_locs.photons, filtered_locs.sx, filtered_locs.bg, em))
NeNA = postprocess.nena(filtered_locs, info)[1]
print("Median localization precision after filtering:", med_loc_prec * pixelsize,"nm.")
print("NeNA precision after filtering:", NeNA * pixelsize,"nm.")

# Save

# Perform first DBSCAN
clusters_1, cluster_locs_1 = postprocess.dbscan(filtered_locs, radius=DBSCAN_1_params[0], min_density=DBSCAN_1_params[1])
df_clustered_1 = pd.DataFrame.from_records(cluster_locs_1)
print("Number of clusters:", len(df_clustered_1.group.unique()))
print("Average number of localizations per cluster:", (df_clustered_1.frame.groupby(df_clustered_1.group).count()).mean().round())

# Save

# MFF
df_clustered_MFF = mean_frame_filter(df_clustered_1)
print("Number of clusters after MFF:", len(df_clustered_MFF.group.unique()))
print("Average number of localizations per cluster:", (df_clustered_MFF.frame.groupby(df_clustered_MFF.group).count()).mean().round())

# Perform second DBSCAN
clusters_2, cluster_locs_2 = postprocess.dbscan(df_clustered_MFF.to_records(index=False), radius=DBSCAN_2_params[0], min_density=DBSCAN_2_params[1])
df_clustered_2 = pd.DataFrame.from_records(cluster_locs_2)

# StDev Filter
df_clustered_stdev = standard_dev_filter(df_clustered_2, min_std=5000, max_std=50000)
print("Numer of clusters after StDev filtering:", len(df_clustered_stdev.group.unique()))
print("Average number of localizations per cluster:", (df_clustered_stdev.frame.groupby(df_clustered_stdev.group).count()).mean().round())

# Create a new dictionary for the new info
clustered_stdev = df_clustered_stdev.to_records(index=False)
new_info = {}
new_info["Generated by"] = "SW_20220805"
new_info["MFF"] = 'True'
new_info["StDev_filter"] = 'True'

info.append(new_info)

base, ext = _ospath.splitext(save_path)

new_path = base+'_DBSCAN_FILTER.hdf5'


io.save_locs(new_path, clustered_stdev, info)

print('{} locs saved to {}.'.format(len(clustered_stdev), new_path))