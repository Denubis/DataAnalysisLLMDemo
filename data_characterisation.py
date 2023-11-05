#!/usr/bin/env python3
# Import necessary libraries
import pandas as pd


import pyproj
import math

# import geopandas as gpd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import os
import matplotlib.dates as mdates
import glob

# from shapely.geometry import Point
from netCDF4 import MFDataset, num2date
from tqdm import tqdm
import sys
import pickle

# from rocksdict import Rdict
import sqlean as sqlite3



import pyproj
import math
#https://www.perplexity.ai/search/for-spatialite-how-MhMo3ufLTrOR47pmoHYrMw?s=c
def calculate_area(lat, lon, grid_size=1):
    # Define the WGS84 ellipsoid
    wgs84 = pyproj.Geod(ellps='WGS84')

    # Calculate the corner coordinates of the grid square
    lat1 = lat - grid_size / 2
    lat2 = lat + grid_size / 2
    lon1 = lon - grid_size / 2
    lon2 = lon + grid_size / 2

    # Calculate the lengths of the sides of the grid square
    _, _, side_ns = wgs84.inv(lon, lat1, lon, lat2)
    _, _, side_ew = wgs84.inv(lon1, lat, lon2, lat)

    # Calculate the area of the grid square
    area = math.fabs(side_ns * side_ew)

    return area



# Global variables or configurations (like paths to datasets, etc.)

# ---------------------------------------
# 1.1 Sea Ice Extent Daily (CSV)
# ---------------------------------------


def load_sea_ice_extent_csv(filepath):
    """
    Load the Sea Ice Extent CSV file and return a DataFrame.
    file is data/S_seaice_extent_daily_v3.0.csv
    The first line is headers, the second line is units.
    Skip the second line.
    """

    df = pd.read_csv(filepath, skiprows=[1])
    # trim spaces from headers
    df.rename(columns=lambda x: x.strip(), inplace=True)
    print(df.keys())
    # concatenate year month day into a single date column
    df = df.assign(Date=lambda x: pd.to_datetime(x[["Year", "Month", "Day"]]))
    print(df)
    return df


def plot_sea_ice_extent_timeseries(dataframe, outfile):
    """
    Plot a time series of the sea ice extent using seaborn.
    """

    # We want only Date, Missing, and Extent in this dataframe
    monthly_df = dataframe[["Date", "Missing", "Extent"]]
    monthly_avg = monthly_df.resample("M", on="Date").mean()

    plt.figure(figsize=(12, 6))

    # add a seaborn time series plot to the fig
    plot = sns.lineplot(data=monthly_avg, x="Date", y="Extent")
    plot.set(xlabel="Date", ylabel="Extent (10^6 sq km)")
    plot.set_title("Monthly Average Sea Ice Extent")

    # rotate dates on x-axis
    plot.set_xticklabels(
        plot.get_xticklabels(), rotation=45, horizontalalignment="right"
    )
    plot.grid(True)

    # exclude the first and last years for the yearly max and min
    # The monthly_df however, has multiple months in 1978 and 2023

    # calculate how many rows to skip. Monthly_df only has Date, not year. We want to see if the date object falls within 1978 or 2023
    # The monthly_df is still in days, so we need to find the Year from Date
    # Then we need to find how many rows are in that year to skip.
    offset_1978 = len(monthly_df[monthly_df["Date"].dt.year == 1978])
    offset_2023 = len(monthly_df[monthly_df["Date"].dt.year == 2023])

    print(offset_1978)
    yearly_max = (
        monthly_df[offset_1978 : offset_2023 * -1].resample("Y", on="Date").max()
    )
    yearly_min = (
        monthly_df[offset_1978 : offset_2023 * -1].resample("Y", on="Date").min()
    )

    plot = sns.lineplot(data=yearly_max, x="Date", y="Extent", color="red")
    plot = sns.lineplot(data=yearly_min, x="Date", y="Extent", color="green")
    # add a legend
    plot.legend(labels=["Monthly Average", "Yearly Maximum", "Yearly Minimum"])

    # Adjust x-axis ticks for better readability
    years = mdates.YearLocator()  # every year
    years_fmt = mdates.DateFormatter("%Y")
    plot.xaxis.set_major_locator(years)
    plot.xaxis.set_major_formatter(years_fmt)

    plt.tight_layout()

    plot.figure.savefig(outfile)


# ---------------------------------------
# 1.2 Gridded Argo Data (netCDF)
# ---------------------------------------


def load_argo_data(directory_path, sqlitedb="data/argo_data.db"):
    """
    Load the Argo netCDF file and return relevant data arrays.
    """

    # connect to the sqlitedb and enable spatialite extension
    with sqlite3.connect(sqlitedb, isolation_level=None) as db:
        db.enable_load_extension(True)
        db.load_extension("mod_spatialite")
        # Note to user, you'll need to install mod-spatalite here

        # initalise the db if necessary to hold the argo data.
        db.execute("create table surface_area(longitude real, latitude real, area real);")
        db.execute(
            """CREATE TABLE IF NOT EXISTS argo_data (
    source TEXT, 
    ARGO_TEMPERATURE_ANOMALY REAL, 
    ARGO_SALINITY_ANOMALY REAL, 
    LONGITUDE REAL, 
    LATITUDE REAL, 
    PRESSURE REAL, 
    TIME TEXT,
    PRIMARY KEY (LONGITUDE, LATITUDE, PRESSURE, TIME));
"""
        )
        # we're going to make a month by month average instead, so let's make an argo_monthly_data which captures average temp and sal for a given time.
        db.execute("""CREATE TABLE IF NOT EXISTS argo_monthly_data (
    source TEXT,
    MEAN_ARGO_TEMPERATURE_ANOMALY REAL,
    MEAN_ARGO_SALINITY_ANOMALY REAL,
    TIME TEXT PRIMARY KEY);""")



        # make the geometry column if it does not exist and only init the spatial metadata if it does not exist
        db.execute("SELECT InitSpatialMetaData(1);")
        # if (
        #     db.execute(
        #         "SELECT COUNT(*) FROM geometry_columns WHERE f_table_name = 'argo_data';"
        #     ).fetchone()[0]
        #     == 0
        # ):
        #     # We'll also need to intialise the spatial functions
        #     db.execute(
        #         "SELECT AddGeometryColumn('argo_data', 'geometry', 4326, 'POINT', 'XY');"
        #     )
        #     db.execute("SELECT CreateSpatialIndex('argo_data', 'geometry');")
        #     db.commit()s

    # Experiment with MFDataset
    # with MFDataset(f"{directory_path}/RG_ArgoClim_2*_2019.nc") as rootgrp:
    #     pass
    # Fails, no aggregation dimension. We'll do this the hard way.

    # let's define the shape of a geodataframe that can hold the data.
    # We will be loading the 2004-2018 data, and then the 2019-2023 data into a singular geodataframe
    # but we need to indicate the data source for each row.

    # To load the data into a dataframe, we first need to make a dictionary. There should be a Geometry column with a Point formed from the Latitude Longitude columns in the import
    # We will need to join the Salinity and Temperature 2004-2018 datasets together, to match the shape of the 2019-2023 data
    # The dictionary should also have the datasource, being either 'original' or 'extension'

    # data_2004 = Rdict("/tmp/data_2004.dd")

    for file in tqdm(
        glob.glob(f"{directory_path}/RG_ArgoClim_2*_2019.nc"),
        desc=f"Loading 2019-current data",
    ):
        with Dataset(file) as rootgrp:
            # These files have both temperature and salinity, so we can extract both.
            # print(rootgrp.variables)
            argo_temp = rootgrp["ARGO_TEMPERATURE_ANOMALY"][:]
            argo_sal = rootgrp["ARGO_SALINITY_ANOMALY"][:]
            longitude = rootgrp["LONGITUDE"][:]
            latitude = rootgrp["LATITUDE"][:]
            pressure = rootgrp["PRESSURE"][:]
            time = rootgrp["TIME"][:]

            for time_idx in tqdm(range(len(time)), desc="Time"):
                data = []
                # cur = db.cursor()
                # cur.execute("BEGIN TRANSACTION")
                for longitude_idx in tqdm(
                    range(len(longitude)), desc=f"Loading longitudes", leave=False
                ):
                    longitude_value = longitude[longitude_idx].item()

                    for latitude_idx in range(len(latitude)):
                        latitude_value = latitude[latitude_idx].item()

                        # if latitude_value > -40:
                        #     # print("skipping latitude", latitude_value)
                        #     continue
                        for pressure_idx in range(len(pressure)):
                            source_value = "extension"
                            # We need to cast these numpy objects to things that sqlite3 will accept
                            time_num = num2date(
                                times=time[time_idx],
                                units="months since 2004-01-01 00:00:00",
                                calendar="360_day",
                            )
                            time_value = time_num.strftime("%Y-%m-%d %H:%M:%S")

                            pressure_value = pressure[pressure_idx].item()

                            temperature_value = argo_temp[
                                time_idx, pressure_idx, latitude_idx, longitude_idx
                            ].item()

                            sal_value = argo_sal[
                                time_idx, pressure_idx, latitude_idx, longitude_idx
                            ].item()
                            # use spatialite or  geopandas to calculate the area of the lat/long cell for future averaging
                            # Spatialite query: 
                            # insert into surface_area (longitude, latitude, area) SELECT distinct longitude, latitude, area(st_expand(MakePoint(longitude, latitude, 4326),0.5),true) from argo_data;
                            # print(f"SELECT area(st_expand(MakePoint({longitude_value}, {latitude_value}, 4326),0.5),true);")
                            # area = db.execute(f"SELECT area(st_expand(MakePoint({longitude_value}, {latitude_value}, 4326),0.5),true)").fetchone()[0]

                            area = calculate_area(lat=latitude_value, lon=longitude_value, grid_size=0.5)

                            new_row = {
                                "source":source_value,
                                "temp":temperature_value,
                                "sal":sal_value,
                                "time":time_value,
                                "area_sqm":area
                            }
                            data.append(new_row)
                # make a dataframe from the data, then calculate the mean temp and sal based on a weighting from area_sqm
                df = pd.DataFrame(data)
                df["temp_weighted"] = df["temp"] * df["area_sqm"]
                df["sal_weighted"] = df["sal"] * df["area_sqm"]
                mean_temp = df["temp_weighted"].sum() / df["area_sqm"].sum()
                mean_sal = df["sal_weighted"].sum() / df["area_sqm"].sum()
                # insert into the argo_monthly_data table
                db.execute("INSERT INTO argo_monthly_data (source, MEAN_ARGO_TEMPERATURE_ANOMALY, MEAN_ARGO_SALINITY_ANOMALY, TIME) VALUES (?, ?, ?, ?);", (source_value, mean_temp, mean_sal, time_value)) 
                db.commit()

                # cur.executemany(
                #     "INSERT INTO argo_data (source, ARGO_TEMPERATURE_ANOMALY, ARGO_SALINITY_ANOMALY, LONGITUDE, LATITUDE, PRESSURE, TIME, geometry) VALUES (?, ?, ?, ?, ?, ?, ?, MakePoint(?, ?, 4326));",
                #     # cast the timevalue as a sqlite compatible date string
                #     data,
                # )
                #             print(
                #                 db.execute(
                #                     "SELECT source, ARGO_TEMPERATURE_ANOMALY, ARGO_SALINITY_ANOMALY, LONGITUDE, LATITUDE, PRESSURE, TIME, AsWKT(geometry) FROM argo_data WHERE LONGITUDE = ? AND LATITUDE = ? AND PRESSURE = ? AND TIME = ?;",
                #                     (longitude_value, latitude_value, pressure_value, time_value),
                #                 ).fetchone()
                #             )
                # print(time[time_idx], time_num, time_value)
                # print(
                #     source_value,
                #     temperature_value,
                #     sal_value,
                #     longitude_value,
                #     latitude_value,
                #     pressure_value,
                #     time_value,
                #     longitude_value,
                #     latitude_value,
                # )
                # let's QA this by querying what we just entered, making sure to use an aswkt to get the geometry back.

                # cur.execute("COMMIT")
            #     break
            # break

            # sys.exit(0)
            # data_2004[
            #     f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
            # ] = new_row
            # print(new_row)
            # sys.exit(0)
    # sys.exit(0)
    # for filetype in ["Temperature", "Salinity"]:
    with Dataset(
        f"{directory_path}/RG_ArgoClim_Temperature_2019.nc",
        "r",
        format="NETCDF3_CLASSIC",
    ) as nc_temp_file:
        with Dataset(
            f"{directory_path}/RG_ArgoClim_Salinity_2019.nc",
            "r",
            format="NETCDF3_CLASSIC",
        ) as nc_sal_file:
            # print(nc_file.variables)
            # We need to iterate over all four dimensions to add the data to the dict
            #  {'LONGITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LONGITUDE', size = 360, 'LATITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LATITUDE', size = 145, 'PRESSURE': <class 'netCDF4._netCDF4.Dimension'>: name = 'PRESSURE', size = 58, 'TIME': <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'TIME', size = 180}
            # if filetype == "Temperature":
            argo_temp = nc_temp_file.variables["ARGO_TEMPERATURE_ANOMALY"][:]
            # else:
            argo_sal = nc_sal_file.variables["ARGO_SALINITY_ANOMALY"][:]
            # argo_temp.shape: (180, 58, 145, 360)
            longitude = nc_temp_file.variables["LONGITUDE"][:]
            # print(longitude.shape)
            # (360,)
            latitude = nc_temp_file.variables["LATITUDE"][:]
            # print(latitude.shape)
            # (145,)
            pressure = nc_temp_file.variables["PRESSURE"][:]
            # print(pressure.shape)
            # (58,)
            time = nc_temp_file["TIME"][:]
            # print(time.shape)
            # (180,)

            # iterate over all of these to populate a dictionary so that we may join two datasources together
            for time_idx in tqdm(range(len(time))):
                data = []
                # cur = db.cursor()
                # cur.execute("BEGIN TRANSACTION")
                for longitude_idx in tqdm(
                    range(len(longitude)), desc=f"Loading {filetype} data", leave=False
                ):
                    longitude_value = longitude[longitude_idx].item()
                # continue if this longitude already exists in the db, with salinity and source original
                # if (
                #     db.execute(
                #         "SELECT COUNT(ARGO_SALINITY_ANOMALY) FROM argo_data WHERE LONGITUDE = ? and source= ?;",
                #         (longitude_value, "original"),
                #     ).fetchone()[0]
                #     > 0
                # ):
                #     continue

                    for latitude_idx in tqdm(range(len(latitude)), leave=False):
                        latitude_value = latitude[latitude_idx].item()

                    # if latitude_value > -40:
                    #     # print("skipping latitude", latitude_value)
                    #     continue
                        for pressure_idx in range(len(pressure)):
                            # print(longitude_idx, latitude_idx, pressure_idx, time_idx)

                            source_value = "original"
                            time_orig = num2date(
                                times=time[time_idx],
                                units="months since 2004-01-01 00:00:00",
                                calendar="360_day",
                            )
                            time_value = time_orig.strftime("%Y-%m-%d %H:%M:%S")
                            # LEt's run the same prep for mean calculation as we did above
                            # use spatialite to calculate the area of the lat/long cell for future averaging

                            # area = db.execute(f"SELECT area(st_expand(MakePoint({longitude_value}, {latitude_value}, 4326),0.5),true)").fetchone()[0]
                            area = calculate_area(lat=latitude_value, lon=longitude_value, grid_size=0.5)

                            temperature_value = argo_temp[
                                time_idx, pressure_idx, latitude_idx, longitude_idx].item()
                            salinity_value = argo_sal[
                                time_idx, pressure_idx, latitude_idx, longitude_idx].item()
                            new_row = {
                                "source":source_value,
                                "temp":temperature_value,
                                "sal":salinity_value,
                                "time":time_value,
                                "area_sqm":area
                            }
                            data.append(new_row)
                # make a dataframe from the data, then calculate the mean temp and sal based on a weighting from area_sqm
                df = pd.DataFrame(data)
                df["temp_weighted"] = df["temp"] * df["area_sqm"]
                df["sal_weighted"] = df["sal"] * df["area_sqm"]
                mean_temp = df["temp_weighted"].sum() / df["area_sqm"].sum()
                mean_sal = df["sal_weighted"].sum() / df["area_sqm"].sum()
                # insert into the argo_monthly_data table
                db.execute("INSERT INTO argo_monthly_data (source, MEAN_ARGO_TEMPERATURE_ANOMALY, MEAN_ARGO_SALINITY_ANOMALY, TIME) VALUES (?, ?, ?, ?);", (source_value, mean_temp, mean_sal, time_value))
                db.commit()

                            # longitude_value = longitude[longitude_idx].item()
                            # latitude_value = latitude[latitude_idx].item()
                            # pressure_value = pressure[pressure_idx].item()

                            # if filetype == "Temperature":
                            #     temperature_value = argo_temp[
                            #         time_idx, pressure_idx, latitude_idx, longitude_idx
                            #     ].item()
                            #     # insert into the sqlite db, making sure that we also insert a wkt for the geometry
                            #     data.append(
                            #         (
                            #             source_value,
                            #             temperature_value,
                            #             longitude_value,
                            #             latitude_value,
                            #             pressure_value,
                            #             time_value,
                            #             longitude_value,
                            #             latitude_value,
                            #         )
                            #     )

                            #     # new_row = {
                            #     #     "source": source_value,
                            #     #     "ARGO_TEMPERATURE_ANOMALY": temperature_value,
                            #     #     "LONGITUDE": longitude_value,
                            #     #     "LATITUDE": latitude_value,
                            #     #     "PRESSURE": pressure_value,
                            #     #     "TIME": time_value,
                            #     # }
                            #     # data_2004[
                            #     #     f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
                            #     # ] = new_row
                            #     # print(new_row)
                            # else:
                            #     sal_value = argo_sal[
                            #         time_idx, pressure_idx, latitude_idx, longitude_idx
                            #     ].item()
                            #     # update the row
                            #     updates.append(
                            #         (
                            #             sal_value,
                            #             longitude_value,
                            #             latitude_value,
                            #             pressure_value,
                            #             time_value,
                            #         ),
                            #     )
                            #     # data_2004[
                            #     #     f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
                            #     # ]["ARGO_SALINITY_ANOMALY"] = sal_value
                # if data:
                #     cur.executemany(
                #         "INSERT INTO argo_data (source, ARGO_TEMPERATURE_ANOMALY, LONGITUDE, LATITUDE, PRESSURE, TIME, geometry) VALUES (?, ?, ?, ?, ?, ?, MakePoint(?, ?, 4326));",
                #         data,
                #     )
                # if updates:
                #     cur.executemany(
                #         "UPDATE argo_data SET ARGO_SALINITY_ANOMALY = ? WHERE LONGITUDE = ? AND LATITUDE = ? AND PRESSURE = ? AND TIME = ?;",
                #         updates,
                #     )
                # cur.execute("COMMIT")
                # print(
                #     db.execute(
                #         "SELECT source, ARGO_TEMPERATURE_ANOMALY, ARGO_SALINITY_ANOMALY, LONGITUDE, LATITUDE, PRESSURE, TIME, AsWKT(geometry) FROM argo_data WHERE LONGITUDE = ? AND LATITUDE = ? AND PRESSURE = ? AND TIME = ?;",
                #         (longitude_value, latitude_value, pressure_value, time_value),
                #     ).fetchone()
                # )
                # break
                # sys.exit(0)
                # data_2004[f"{time}, {pressure}, {latitude}, {longitude}"] = {"source": "original",
                # "ARGO_TEMPERATURE_ANOMALY": rootgrp["ARGO_TEMPERATURE_ANOMALY"][longitude, latitude, pressure, time]}
    # Iterate over the rocksdict and make a geodataframe
    # data_list = []
    # for k, v in tqdm(data_2004.items(), desc="Making geodataframe"):
    #     data_dict = v  # this is your original data dictionary
    #     data_list.append(data_dict)

    # df = pd.DataFrame(data_list)

    # # make a geodataframe from data_2004
    # # df = pd.DataFrame.from_dict(dict_data_2004, orient="index")
    # # make a geometry column of Points from the Long Lat
    # df["geometry"] = df.apply(
    #     lambda row: Point(row["LONGITUDE"], row["LATITUDE"]), axis=1
    # )
    # # make a geodataframe from the dataframe
    # gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
    # return gdf


def plot_temperature_distribution(sqlitedb="data/argo_data.db", output_file="output/temperature_distribution.png""):
    """
    Plot the distribution of the ocean temperature data.

    recreate "Timeseries of 0-2000dbar Global-Average Open Ocean Potential Temperature (oC) (Spatial mask defined within RG netCDF file. Marginal Seas and the Arctic Ocean are excluded from the calculation) over the 180 months of the RG Climatology (2004-2018) [black line] and 2019-2020 monthly extensions [red line]. The 12-Month box-car average is shown in green."
    """
    # make a dataframe from the sqlite db in data_directory
    # CREATE TABLE argo_monthly_data (
    # source TEXT,
    # MEAN_ARGO_TEMPERATURE_ANOMALY REAL,
    # MEAN_ARGO_SALINITY_ANOMALY REAL,
    # TIME TEXT PRIMARY KEY);
    with sqlite3.connect(sqlitedb) as db:
        query = "SELECT source, MEAN_ARGO_TEMPERATURE, MEAN_ARGO_SALINITY, TIME from argo_monthly_data"

        df = pd.read_sql(query, conn)

        # make sure that the time output from sqlite is mapped properly to python
        df["TIME"] = pd.to_datetime(df["TIME"])

        # we are now going to make two plots, mean temperature over time, and mean salinity over time using seaborn

        # make a figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

        # plot the mean temperature over time
        sns.lineplot(data=df, x="TIME", y="MEAN_ARGO_TEMPERATURE", ax=ax1)
        ax1.set_title("Mean Temperature over Time")
        ax1.set_ylabel("Temperature (C)")

        # plot the mean salinity over time
        sns.lineplot(data=df, x="TIME", y="MEAN_ARGO_SALINITY", ax=ax2)
        ax2.set_title("Mean Salinity over Time")
        ax2.set_ylabel("Salinity (PSU)")

        # rotate dates on x-axis
        ax1.set_xticklabels(
            ax1.get_xticklabels(), rotation=45, horizontalalignment="right"
        )

        ax2.set_xticklabels(
            ax2.get_xticklabels(), rotation=45, horizontalalignment="right"
        )

        plt.tight_layout()
        
        # Save fig to output file
        fig.savefig(output_file)
        
    return df




# ---------------------------------------
# 1.3 Sea Ice Concentration Data (netCDF)
# ---------------------------------------


def load_sea_ice_concentration_data(filepath):
    """
    Load the Sea Ice Concentration netCDF file and return the data array.
    """
    pass


def plot_sea_ice_concentration_timeseries(data_array):
    """
    Plot a time series of the sea ice concentration.
    """
    pass


# ---------------------------------------
# 1.4 Sea Ice Extent Data (GeoTIFF, shapefile, CSV)
# ---------------------------------------


def load_sea_ice_extent_data(filepath):
    """
    Load the Sea Ice Extent data and return it.
    (Choose the most suitable format, e.g., CSV)
    """
    pass


def plot_sea_ice_extent_trends(data):
    """
    Plot the trend of sea ice extent over time.
    """
    pass


# ---------------------------------------
# Main execution
# ---------------------------------------

if __name__ == "__main__":
    OUTPUT = "output"
    # delete and recreate output
    if os.path.exists(OUTPUT):
        shutil.rmtree(OUTPUT)
    os.mkdir(OUTPUT)

    # Load datasets
    sea_ice_extent_df = load_sea_ice_extent_csv("data/S_seaice_extent_daily_v3.0.csv")

    # try to load the argo data geodataframe from parquet. If not, regenerate it.
    # try:
    #     argo_data = pd.read_parquet("data/argo_data.parquet")
    # except:

    argo_data = load_argo_data(directory_path="data/ArgoClim")
    # run once.

    # we're going to need to use a db, not a dataframe.
    # argo_data.to_parquet("data/argo_data.parquet")

    # sea_ice_concentration_data = load_sea_ice_concentration_data(
    #     "path_to_sea_ice_concentration_data.nc"
    # )
    # sea_ice_extent_data = load_sea_ice_extent_data(
    #     "path_to_sea_ice_extent_data.csv_or_other_format"
    # )

    # Preliminary visualizations
    plot_sea_ice_extent_timeseries(
        directory_path="data/ArgoClim", 
        output_dir=f"{OUTPUT}/argo_demo.png"
    )
    # plot_temperature_distribution(argo_data)
    plot_sea_ice_concentration_timeseries(sea_ice_concentration_data)
    plot_sea_ice_extent_trends(sea_ice_extent_data)
