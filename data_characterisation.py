#!/usr/bin/env python3
# Import necessary libraries
import pandas as pd

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
from rocksdict import Rdict

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


def load_argo_data(directory_path):
    """
    Load the Argo netCDF file and return relevant data arrays.
    """
    # let's see what the monthly data looks like

    def walktree(top):
        yield top.groups.values()
        for value in top.groups.values():
            yield from walktree(value)

    # test_path = "data/ArgoClim/RG_ArgoClim_201901_2019.nc"
    # with Dataset(test_path, "r", format="NETCDF3_CLASSIC") as rootgrp:
    #     # print(rootgrp.data_model)
    #     # https://unidata.github.io/netcdf4-python/#creatingopeningclosing-a-netcdf-file
    #     # for children in walktree(rootgrp):
    #     #     for child in children:
    #     #         print(child)
    #     # there are no groups in this data...
    #     # print(rootgrp.dimensions)
    #     # {'LONGITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LONGITUDE', size = 360,
    #     #  'LATITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LATITUDE', size = 145,
    #     #  'PRESSURE': <class 'netCDF4._netCDF4.Dimension'>: name = 'PRESSURE', size = 58,
    #     #  'TIME': <class 'netCDF4._netCDF4.Dimension'>: name = 'TIME', size = 1}
    #     # print(rootgrp.variables)
    #     # float32 LONGITUDE(LONGITUDE)
    #     #     units: degrees_east
    #     #     modulo: 360.0
    #     #     point_spacing: even
    #     #     axis: X
    #     # unlimited dimensions:
    #     # current shape = (360,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'LATITUDE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 LATITUDE(LATITUDE)
    #     #     units: degrees_north
    #     #     point_spacing: even
    #     #     axis: Y
    #     # unlimited dimensions:
    #     # current shape = (145,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'PRESSURE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 PRESSURE(PRESSURE)
    #     #     units: dbar
    #     #     positive: down
    #     #     point_spacing: uneven
    #     #     axis: Z
    #     # unlimited dimensions:
    #     # current shape = (58,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'TIME': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 TIME(TIME)
    #     #     units: months since 2004-01-01 00:00:00
    #     #     time_origin: 01-JAN-2004 00:00:00
    #     #     axis: T
    #     # unlimited dimensions:
    #     # current shape = (1,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'ARGO_TEMPERATURE_ANOMALY': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_TEMPERATURE_ANOMALY(TIME, PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: degree celcius (ITS-90)
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO TEMPERATURE ANOMALY defined by 2004 - 2018 RG CLIMATOLOGY
    #     # unlimited dimensions:
    #     # current shape = (1, 58, 145, 360)
    #     # filling on, 'ARGO_SALINITY_ANOMALY': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_SALINITY_ANOMALY(TIME, PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: Practical Salinity Scale 78
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO SALINITY ANOMALY defined by 2004 - 2018 RG CLIMATOLOGY
    #     # unlimited dimensions:
    #     # current shape = (1, 58, 145, 360)
    #     # filling on}
    #     # print("attrs", rootgrp.__dict__)
    #     # for name in rootgrp.ncattrs():
    #     #     print("Global attr {} = {}".format(name, getattr(rootgrp, name)))

    #     # print the first row of data

    # print("salinity")
    # test_path3 = "data/ArgoClim/RG_ArgoClim_Salinity_2019.nc"
    # with Dataset(test_path3, "r", format="NETCDF3_CLASSIC") as rootgrp3:
    #     pass
    #     # print( rootgrp3.variables)

    #     # {'LONGITUDE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 LONGITUDE(LONGITUDE)
    #     #     units: degrees_east
    #     #     modulo: 360.0
    #     #     point_spacing: even
    #     #     axis: X
    #     # unlimited dimensions:
    #     # current shape = (360,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'LATITUDE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 LATITUDE(LATITUDE)
    #     #     units: degrees_north
    #     #     point_spacing: even
    #     #     axis: Y
    #     # unlimited dimensions:
    #     # current shape = (145,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'PRESSURE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 PRESSURE(PRESSURE)
    #     #     units: dbar
    #     #     positive: down
    #     #     point_spacing: uneven
    #     #     axis: Z
    #     # unlimited dimensions:
    #     # current shape = (58,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'TIME': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 TIME(TIME)
    #     #     units: months since 2004-01-01 00:00:00
    #     #     time_origin: 01-JAN-2004 00:00:00
    #     #     axis: T
    #     # unlimited dimensions: TIME
    #     # current shape = (180,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'ARGO_SALINITY_MEAN': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_SALINITY_MEAN(PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: Practical Salinity Scale 78
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO SALINITY MEAN Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on, 'ARGO_SALINITY_ANOMALY': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_SALINITY_ANOMALY(TIME, PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: Practical Salinity Scale 78
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO SALINITY ANOMALY defined by Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY
    #     # unlimited dimensions: TIME
    #     # current shape = (180, 58, 145, 360)
    #     # filling on, 'BATHYMETRY_MASK': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 BATHYMETRY_MASK(PRESSURE, LATITUDE, LONGITUDE)
    #     #     missing_value: -9.0
    #     #     _FillValue: -9.0
    #     #     long_name: BATHYMETRY MASK
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on, 'MAPPING_MASK': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 MAPPING_MASK(PRESSURE, LATITUDE, LONGITUDE)
    #     #     missing_value: -9.0
    #     #     _FillValue: -9.0
    #     #     long_name: MAPPING MASK: pressure limits of mapping can be shallower than 2000dbar in marginal seas
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on}

    # print("\n\n***\n\n")
    # test_path4 = "data/ArgoClim/RG_ArgoClim_Temperature_2019.nc"
    # with Dataset(test_path4, "r", format="NETCDF3_CLASSIC") as rootgrp4:
    #     # print(rootgrp4.variables)
    #     pass

    #     # {'LONGITUDE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 LONGITUDE(LONGITUDE)
    #     #     units: degrees_east
    #     #     modulo: 360.0
    #     #     point_spacing: even
    #     #     axis: X
    #     # unlimited dimensions:
    #     # current shape = (360,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'LATITUDE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 LATITUDE(LATITUDE)
    #     #     units: degrees_north
    #     #     point_spacing: even
    #     #     axis: Y
    #     # unlimited dimensions:
    #     # current shape = (145,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'PRESSURE': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 PRESSURE(PRESSURE)
    #     #     units: dbar
    #     #     positive: down
    #     #     point_spacing: uneven
    #     #     axis: Z
    #     # unlimited dimensions:
    #     # current shape = (58,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'TIME': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 TIME(TIME)
    #     #     units: months since 2004-01-01 00:00:00
    #     #     time_origin: 01-JAN-2004 00:00:00
    #     #     axis: T
    #     # unlimited dimensions: TIME
    #     # current shape = (180,)
    #     # filling on, default _FillValue of 9.969209968386869e+36 used, 'ARGO_TEMPERATURE_MEAN': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_TEMPERATURE_MEAN(PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: degree celcius (ITS-90)
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO TEMPERATURE MEAN Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on, 'ARGO_TEMPERATURE_ANOMALY': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 ARGO_TEMPERATURE_ANOMALY(TIME, PRESSURE, LATITUDE, LONGITUDE)
    #     #     units: degree celcius (ITS-90)
    #     #     missing_value: -999.0
    #     #     _FillValue: -999.0
    #     #     long_name: ARGO TEMPERATURE ANOMALY defined by Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY
    #     # unlimited dimensions: TIME
    #     # current shape = (180, 58, 145, 360)
    #     # filling on, 'BATHYMETRY_MASK': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 BATHYMETRY_MASK(PRESSURE, LATITUDE, LONGITUDE)
    #     #     missing_value: -9.0
    #     #     _FillValue: -9.0
    #     #     long_name: BATHYMETRY MASK
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on, 'MAPPING_MASK': <class 'netCDF4._netCDF4.Variable'>
    #     # float32 MAPPING_MASK(PRESSURE, LATITUDE, LONGITUDE)
    #     #     missing_value: -9.0
    #     #     _FillValue: -9.0
    #     #     long_name: MAPPING MASK: pressure limits of mapping can be shallower than 2000dbar in marginal seas
    #     # unlimited dimensions:
    #     # current shape = (58, 145, 360)
    #     # filling on}

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

    data_2004 = Rdict("/tmp/data_2004.dd")

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
            time = num2date(
                times=rootgrp["TIME"][:],
                units="months since 2004-01-01 00:00:00",
                calendar="360_day",
            )

            for longitude_idx in tqdm(
                range(len(longitude)), desc=f"Loading longitudes", leave=False
            ):
                longitude_value = longitude[longitude_idx]

                for latitude_idx in range(len(latitude)):
                    latitude_value = latitude[latitude_idx]

                    if latitude_value > -40:
                        # print("skipping latitude", latitude_value)
                        continue
                    for pressure_idx in range(len(pressure)):
                        for time_idx in range(len(time)):
                            source_value = "extension"
                            time_value = time[time_idx]
                            pressure_value = pressure[pressure_idx]
                            temperature_value = argo_temp[
                                time_idx, pressure_idx, latitude_idx, longitude_idx
                            ]

                            sal_value = argo_sal[
                                time_idx, pressure_idx, latitude_idx, longitude_idx
                            ]

                            new_row = {
                                "source": source_value,
                                "ARGO_TEMPERATURE_ANOMALY": temperature_value,
                                "ARGO_SALINITY_ANOMALY": sal_value,
                                "LONGITUDE": longitude_value,
                                "LATITUDE": latitude_value,
                                "PRESSURE": pressure_value,
                                "TIME": time_value,
                            }
                            data_2004[
                                f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
                            ] = new_row
                            # print(new_row)
                            # sys.exit(0)

    for filetype in ["Temperature", "Salinity"]:
        with Dataset(
            f"{directory_path}/RG_ArgoClim_{filetype}_2019.nc",
            "r",
            format="NETCDF3_CLASSIC",
        ) as nc_file:
            print(nc_file.variables)
            # We need to iterate over all four dimensions to add the data to the dict
            #  {'LONGITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LONGITUDE', size = 360, 'LATITUDE': <class 'netCDF4._netCDF4.Dimension'>: name = 'LATITUDE', size = 145, 'PRESSURE': <class 'netCDF4._netCDF4.Dimension'>: name = 'PRESSURE', size = 58, 'TIME': <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'TIME', size = 180}
            if filetype == "Temperature":
                argo_temp = nc_file.variables["ARGO_TEMPERATURE_ANOMALY"][:]
            else:
                argo_sal = nc_file.variables["ARGO_SALINITY_ANOMALY"][:]
            # argo_temp.shape: (180, 58, 145, 360)
            longitude = nc_file.variables["LONGITUDE"][:]
            # print(longitude.shape)
            # (360,)
            latitude = nc_file.variables["LATITUDE"][:]
            print(latitude.shape)
            # (145,)
            pressure = nc_file.variables["PRESSURE"][:]
            print(pressure.shape)
            # (58,)
            time = num2date(
                times=nc_file["TIME"][:],
                units="months since 2004-01-01 00:00:00",
                calendar="360_day",
            )
            print(time.shape)
            # (180,)

            # iterate over all of these to populate a dictionary so that we may join two datasources together
            for longitude_idx in tqdm(
                range(len(longitude)), desc=f"Loading {filetype} data"
            ):
                longitude_value = longitude[longitude_idx]

                for latitude_idx in tqdm(range(len(latitude)), leave=False):
                    latitude_value = latitude[latitude_idx]

                    if latitude_value > -40:
                        # print("skipping latitude", latitude_value)
                        continue
                    for pressure_idx in range(len(pressure)):
                        for time_idx in range(len(time)):
                            # print(longitude_idx, latitude_idx, pressure_idx, time_idx)

                            source_value = "original"
                            time_value = time[time_idx]
                            longitude_value = longitude[longitude_idx]
                            latitude_value = latitude[latitude_idx]
                            pressure_value = pressure[pressure_idx]

                            if filetype == "Temperature":
                                temperature_value = argo_temp[
                                    time_idx, pressure_idx, latitude_idx, longitude_idx
                                ]
                                new_row = {
                                    "source": source_value,
                                    "ARGO_TEMPERATURE_ANOMALY": temperature_value,
                                    "LONGITUDE": longitude_value,
                                    "LATITUDE": latitude_value,
                                    "PRESSURE": pressure_value,
                                    "TIME": time_value,
                                }
                                data_2004[
                                    f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
                                ] = new_row
                                print(new_row)
                            else:
                                sal_value = argo_sal[
                                    time_idx, pressure_idx, latitude_idx, longitude_idx
                                ]
                                data_2004[
                                    f"{time_value}, {pressure_value}, {latitude_value}, {longitude_value}"
                                ]["ARGO_SALINITY_ANOMALY"] = sal_value
                        # sys.exit(0)
                        # data_2004[f"{time}, {pressure}, {latitude}, {longitude}"] = {"source": "original",
                        # "ARGO_TEMPERATURE_ANOMALY": rootgrp["ARGO_TEMPERATURE_ANOMALY"][longitude, latitude, pressure, time]}

    # make a geodataframe from data_2004
    df = pd.DataFrame.from_dict(data_2004, orient="index")
    # make a geometry column of Points from the Long Lat
    df["geometry"] = df.apply(
        lambda row: Point(row["LONGITUDE"], row["LATITUDE"]), axis=1
    )
    # make a geodataframe from the dataframe
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
    return gdf


def plot_temperature_distribution(data_array):
    """
    Plot the distribution of the ocean temperature data.

    recreate "Timeseries of 0-2000dbar Global-Average Open Ocean Potential Temperature (oC) (Spatial mask defined within RG netCDF file. Marginal Seas and the Arctic Ocean are excluded from the calculation) over the 180 months of the RG Climatology (2004-2018) [black line] and 2019-2020 monthly extensions [red line]. The 12-Month box-car average is shown in green."
    """
    pass


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
    try:
        argo_data = pd.read_parquet("data/argo_data.parquet")
    except:
        argo_data = load_argo_data(directory_path="data/ArgoClim")
        argo_data.to_parquet("data/argo_data.parquet")

    sea_ice_concentration_data = load_sea_ice_concentration_data(
        "path_to_sea_ice_concentration_data.nc"
    )
    sea_ice_extent_data = load_sea_ice_extent_data(
        "path_to_sea_ice_extent_data.csv_or_other_format"
    )

    # Preliminary visualizations
    plot_sea_ice_extent_timeseries(
        sea_ice_extent_df, f"{OUTPUT}/sea_ice_extent_timeseries.png"
    )
    plot_temperature_distribution(argo_data)
    plot_sea_ice_concentration_timeseries(sea_ice_concentration_data)
    plot_sea_ice_extent_trends(sea_ice_extent_data)
