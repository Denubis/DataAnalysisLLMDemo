#!/usr/bin/env bash

set -euo pipefail

mkdir -p data
cp original_data/S_seaice_extent_daily_v3.0.csv data/

cd data

# If this doesn't exist, download it
if [ ! -f RG_ArgoClim_Temperature_2019.nc ]; then
    wget https://sio-argo.ucsd.edu/pub/www-argo/RG/RG_ArgoClim_Temperature_2019.nc.gz 
    gunzip RG_ArgoClim_Temperature_2019.nc.gz
fi

if [ ! -f RG_ArgoClim_Salinity_2019.nc ]; then
    wget https://sio-argo.ucsd.edu/pub/www-argo/RG/RG_ArgoClim_Salinity_2019.nc.gz
    gunzip RG_ArgoClim_Salinity_2019.nc.gz
fi



# we need to use wget to get all .nc files in: https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly/
mkdir -p G02202_V4
cd G02202_V4
wget --no-check-certificate --reject "index.html*" -np -e robots=off -r -nH --cut-dirs=3 -R index.html* https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly/
cd ..
# Let's get the montly south data recursively, four folders
mkdir -p G02135
cd G02135
wget --no-check-certificate --reject "index.html*" -np -e robots=off  -r -np -nH --cut-dirs=3 -R index.html* https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/
cd ..
# this data goes from 201901 to 202309. Please iterate over that range and wget, then ungzip each one.
# wget https://sio-argo.ucsd.edu/pub/www-argo/RG/RG_ArgoClim_201901_2019.nc.gz

# derive the current yearmonth from the current date

export maxyearmonth=$(date +%Y%m)

# make a range from 2019 01 to maxyearmonth for the bash for loop.
# this is year 2019 through current year
# month 01 through 12
# make a oneline python script to print the range

export range=$(python -c "print(' '.join(['%s%02d' % (y, m) for y in range(2019, int('$maxyearmonth') // 100 + 1) for m in range(1, 13)]))")





for i in $range; do
    echo "Downloading $i"
    if [ ! -f RG_ArgoClim_${i}_2019.nc ]; then
        wget https://sio-argo.ucsd.edu/pub/www-argo/RG/RG_ArgoClim_${i}_2019.nc.gz
        gunzip RG_ArgoClim_${i}_2019.nc.gz
    fi
done
