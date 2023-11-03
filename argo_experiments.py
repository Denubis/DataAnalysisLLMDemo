
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