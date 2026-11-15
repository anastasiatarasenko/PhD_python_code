import datetime

def write_array_to_nc(data_array, lats, lons, nc_path):    

    dataset = netCDF4.Dataset(nc_path , 'w', format='NETCDF4')

    dataset.createDimension('r', lats.shape[0])
    dataset.createDimension('c', lons.shape[1])
    
    # attributes
    dataset.title = 'Sentinel-1'
    dataset.title_short_name = 'Sentinel-1'
    dataset.institution = 'Satellite Oceanography Laboratory (SOLab) of the Russian State Hydrometeorological University'
    dataset.creator_url = 'http://solab.rshu.ru/'
    dataset.source = ''
    dataset.source_provider = ''
    dataset.source_filename = ''
    dataset.history = unicode(datetime.datetime.now())[:-7]
    dataset.references = ''
    dataset.comment = ''
    dataset.conventions = 'CF-1.7'
    dataset.start_date = ''

    lat_var = dataset.createVariable('lat', 'f4', ('r','c'),
                            zlib=True, complevel=6)
    lat_var[:] = lats[:]
    lat_var.units = 'degrees_north'
    lat_var.scale_factor = 9.99999975e-06

    lon_var = dataset.createVariable('lon', 'f4', ('r','c'),
                            zlib=True, complevel=6)
    lon_var[:] = lons[:]
    lon_var.units = 'degrees_east'
    lon_var.scale_factor = 9.99999975e-06
    
    if isinstance(data_array, dict):
        for key, data in data_array.iteritems():
            key = dataset.createVariable(key, 'f4', ('r', 'c'),
                                         zlib=True, complevel=6)
            key[:] = data[:]
            
    else:    
        data_var = dataset.createVariable('data', 'f4', ('r','c'),
                             zlib=True, complevel=6)
        data_var[:] = data_array[:]
    
    

    dataset.close()