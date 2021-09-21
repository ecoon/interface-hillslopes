import os
import numpy as np

import scipy.stats
import collections
import fiona, rasterio, shapely
import rasterio.warp

import workflow
import workflow.crs
import workflow.warp

import landcover
import datetime


def get_filenames(huc, huc_directory, raster_extension='tif'):
    """Set up package directory for one huc to return a dictionary of filenames"""
    ppm_dir = os.path.join(huc_directory,'data_preprocessed-meshing')
    
    filenames = dict()
    def register_shapefile(key, filename):
        filenames[key] = os.path.join(ppm_dir, filename)+'.shp'

    def register_raster(key, filename):
        filenames[key] = os.path.join(ppm_dir, filename)+'.'+raster_extension

    # the HUC shapefile, in native and projected CRS
    register_shapefile('huc', f'huc_{huc}')
    register_shapefile('huc_proj', f'huc_{huc}_proj')

    # the HUC DEM, land_cover, slope, and aspect rasters
    register_raster('dem', f'huc_{huc}_dem')  # in units of [m] above sea level
    #register_raster('dem_filled', f'huc_{huc}_dem_filled')
    #register_raster('d8', f'huc_{huc}_d8')
    #register_raster('d8i', f'huc_{huc}_d8i')
    #register_raster('weights', f'huc_{huc}_flowpath_weights')
    register_raster('land_cover', f'huc_{huc}_landcover')
    register_raster('slope', f'huc_{huc}_slope')  # in units of [-], e.g. rise over run, NOT % slope
    register_raster('aspect', f'huc_{huc}_aspect') # in units of degrees clockwise from North (0 = N, 90 = E, 180 = S)

    # raster and shapefiles of the stream network
    register_raster('streams_raster', f'huc_{huc}_streams')
    register_shapefile('streams', f'huc_{huc}_streams')  # the shapefile extracted from the above raster
    register_shapefile('streams_network', f'huc_{huc}_streams_network') # simplified to a network for MOSART
    register_shapefile('streams_network_proj', f'huc_{huc}_streams_network_proj') # projected form of the above

    # delineated subcatchments within the HUC
    register_shapefile('subcatchments', f'huc_{huc}_subcatchments')  # delineated subcatchments
    
    filenames['flowpaths'] = os.path.join(ppm_dir, 'flowpaths', f'hs_{{}}_flowpaths.pkl')
    register_raster('flowpath_length', f'huc_{huc}_flowpath_lengths') # flowpaths within the delineated subcatchments
    register_raster('elev_above_streams', f'huc_{huc}_elev_above_streams')

    filenames['daymet'] = os.path.join(huc_directory, 'daymet', f'huc_{huc}_subcatchment{{}}_1980_2020.h5')
    filenames['mesh'] = os.path.join(huc_directory, 'mesh', f'sag_hillslope{{}}.exo')
    
    
    return filenames


def save_shapefile(filename, shps, crs, extra_properties=None):
        if len(shps) == 0:
            return
        
        schema = dict()
        if type(shps[0]) is shapely.geometry.Polygon:
            schema['geometry'] = 'Polygon'
        elif type(shps[0]) is shapely.geometry.LineString:
            schema['geometry'] = 'LineString'
        else:
            raise RuntimeError('Currently this function only writes Polygon or LineString types')
            
        # set up the properties' schema, used for open and save geodata by fiona
        schema['properties'] = collections.OrderedDict()
        def register_type(key,atype):
            if atype is int:
                schema['properties'][key] = 'int'
            elif atype is str:
                schema['properties'][key] = 'str'
            elif atype is float:
                schema['properties'][key] = 'float'
            else:
                pass
        if extra_properties is None:
            extra_properties = dict()
        for key, val in extra_properties.items():
            register_type(key, type(val))
            
        try:
            shp_property = shps[0].properties
        except AttributeError:
            pass
        else:
            for key, val in shp_property.items():
                register_type(key, type(val))
                
        
        with fiona.open(filename, 'w',
                        driver='ESRI Shapefile',
                        schema=schema,
                        crs=workflow.crs.to_fiona(crs),
                        crs_wkt=workflow.crs.to_wkt(crs)) as c:
            for shp in shps:
                props = extra_properties.copy()
                try:
                    props.update(shp.properties)
                except AttributeError:
                    pass
                
                for key in list(props.keys()):
                    if key not in schema['properties']:
                        props.pop(key)
                        
                c.write({'geometry': shapely.geometry.mapping(shp),
                         'properties': props})
                
                
def save_demfile(filename, dem_profile, dem_raster):
    rio_profile = dict(dem_profile).copy()
    rio_profile.pop('blockxsize')
    rio_profile.pop('blockysize')
    rio_profile.pop('tiled')
    rio_profile['nodata'] = -9999.0
    
    rio_dem = np.where(np.isnan(dem_raster), rio_profile['nodata'], dem_raster).astype(rio_profile['dtype'])
    with rasterio.open(filename, 'w', **rio_profile) as dst:
        dst.write(rio_dem,1)
        
        
# Average aspect across the domain
def meanAspect(aspect):
    '''
    # tests...
      assert(meanAspect(np.array([90,90,90])) == 90)
      assert(meanAspect(np.array([0,0,0])) == 0)
      assert(meanAspect(np.array([180,180,180])) == 180)
      assert(meanAspect(np.array([270,270,270])) == 270)
      assert(meanAspect(np.array([89,90,91])) == 90)
      assert(meanAspect(np.array([179,180,181])) == 180)
      assert(meanAspect(np.array([269,270,271])) == 270)
      assert(meanAspect(np.array([359,0,1])) == 0)
    '''
    
    a = aspect[~np.isnan(aspect)]
    a = np.where(a > 180, a - 360, a)
    sina = np.sin(np.radians(a))
    cosa = np.cos(np.radians(a))
    avg_aspect_radians = np.arctan2(sum(sina), sum(cosa))
    if avg_aspect_radians < 0:
        avg_aspect_radians = 2*np.pi + avg_aspect_radians
    return np.degrees(avg_aspect_radians)
    

    
# Load a subcatchment shape
def loadSubcatchmentShape(filename, subcatch_id):
    subcatch_crs, subcatchments = workflow.get_shapes(filename)
    matches = [sc for sc in subcatchments if sc.properties['hs_id'] == subcatch_id]
    if len(matches) != 1:
        raise RuntimeError("Found invalid subcatchments file or no subcatchment of this id.")
    return subcatch_crs, matches[0]


# Load a subcatchment raster
def loadSubcatchmentRaster(filename, subcatch_shape, subcatch_crs, nanit=True):
    profile, raster = workflow.get_raster_on_shape(filename, subcatch_shape, subcatch_crs, mask=True)
    if nanit:
        raster[raster==profile['nodata']] = np.nan
    return profile, raster


# Determine hillslope properties for a given subcatchment.
def parameterizeSubcatchment(filenames, huc, subcatch_id,
                             target_crs=workflow.crs.default_alaska_crs(),
                             hillslope_keep_fraction=0.95,
                             hillslope_bin_dx=100):
    
    # Find the given subcatchment shape
    subcatch_crs, subcatch_shape = loadSubcatchmentShape(filenames['subcatchments'], subcatch_id)

    # Create a dictionary icluding all hillslope info and parameters
    hillslope = dict()
    hillslope['huc'] = huc
    hillslope['subcatchment_id'] = subcatch_id
    hillslope['subcatchment'] = subcatch_shape
    hillslope['centroid'] = subcatch_shape.centroid.coords[0]     # lat/long, required to get meteorological data
    hillslope['total_area'] = workflow.warp.shply(subcatch_shape, subcatch_crs, target_crs).area   # m^2
    
    # Procedures to determine a single hillslope profile geometry for each subcatchment 
    # Most of the parameters are based on bins in length of flowpath to the stream network
    
    # 1. Load raster of flowpath lengths for the given subcatchment
    fp_profile, fp_lengths = loadSubcatchmentRaster(filenames['flowpath_length'], subcatch_shape, subcatch_crs)
    subcatch_mask = ~np.isnan(fp_lengths)
    hillslope['raster_profile'] = fp_profile
    
    # 2. Ensure raster of flowpath lengths >= 0, and sort flowpath lengths
    assert(fp_lengths[subcatch_mask].min() >= -1.e-3)
    fp_lengths[fp_lengths < 0] = 0.
    fp_lengths_masked = fp_lengths[subcatch_mask]
    fp_lengths_masked_sorted = np.sort(fp_lengths_masked)
    
    # 3. Set bin parameters based on flowpath length
    hillslope['total_length'] = fp_lengths_masked_sorted[-1] * hillslope_keep_fraction
    hillslope['num_bins'] = int(np.round(hillslope['total_length'] / hillslope_bin_dx))
    hillslope['bin_dx'] = hillslope['total_length'] / hillslope['num_bins']
    assert(hillslope['num_bins'] > 3)
    
    # 4. Bin by flowpath length and save to dictionary
    pixel_bin_raster = np.nan * np.ones(fp_lengths.shape, 'd')
    pixel_bin_counts = np.zeros((hillslope['num_bins'],),'i')
    
    i = 0
    interval = 1
    while i < hillslope['num_bins']:
        bin_start = i * hillslope['bin_dx']
        bin_end = (i+interval) * hillslope['bin_dx']
        local_mask = (fp_lengths >= bin_start) & (fp_lengths < bin_end)
        pixel_bin_raster[local_mask] = i
        pixel_bin_counts[i] = np.count_nonzero(local_mask)
    
        if pixel_bin_counts[i] > 0:
            i += 1
            if i-1+interval>hillslope['num_bins']:
                break
        else:
            interval += 1
            continue
    
    hillslope['bin_counts'] = pixel_bin_counts[np.nonzero(pixel_bin_counts)]
    hillslope['num_bins'] = len(hillslope['bin_counts'])
    hillslope['bins'] = pixel_bin_raster
    
#     assert(hillslope['bin_counts'].min() > 0)
#     hillslope['bin_counts'] = pixel_bin_counts   # THE OLD
        
    
    
    # 5. Average elevation (height over stream) for each bin and save to dictionary
    _, elevs = loadSubcatchmentRaster(filenames['elev_above_streams'], subcatch_shape, subcatch_crs)
    elev_bins = np.zeros((hillslope['num_bins'],),'d')
    for i in range(hillslope['num_bins']):
        elev_bins[i] = elevs[(hillslope['bins'] == i) & (~np.isnan(elevs))].mean()

    hillslope['elevation'] = elev_bins
    
    # 6. Average aspect across the entire subcatchment and save to dictionary
    _, aspects = loadSubcatchmentRaster(filenames['aspect'], subcatch_shape, subcatch_crs)
    hillslope['aspect'] = meanAspect(aspects)
    
    # 7. Get land cover using mode for each bin and save to dictionary
    lc_profile, lc = loadSubcatchmentRaster(filenames['land_cover'], subcatch_shape, subcatch_crs, False) # don't nan-it
    assert(lc_profile['nodata'] == 255) # uint8, nan is -1 == 255
    # classify land cover into veg classes
    hillslope['lc'] = lc
    num_missing, lc = landcover.classifyVegetation(lc)
    lc_bins = np.zeros((hillslope['num_bins'],),'i')
    
    for i in range(hillslope['num_bins']):
        bin_lc = lc[hillslope['bins'] == i]
        assert(len(bin_lc) > 0)
        lc_bins[i] = scipy.stats.mode(bin_lc.ravel())[0][0]
        
    hillslope['land_cover'] = lc_bins
    hillslope['land_cover_raster'] = landcover.veg2img(lc)

    return hillslope




# Download DayMet
def downloadDaymet(hillslope_pars, raw_directory, save_filename, 
                   start=datetime.date(1980,1,1),
                   end=datetime.date(2020,12,31)):
    import daymet_to_ats
    start, end = daymet_to_ats.validate_start_end(start, end)
    lon, lat = hillslope_pars['centroid']
    daymet = daymet_to_ats.download_daymet(raw_directory, lat, lon, start, end)
    ats = daymet_to_ats.daymet_to_ats(daymet)
    attrs = daymet_to_ats.daymet_attrs(lat, lon, start, end)
    
    daymet_to_ats.write_ats(ats, attrs,save_filename)
                          
    