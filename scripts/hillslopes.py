"""Functions for parameterizing hillslope simulations from catchments."""

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import scipy.optimize, scipy.signal, scipy.stats, scipy.integrate
import collections
import logging
import fiona, rasterio, shapely
import rasterio.warp
import attr

import workflow
import workflow.crs
import workflow.warp
import workflow.source_list
import workflow.utils
import workflow.ui
import workflow.conf
import workflow.mesh

import land_cover

def get_filenames(huc, package_directory, raster_extension='tif'):
    """Set up the package directory and return a dictionary of filenames"""
    huc_directory = os.path.join(package_directory, huc)
    if not os.path.isdir(huc_directory):
        os.mkdir(huc_directory)

    filenames = dict()
    def register_shapefile(key, filename):
        filenames[key] = os.path.join(huc_directory, filename)+'.shp'

    def register_raster(key, filename):
        filenames[key] = os.path.join(huc_directory, filename)+'.'+raster_extension

    # the HUC shapefile, in native and projected CRS
    register_shapefile('huc', f'huc_{huc}')
    register_shapefile('huc_proj', f'huc_{huc}_proj')

    # the HUC DEM, slope, and aspect rasters
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
    register_shapefile('streams_network_proj', f'huc_{huc}_streams_network_proj') # projected form of above

    # delineated subcatchments within the HUC
    register_shapefile('subcatchments', f'huc_{huc}_subcatchments')  # delineated subcatchments
    
    filenames['flowpaths'] = os.path.join(huc_directory, 'flowpaths', 'hs_{}_flowpaths.pkl')
    register_raster('flowpath_length', f'huc_{huc}_flowpath_lengths') # flowpaths within the delineated subcatchments
    register_raster('elev_above_streams', f'huc_{huc}_elev_above_streams')

    filenames['mesh'] = os.path.join(package_directory, 'meshes', f'huc_{huc}_subcatchment_{{}}.exo')
    filenames['daymet'] = os.path.join(package_directory, 'daymet', f'huc_{huc}_subcatchment_{{}}.h5')
    
    return filenames

def get_huc_and_dem(huc, filenames, huc_source=None, dem_source=None, clobber=False):
    """Downloads HUC and DEM files, saving to disk."""

    # get sources
    if huc_source is None or dem_source is None:
        sources = workflow.source_list.get_default_sources()
        if huc_source is None:
            huc_source = sources['HUC']
        if dem_source is None:
            dem_source = sources['DEM']

    # get shape, from file or from watershed workflow
    if not clobber and os.path.isfile(filenames['huc']):
        huc_crs, huc_shape = workflow.get_shape(filenames['huc'])
    else:
        # download (if necessary) the HUC shapefile
        huc_crs, huc_shape = workflow.get_huc(huc_source, huc)

    # get dem, from file or from watershed workflow
    if not clobber and os.path.isfile(filenames['dem']):
        dem_profile, dem = workflow.get_raster_on_shape(filenames['dem'], huc_shape, huc_crs,
                                                        mask=True)
    else:
        # download (if necessary) the DEM
        dem_profile, dem = workflow.get_raster_on_shape(data_sources['dem'], huc_shape, huc_crs,
                                                        mask=True)

    dem[dem == dem_profile['nodata']] = np.nan
    native_crs = workflow.crs.from_rasterio(dem_profile['crs'])

    # project the shapefile into the native CRS
    huc_shape = workflow.warp.shply(huc_shape, huc_crs, native_crs)

    # save shapefile to disk
    if clobber or not os.path.isfile(filenames['huc']):
        # -- set up the schema -- NOTE: really need to move these to workflow...
        schema = dict()
        if type(shps[0]) is shapely.geometry.Polygon:
            schema['geometry'] = 'Polygon'
        elif type(shps[0]) is shapely.geometry.LineString:
            schema['geometry'] = 'LineString'
        else:
            raise RuntimeError('Currently this function only writes Polygon or LineString types')
        schema['properties'] = collections.OrderedDict()

        # -- set up the properties schema
        def register_type(key, atype):
            if atype is int:
                schema['properties'][key] = 'int'
            elif atype is str:
                schema['properties'][key] = 'str'
            elif atype is float:
                schema['properties'][key] = 'float'
            else:
                pass

        try:
            item_properties = shps[0].properties
        except AttributeError:
            pass
        else:
            for key, val in item_properties.items():
                register_type(key, type(val))

        # -- write to disk
        with fiona.open(filename, 'w', 
                    driver='ESRI Shapefile', 
                    crs=workflow.crs.to_fiona(native_crs), 
                    crs_wkt=workflow.crs.to_wkt(native_crs),
                    schema=schema) as c:
            try:
                props.update(huc_shape.properties)
            except AttributeError:
                pass
            
            for key in list(props.keys()):
                if key not in schema['properties']:
                    props.pop(key)
                  
            c.write({'geometry': shapely.geometry.mapping(huc_shape),
                     'properties': props })

    # save DEM to disk
    if clobber or not os.path.isfile(filenames['dem']):
        # write the raster
        rio_profile = dict(dem_profile).copy()
        rio_profile.pop('blockxsize')
        rio_profile.pop('blockysize')
        rio_profile.pop('tiled')
        rio_profile['nodata'] = -9999.
        rio_profile['driver'] = 'GTiff'

        rio_dem = np.where(np.isnan(dem), rio_profile['nodata'], dem).astype(rio_profile['dtype'])
        with rasterio.open(filenames['dem'], 'w', **rio_profile) as dst:
            dst.write(rio_dem,1)

    return dem_profile, dem, huc_shape


def assertJonIsDone(filenames):
    """A simple function that checks that all filenames needed from Jon are present."""
    # NOTE: this needs to be added by Jon!
    #
    # At this point, we need:
    assert(os.path.isfile(filenames['subcatchments'])) # subcatchment shapefile
    assert(os.path.isfile(filenames['streams_raster'])) # streams raster
    assert(os.path.isfile(filenames['aspect'])) # aspect raster
    assert(os.path.isfile(filenames['slope'])) # slope raster
    assert(os.path.isfile(filenames['flowpath_length'])) # raster of each pixel's distance to the stream
    assert(os.path.isfile(filenames['elev_above_streams'])) # raster of HAND

    
def meanAspect(aspect):
    """Aspects are hard to take the mean of because of the angular aspect.

    # tests...
    assert(meanAspect(np.array([90,90,90])) == 90)
    assert(meanAspect(np.array([0,0,0])) == 0)
    assert(meanAspect(np.array([180,180,180])) == 180)
    assert(meanAspect(np.array([270,270,270])) == 270)
    assert(meanAspect(np.array([89,90,91])) == 90)
    assert(meanAspect(np.array([179,180,181])) == 180)
    assert(meanAspect(np.array([269,270,271])) == 270)
    assert(meanAspect(np.array([359,0,1])) == 0)
    """
    a = aspect[~np.isnan(aspect)]
    a = np.where(a > 180, a - 360, a)
    sina = np.sin(np.radians(a))
    cosa = np.cos(np.radians(a))
    avg_aspect_radians = np.arctan2(sum(sina), sum(cosa))
    if avg_aspect_radians < 0:
        avg_aspect_radians = 2*np.pi + avg_aspect_radians
    return np.degrees(avg_aspect_radians)

def loadSubcatchmentShape(filename, subcatch_id, crs=None):
    """Loads a shapefile to get the subcatchment."""
    crs_new, subcatchments = workflow.get_shapes(filename, crs=crs)
    matches = [sc for sc in subcatchments if sc.properties['hs_id'] == subcatch_id]
    if len(matches) != 1:
        raise RuntimeError("Found invalid subcatchments file or no subcatchment of this id.")
    return crs_new, matches[0]

def loadSubcatchmentRaster(filename, subcatch, subcatch_crs, nanit=True):
    """Load a raster on the subcatchment."""    
    profile, raster = workflow.get_raster_on_shape(filename, subcatch, subcatch_crs, mask=True)
    if nanit:
        raster[raster==profile['nodata']] = np.nan
    return profile, raster

def numSubcatchments(filenames):
    """Returns the number of subcatchments found in the file."""
    _, subcatches = workflow.get_shapes(filenames['subcatchments'])
    return len(subcatches)

def parameterizeSubcatchment(huc, sc, filenames,
                             target_crs=None,
                             hillslope_keep_fraction=0.95,
                             hillslope_bin_dx=100,
                             ):
    """Determine all hillslope properties for a given subcatchment."""

    if target_crs is None:
        target_crs = workflow.crs.default_alaska_crs()

    # get shape
    subcatch_crs, subcatch = loadSubcatchmentShape(filenames['subcatchments'], sc)


    # create a dictionary for hillslope parameters
    hillslope = dict(huc=huc, subcatchment_id=sc)

    # area in m^2
    hillslope['total_area'] = workflow.warp.shply(subcatch, subcatch_crs, target_crs).area

    # centroid in lat/long
    hillslope['centroid'] = workflow.warp.shply(subcatch, subcatch_crs, workflow.crs.latlon_crs()).centroid.coords[0]
    
    # Most of the parameters are based on bins in length of flowpath to the stream network
    # -- load raster of flowpath lengths for a single subcatchment
    fp_profile, fp_lengths = loadSubcatchmentRaster(filenames['flowpath_length'], subcatch, subcatch_crs)
    subcatch_mask = ~np.isnan(fp_lengths)

    # -- ensure >= 0, sort
    assert(fp_lengths[subcatch_mask].min() >= -1.e-3)
    fp_lengths[fp_lengths < 0] = 0.
    fp_lengths_masked = fp_lengths[subcatch_mask]
    fp_lengths_masked_sorted = np.sort(fp_lengths_masked)
    
    # -- set hillslope geometry parameters that are based on flowpath length
    # skip back by hillslope_keep_fraction
    longest_pixel = int(np.round(len(fp_lengths_masked_sorted) * hillslope_keep_fraction))
    hillslope['total_length'] = fp_lengths_masked_sorted[longest_pixel]
    hillslope['num_bins'] = int(np.round(hillslope['total_length'] / hillslope_bin_dx))
    hillslope['bin_dx'] = hillslope['total_length'] / hillslope['num_bins']
    assert(hillslope['num_bins'] > 3)
    
    # -- bin by flowpath length
    pixel_bin_raster = -1 * np.ones(fp_lengths.shape, 'i')
    pixel_bin_counts = np.zeros((hillslope['num_bins'],),'i')

    for i in range(hillslope['num_bins']):
        bin_start = i * hillslope['bin_dx']
        bin_end = (i+1) * hillslope['bin_dx']
        local_mask = (fp_lengths >= bin_start) & (fp_lengths < bin_end)
        pixel_bin_raster[local_mask] = i
        pixel_bin_counts[i] = np.count_nonzero(local_mask)

    assert(pixel_bin_counts.min() > 0)
    hillslope['bin_counts'] = pixel_bin_counts

    # mean height over stream provides elev for each bin
    _, elevs = loadSubcatchmentRaster(filenames['elev_above_streams'], subcatch, subcatch_crs)
    elev_bins = np.zeros((hillslope['num_bins'],),'d')
    for i in range(hillslope['num_bins']):
        elev_bins[i] = elevs[(pixel_bin_raster == i) & (~np.isnan(elevs))].mean()

    hillslope['elevation'] = elev_bins
    
    # average aspect across the entire subcatchment
    _, aspects = loadSubcatchmentRaster(filenames['aspect'], subcatch, subcatch_crs)
    hillslope['aspect'] = meanAspect(aspects)

    # land cover, binned
    lc_profile, lc = loadSubcatchmentRaster(filenames['land_cover'], subcatch, subcatch_crs, False) # don't nan-it
    assert(lc_profile['nodata'] == 255) # uint8, nan is -1 == 255
    lc_bins = np.zeros((hillslope['num_bins'],),'i')

    # -- classify the land cover into veg classes
    land_cover.classifyVegetation(lc)    
    
    for i in range(hillslope['num_bins']):
        pixel_vals = lc[pixel_bin_raster == i]
        assert(len(pixel_vals) > 0)
        lc_bins[i] = scipy.stats.mode(pixel_vals.ravel())[0][0]

    hillslope['land_cover'] = lc_bins
    print('hillslope land cover = ', lc_bins)
    return hillslope


def parameterizeMesh(hillslope, dx,
                     riparian_slope_min=0.01,
                     hillslope_slope_min=0.1,
                     min_area_ratio=0.1
                     ):
    """Given the hillslope parameters, calculate parameters that define the discrete mesh."""
    mesh = dict()
    
    # x-coordinate
    mesh['huc'] = hillslope['huc']
    mesh['subcatchment_id'] = hillslope['subcatchment_id']
    mesh['dx'] = dx
    mesh['num_cells'] = int(np.round(hillslope['total_length'] / dx))
    mesh['length'] = mesh['dx'] * mesh['num_cells']
    x = np.linspace(0, mesh['length'], mesh['num_cells'])
    mesh['x'] = x

    # z-coordinate
    # -- interpolate from z=0 at x=0, and z=bin_average at x=bin_centroid
    x_bin = np.concatenate([np.array([0.,]), (np.arange(0,hillslope['num_bins']) + 0.5)*hillslope['bin_dx']])
    z_bin = np.concatenate([np.array([0.,]), hillslope['elevation']])
    z_native = np.interp(x, x_bin, z_bin)
    z = scipy.signal.savgol_filter(z_native, 11, 3, mode='nearest')
    z = z - z[0]
    mesh['z'] = z

    # -- set minimum slope and determine what is riparian and what is hillslope
    riparian = np.zeros(len(x)-1, 'bool')
    i = 1
    slope = (z[i] - z[i-1])/(x[i] - x[i-1])
    while i < len(z) and slope < 0.1: # in the riparian zone
        riparian[i-1] = True
        if slope < riparian_slope_min:
            z[i] = riparian_slope_min * (x[i] - x[i-1]) + z[i-1]
        i += 1
        if i < len(z):
            slope = (z[i] - z[i-1])/(x[i] - x[i-1])

    while i < len(z): # in the hillslope zone
        riparian[i-1] = False
        if slope < hillslope_slope_min:
            z[i] = hillslope_slope_min * (x[i] - x[i-1]) + z[i-1]
        i += 1
        if i < len(z):
            slope = (z[i] - z[i-1])/(x[i] - x[i-1])

    # -- smooth once more to deal with discontinuities, but only a little...
    z = scipy.signal.savgol_filter(z, 5, 3)
    z = z - z[0]
    assert((z[1:] - z[:-1]).min() > 0)

    # y-coordinate
    # -- interpolate bins to create a bin-consistent profile
    y = np.interp(x, x_bin[1:], hillslope['bin_counts'])

    # -- smooth
    y = scipy.signal.savgol_filter(y, 51, 3, mode='nearest') # window size, poly order

    # -- don't let individual areas get too small -- 10% mean as a min value?
    min_y = min_area_ratio * y.mean()
    y = np.maximum(y, min_y)

    # -- scale by area ratios to ensure that the final mesh has the identical surface area as the
    #    subcatchment it represents
    y_factor = hillslope['total_area'] / np.trapz(y, x)
    y = y_factor * y
    mesh['width'] = y

    # resample land cover onto mesh
    x_bin_nodes = hillslope['bin_dx'] * np.array(range(0,hillslope['num_bins']+1))
    land_cover_mesh = np.zeros(len(x)-1, 'i')
    def lc_index(lc_type, is_riparian):
        if is_riparian:
            return lc_type + 10
        else:
            return lc_type
    for i in range(len(land_cover_mesh)):
        x = (mesh['x'][i] + mesh['x'][i+1])/2
        j_bin = int(np.round(x / hillslope['bin_dx'] - 0.5))
        land_cover_mesh[i] = lc_index(hillslope['land_cover'][j_bin], riparian[i])
    mesh['land_cover'] = land_cover_mesh

    # aspect is still needed
    mesh['aspect'] = hillslope['aspect']
    return mesh


def createMesh2D(mesh_pars):
    """Take mesh parameters and turn those into a 2D surface transect mesh."""
    labeled_sets = list()
    
    for i,vtype in zip(range(100, 104), land_cover.vegClasses()):
        labeled_sets.append(workflow.mesh.LabeledSet(f'hillslope {vtype}', i, 'CELL', 
                                                 [int(c) for c in np.where(mesh_pars['land_cover'] == i)[0]]))
    for i,vtype in zip(range(110, 114), land_cover.vegClasses()):
        labeled_sets.append(workflow.mesh.LabeledSet(f'riparian {vtype}', i, 'CELL', 
                                                 [int(c) for c in np.where(mesh_pars['land_cover'] == i)[0]]))
    
    assert(min(mesh_pars['width']) > 0)
    m2 = workflow.mesh.Mesh2D.from_Transect(mesh_pars['x'], mesh_pars['z'], mesh_pars['width'],
                                            labeled_sets=labeled_sets)

    rotation = 90 + mesh_pars['aspect'] # 90 makes the aspect due north
    m2.transform(mat=workflow.mesh.transform_rotation(np.radians(rotation)))
    return m2

def layeringStructure():
    """For now this is hard-coded, becuase it required significant experimentation."""
    # preparing layer extrusion data -- true for all transects
    #
    # Meshes are extruded in the vertical by "layer", where a layer may 
    # consist of multiple cells in the z direction.  These layers are 
    # logical unit to make construction easier, and may or may not 
    # correspond to material type (organic/mineral soil).
    # 
    # The extrusion process is then given four lists, each of length
    # num_layers.
    #
    layer_types = []  # a list of strings that tell the extruding 
                      # code how to do the layers.  See meshing_ats 
                      # documentation for more, but here we will use
                      # only "constant", which means that dz within
                      # the layer is constant.

    layer_data = []   # this data depends upon the layer type, but
                      # for constant is the thickness of the layer

    layer_ncells = [] # number of cells (in the vertical) in the layer.
                      # The dz of each cell is the layer thickness / number of cells.

    # 30 layers at 2cm makes 60cm of cells, covering the deepest organic layer and getting into mineral soil
    n_top = 30
    dz = 0.02
    current_depth = 0
    for i in range(n_top):
        layer_types.append('constant')
        layer_data.append(dz)
        layer_ncells.append(1)

    # telesecop from 2cm to 2m, in 9.4m, with 20 cells
    dzs, res = workflow.mesh.optimize_dzs(0.02, 2.0, 9.4, 20)

    for dz in dzs:
        layer_types.append('constant')
        layer_data.append(dz)
        layer_ncells.append(1)

    num_at_2m = int(np.round((45 - sum(layer_data)) / 2.0))
    for i in range(num_at_2m):
        layer_types.append('constant')
        layer_data.append(2)
        layer_ncells.append(1)

    return layer_types, layer_data, layer_ncells


#
# Construct the 3D mesh
#
def createMesh3D(m2, mesh_pars, layer_info, fname, clobber=False):
    """Create the 3D mesh via extrusion"""
    if os.path.isfile(fname) and clobber:
        os.remove(fname)

    if not os.path.isfile(fname) or clobber:
        layer_types, layer_data, layer_ncells = layer_info

        layer_mat_ids_near_surface = np.array([land_cover.soilStructure(layer_data, mesh_pars['land_cover'][c]) 
                                               for c in range(m2.num_cells())]).transpose()
        layer_mat_ids = list(layer_mat_ids_near_surface)
    
        # make the mesh, save it as an exodus file
        m3 = workflow.mesh.Mesh3D.extruded_Mesh2D(m2, layer_types, layer_data, layer_ncells, layer_mat_ids)
        m3.write_exodus(fname)

#
# column mesh for spinup
#
def createColumnMesh(layer_info, soil_horizons, filename):
    """Creates a columnar mesh with a given layering profile and soil horizons."""
    m2 = workflow.mesh.Mesh2D.from_Transect(np.array([-0.5, 0.5]), np.array([0., 0.]), 1.0)

    layer_types, layer_data, layer_ncells = layer_info
    layer_mat_ids = list(land_cover.soilStructure(layer_data, -999))

    m3 = workflow.mesh.Mesh3D.extruded_Mesh2D(m2, layer_types, layer_data, layer_ncells, layer_mat_ids)
    m3.write_exodus(filename)

    

#
# Download DayMet
#
def downloadMeteorologyDayMet(hillslope_pars, raw_directory, filename_out):
    """Downloads met data from DayMet"""
    import daymet_to_ats
    start, end = daymet_to_ats.validate_start_end(None, None)
    lon, lat = hillslope_pars['centroid']
    daymet = daymet_to_ats.download_daymet(raw_directory, lat, lon, start, end)

    ats = daymet_to_ats.daymet_to_ats(daymet)
    attrs = daymet_to_ats.daymet_attrs(lat, lon, start, end)        

    daymet_to_ats.write_ats(ats, attrs, filename_out)
    

if __name__ == "__main__":
    huc = '190604020404'
    mesh_dx = 20
    package_directory = '/Users/uec/research/interface/problems/interface-hillslopes/'
    filenames = get_filenames(huc, package_directory)

    # set up directory structure for raw data
    workflow.ui.setup_logging(0)
    raw_data_dir = os.path.join(package_directory, 'data_preprocessed-meshing')
    if not os.path.isdir(raw_data_dir):
        os.mkdir(raw_data_dir)

    # set up directory structure for meshes
    mesh_dir = os.path.join(package_directory, 'meshes')
    if not os.path.isdir(mesh_dir):
        os.mkdir(mesh_dir)

    # column mesh for spinup
    column_mesh_filename = os.path.join(package_directory, 'meshes', 'column.exo')
    if not os.path.isfile(column_mesh_filename):
        layering = layeringStructure()
        createColumnMesh(layering, land_cover.soil_horizons['average'], column_mesh_filename)

    # set up directory structure for forcing datasets
    daymet_dir = os.path.join(package_directory, 'daymet')
    if not os.path.isdir(daymet_dir):
        os.mkdir(daymet_dir)

    # column forcing
    column_met_forcing_filename = os.path.join(package_directory, 'daymet', f'huc_{huc}.h5')
    if not os.path.isfile(column_met_forcing_filename):
        _, huc_shape = workflow.get_shapes(filenames['huc'], crs=workflow.crs.latlon_crs())
        centroid = huc_shape[0].centroid.coords[0]
        downloadMeteorologyDayMet(dict(centroid=centroid), raw_data_dir, column_met_forcing_filename)

    # setup
    workflow.conf.rcParams['DEFAULT']['data_directory'] = raw_data_dir
    assertJonIsDone(filenames)
    alaska_crs = workflow.crs.default_alaska_crs()

    for sc in range(1,numSubcatchments(filenames)+1):
        hillslope_pars = parameterizeSubcatchment(huc, sc, filenames, alaska_crs)
        raise RuntimeError()

        # make the mesh
        mesh_filename = filenames['mesh'].format(sc)
        if not os.path.isfile(mesh_filename):
            print(f'Generating mesh: {mesh_filename}')
            mesh_pars = parameterizeMesh(hillslope_pars, mesh_dx)
            m2 = createMesh2D(mesh_pars)
            layering = layeringStructure()
            createMesh3D(m2, mesh_pars, layering, filenames['mesh'].format(sc), clobber=True)

        # download daymet
        daymet_filename = filenames['daymet'].format(sc)
        if not os.path.isfile(daymet_filename):
            print(f'Download meteorology: {daymet_filename}')
            downloadMeteorologyDayMet(hillslope_pars, raw_data_dir, daymet_filename)
