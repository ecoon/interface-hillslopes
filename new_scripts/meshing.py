import numpy as np
import scipy
import hillslopes
import landcover
import workflow
import os

# Given a hillslope's parameters, generate the mesh parameters
def parameterizeMesh(hillslope, dx,
                     riparian_slope_min=0.01,
                     hillslope_slope_min=0.1,
                     min_area_ratio=0.1):
    mesh = dict()
    
    # 1. Determine x coordinate
    mesh['huc'] = hillslope['huc']
    mesh['subcatchment_id'] = hillslope['subcatchment_id']
    mesh['dx'] = dx
    mesh['num_cells'] = int(np.round(hillslope['total_length'] / dx))
    mesh['length'] = mesh['dx'] * mesh['num_cells']
    mesh['x'] = np.linspace(0, mesh['length'], mesh['num_cells'])
    
    # 2. Determine z coordinate
    # 2.1 Interpolate from (x,z)=(0,0), and (x,z)=(bin_centroid,bin_average)
    x_bin = np.concatenate([np.array([0.,]), (np.arange(0,hillslope['num_bins']) + 0.5)*hillslope['bin_dx']])
    z_bin = np.concatenate([np.array([0.,]), hillslope['elevation']])
    z_native = np.interp(mesh['x'], x_bin, z_bin)
    z = scipy.signal.savgol_filter(z_native, window_length=11, polyorder=3)
    z = z - z[0]
    
    # 2.2 Determine riparian area and hillslope area according to slope
    slope = (z[1:] - z[0:-1]) / (mesh['x'][1:] - mesh['x'][0:-1])
    riparian = np.zeros(slope.shape, 'i')
    i = 0
    while i < len(slope) and slope[i] < 0.1:     # in the riparian zone
        riparian[i] = 1
        if slope[i] < riparian_slope_min:
            z[i+1] = riparian_slope_min * (mesh['x'][i+1] - mesh['x'][i]) + z[i]
        i += 1
        
    if i == 1:
        mesh['riparian_width'] = 0
    else:
        mesh['riparian_width'] = (mesh['x'][i-1] + mesh['x'][i])/2.
    
    while i < len(slope):   # hillslope starts from here 
        if slope[i] < hillslope_slope_min:
            z[i+1] = hillslope_slope_min * (mesh['x'][i+1] - mesh['x'][i]) + z[i]
        i += 1
        
    mesh['riparian'] = riparian
    
    # smooth z once more to deal with discontinuities
    z = scipy.signal.savgol_filter(z, window_length=5, polyorder=3)
    z = z - z[0]
    mesh['z'] = z
    
    # 3. Determine y coordinate
    # 3.1 Interpolate bins to create a bin-consistent profile
    y = np.interp(mesh['x'], x_bin[1:], hillslope['bin_counts'])
    
    # 3.2 Smooth y
    y = scipy.signal.savgol_filter(y, window_length=51, polyorder=3, mode='nearest')
    
    # 3.3 Don't let individual areas get too small -- 10% mean as a min value?
    min_y = min_area_ratio * y.mean()
    y = np.maximum(y, min_y)
    
    # 3.4 Scale by area ratios to ensure that the final mesh has the identical surface area as the subcatchment it represents
    y_factor = hillslope['total_area'] / np.trapz(y, mesh['x'])
    y = y_factor * y
    mesh['y'] = y
    
    # 4. Resample land cover onto mesh
    land_cover_mesh = np.zeros(slope.shape, 'i')
    def lc_index(lc_type, is_riparian):
        if is_riparian:
            return lc_type + 10
        else:
            return lc_type
        
    for i in range(len(land_cover_mesh)):
        x = (mesh['x'][i] + mesh['x'][i+1])/2
        j_bin = np.where(x_bin>x,x_bin,np.inf).argmin()-1
#         j_bin = int(np.round(x / hillslope['bin_dx'] - 0.5))   # THE OLD
        land_cover_mesh[i] = lc_index(hillslope['land_cover'][j_bin], riparian[i])
    mesh['land_cover'] = land_cover_mesh
    
    # 5. Add aspect to mesh
    mesh['aspect'] = hillslope['aspect']
    
    return mesh
                
                

# Take mesh parameters and turn those into a 2D surface transect mesh
def createHillslopeMesh2D(mesh_pars):
    labeled_sets = list()
    for i,vtype in zip(range(100, 104), landcover.vegClasses()):
        labeled_sets.append(workflow.mesh.LabeledSet(f'hillslope {vtype}', i, 'CELL', 
                                                     [int(c) for c in np.where(mesh_pars['land_cover'] == i)[0]]))
    for i,vtype in zip(range(110, 114), landcover.vegClasses()):
        labeled_sets.append(workflow.mesh.LabeledSet(f'riparian {vtype}', i, 'CELL',
                                                     [int(c) for c in np.where(mesh_pars['land_cover'] == i)[0]]))
    assert(min(mesh_pars['y']) > 0)
    m2 = workflow.mesh.Mesh2D.from_Transect(mesh_pars['x'], mesh_pars['z'], mesh_pars['y'],
                                            labeled_sets=labeled_sets)
    rotation = 90 + mesh_pars['aspect'] # 90 makes the aspect due north
    m2.transform(mat=workflow.mesh.transform_rotation(np.radians(rotation)))
    return m2
    

# Preparing layer extrusion data
#
# Meshes are extruded in the vertical by "layer", where a layer may consist of multiple cells in the z direction.
# These layers are logical unit to make construction easier, 
# and may or may not correspond to material type (organic/mineral soil).
#
# The extrusion process is then given four lists, each of length num_layers.
def layeringStructure(organic_cells=30, organic_cell_dz=0.02, 
                      increase2depth=9.4, increase_cells=20, largest_dz=2.0,
                      bottom_depth=45):
    layer_types = []  # a list of strings that tell the extruding code how to do the layers.  
                      # See meshing_ats documentation for more, but here we will use only "constant",
                      # which means that dz within the layer is constant.

    layer_data = []   # this data depends upon the layer type, but for constant is the thickness of the layer
    
    layer_ncells = [] # number of cells (in the vertical) in the layer.
                      # The dz of each cell is the layer thickness / number of cells.
        
    # 30 layers at 2cm makes 60cm of cells, covering the deepest organic layer and getting into mineral soil
    n_top = organic_cells
    dz = organic_cell_dz
    current_depth = 0
    for i in range(n_top):
        layer_types.append('constant')
        layer_data.append(dz)
        layer_ncells.append(1)
        
    # telescope from 2cm to 2m
    dzs, res = workflow.mesh.optimize_dzs(organic_cell_dz, largest_dz, increase2depth, increase_cells)
    for dz in dzs:
        layer_types.append('constant')
        layer_data.append(dz)
        layer_ncells.append(1)
        
    num_at_2m = int(np.round((bottom_depth - sum(layer_data)) / largest_dz))
    for i in range(num_at_2m):
        layer_types.append('constant')
        layer_data.append(2)
        layer_ncells.append(1)
        
    return layer_types, layer_data, layer_ncells
        
    
    
# Exstrude the 2D hillslope surface mesh downward to construct 3D hillslope mesh
def createHillslopeMesh3D(m2, mesh_pars, layer_info, mesh_fname):
    if os.path.isfile(mesh_fname):
        os.remove(mesh_fname)
        
    layer_types, layer_data, layer_ncells = layer_info
    layer_mat_ids_near_surface = np.array([landcover.soilStructure(layer_data, mesh_pars['land_cover'][c]) 
                                          for c in range(m2.num_cells())]).transpose()
    layer_mat_ids = list(layer_mat_ids_near_surface)
    m3 = workflow.mesh.Mesh3D.extruded_Mesh2D(m2, layer_types, layer_data, layer_ncells, layer_mat_ids)
    m3.write_exodus(mesh_fname)
        
    

# Create column mesh for spinup
def createColumnMesh(layer_info, filename):
    m2 = workflow.mesh.Mesh2D.from_Transect(np.array([-0.5, 0.5]), np.array([0., 0.]), 1.0)

    layer_types, layer_data, layer_ncells = layer_info
    layer_mat_ids = list(landcover.soilStructure(layer_data, -999))

    m3 = workflow.mesh.Mesh3D.extruded_Mesh2D(m2, layer_types, layer_data, layer_ncells, layer_mat_ids)
    m3.write_exodus(filename)    
    
    
    

    
    
    