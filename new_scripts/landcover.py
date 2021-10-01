import numpy as np
import rasterio
import attr

_lc_ids = [1,2,3,5,6,7,8,9,10,11,14,16,18,19,20,21,22,23,50,51,91,92,95,96,99]
_lc_names = ['Bare Ground','Sparsely Vegetated','Open Water','FWM: Arctophila Fulva',
             'FWM: Carex Aquatillis','Wet Sedge','Wet Sedge - Sphagnum','Mesic Herbaceous','Tussock Tundra',
             'Tussock Shrub Tundra','Mesic Sedge-Dwarf Shrub Tundra','Dwarf Shrub - Dryas','Dwarf Shrub - Other',
             'Birch Ericaceous Low Shrub','Low-Tall Willow','Alder','Marine Beach/Beach Meadow','Coastal Marsh',
             'Ice / Snow','Burned Area','Open Needleleaf','Woodland Needleleaf','Open Mixed Needleleaf/Deciduous',
             'Deciduous','Unclassified (cloud, terrain shadow, etc.)']

name_to_id = dict(zip(_lc_names, _lc_ids))
id_to_name = dict(zip(_lc_ids, _lc_names))

class_names, class_ids = dict(), dict()
class_names['shrub'] = ['Birch Ericaceous Low Shrub','Low-Tall Willow','Alder']
class_names['sedge'] = ['FWM: Carex Aquatillis', 'Mesic Sedge-Dwarf Shrub Tundra', 'Wet Sedge','Wet Sedge - Sphagnum']
class_names['tussock'] = ['Tussock Tundra','Tussock Shrub Tundra', 'Dwarf Shrub - Dryas']
class_names['sparse_veg'] = ['Bare Ground','Sparsely Vegetated','Dwarf Shrub - Other','Open Water','Ice / Snow']
for group, names in class_names.items():
    class_ids[group] = [name_to_id[name] for name in names]

    
def vegClasses():
    return list(class_ids.keys())
    
        
# Remaps from NSSI IDs to lumped ATS IDs around shrub,sedge,tussock,and sparse_veg
def classifyVegetation(lc):
    lc_remap = np.zeros(lc.shape, dtype=np.uint8)
    lc_remap[lc == 255] = 255
    for i, (group,ids) in enumerate(class_ids.items()):
        ats_id = 100+i
        for veg_id in ids:
            lc_remap[lc == veg_id] = ats_id
    missing_id = set(lc[[lc_remap == 0]])
    for mid in missing_id:
        print(f'Missing vegetation id {mid} type: ', id_to_name[mid])
        
    return len(missing_id), lc_remap


# Return a float version with NaNs to improve plotting
def veg2img(lc):
    new_lc = np.nan * np.ones(lc.shape, 'd')
    for i in set(lc.ravel()):
        new_lc[lc == i] = i
    new_lc[lc == 255] = np.nan
    return new_lc
    
    
# Reprojects land cover from the NSSI default CRS into (ugh) Lat/Long the CRS of the bands
def reprojectLandCover(dem_profile, nssiImg_filename, lc_filename):
    rio_profile = dem_profile.copy()
    rio_profile.pop('blockxsize')
    rio_profile.pop('blockysize')
    rio_profile.pop('tiled')
    rio_profile['nodata'] = 255
    rio_profile['driver'] = 'GTiff'
    
    with rasterio.open(nnsiImg_filename, 'r') as fin:
        with rasterio.open(lc_filename, 'w', **rio_profile) as fout:
            rasterio.warp.reproject(
                source=rasterio.band(fin, 1),
                destination=rasterio.band(fout, 1),
                src_transform=fin.transform,
                src_crs=fin.crs,
                dst_transform=rio_profile['transform'],
                dst_crs=rio_profile['crs'],
                resampling=rasterio.enums.Resampling.nearest)
            


    
# Determine soil structure according to land cover
@attr.s
class SoilHorizons:
    acrotelm = attr.ib()
    catotelm = attr.ib()
    
soil_horizons = dict()
soil_horizons['riparian sparse_veg'] = SoilHorizons(acrotelm=0, catotelm=0)
soil_horizons['hillslope sparse_veg'] = SoilHorizons(acrotelm=0, catotelm=0)

soil_horizons['riparian sedge'] = SoilHorizons(acrotelm=0.10, catotelm=0.30)
soil_horizons['hillslope sedge'] = SoilHorizons(acrotelm=0.08, catotelm=0.20)

soil_horizons['riparian shrub'] = SoilHorizons(acrotelm=0.12, catotelm=0.24)
soil_horizons['hillslope shrub'] = SoilHorizons(acrotelm=0.14, catotelm=0.14)

soil_horizons['riparian tussock'] = SoilHorizons(acrotelm=0.10, catotelm=0.14)
soil_horizons['hillslope tussock'] = SoilHorizons(acrotelm=0.10, catotelm=0.14)

soil_horizons['average'] = SoilHorizons(acrotelm=0.1, catotelm=0.2)
    


# Given a soil cell thickness dz, determine the material label of each cell            
def soilStructure(dz, land_cover):
    '''
    1000 = mineral
    1001 = acrotelm
    1002 = catotelm
    '''
    if land_cover == -999:
        lc_name = 'average'
    elif land_cover < 110:
        lc_name = 'hillslope ' + vegClasses()[land_cover - 100]
    else:
        lc_name = 'riparian ' + vegClasses()[land_cover - 110]
        
    horizons = soil_horizons[lc_name]
    z_depth = np.cumsum(dz)
    soil_type = np.where(z_depth < horizons.acrotelm, 1001,
                         np.where(z_depth < horizons.acrotelm + horizons.catotelm, 1002, 1000))
    
    return soil_type