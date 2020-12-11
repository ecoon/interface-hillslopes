"""Namespace for dealing with NSSI Land Cover dataset."""

import attr
import numpy as np
import rasterio

#
# Interpretation for land cover map
#
_lc_ids = [1,2,3,5,6,7,8,9,10,11,14,16,18,19,20,21,22,23,50,51,91,92,95,96,99]
_lc_names = ['Bare Ground','Sparsely Vegetated','Open Water','FWM: Arctophila Fulva',
            'FWM: Carex Aquatillis','Wet Sedge','Wet Sedge - Sphagnum','Mesic Herbaceous','Tussock Tundra',
            'Tussock Shrub Tundra','Mesic Sedge-Dwarf Shrub Tundra','Dwarf Shrub - Dryas','Dwarf Shrub - Other',
            'Birch Ericaceous Low Shrub','Low-Tall Willow','Alder','Marine Beach/Beach Meadow','Coastal Marsh',
            'Ice / Snow','Burned Area','Open Needleleaf','Woodland Needleleaf','Open Mixed Needleleaf/Deciduous',
            'Deciduous','Unclassified (cloud, terrain shadow, etc.)']

# maps from name to id and back
name_to_id = dict(zip(_lc_names, _lc_ids))
id_to_name = dict(zip(_lc_ids, _lc_names))

# maps from coarse-binned class name to a list of names
class_to_names = dict(
    sparse_veg=['Bare Ground','Sparsely Vegetated','Dwarf Shrub - Other','Open Water','Ice / Snow'],
    sedge=['FWM: Carex Aquatillis', 'Mesic Sedge-Dwarf Shrub Tundra', 'Wet Sedge','Wet Sedge - Sphagnum'],
    shrub=['Birch Ericaceous Low Shrub','Low-Tall Willow','Alder'],
    tussock=['Tussock Tundra','Tussock Shrub Tundra', 'Dwarf Shrub - Dryas']
    )

# maps from coarse-binned class name to a list of ids
class_to_ids = dict((k, [name_to_id[name] for name in val]) for k,val in class_to_names.items())

# all ids currently handled
all_ids = [v for val in class_to_ids.values() for v in val]

# soil structure
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

def vegClasses():
    return list(class_to_ids.keys())


def classifyVegetation(lc):
    """Remaps from NSSI IDs to classes: sparse_veg, sedge, shrub, and tussock.  INPLACE"""
    lc_vals = set(lc.ravel())
    if 255 in lc_vals:
        lc_vals.remove(255) # nan

    for i,(vclass,ids) in enumerate(class_to_ids.items()):
        ats_id = 100 + i
        for veg_id in ids:
            try:
                lc_vals.remove(veg_id)
            except KeyError:
                pass
            lc[lc == veg_id] = ats_id
        
    if not len(lc_vals) == 0:
        raise RuntimeError(f'Unclassified vegetation ids {lc_vals} type {[lc_keys_reversed[m] for m in lc_vals]}')

def vegToImage(lc):
    """Returns a float version with NaNs to improve plotting"""
    new_lc = np.nan * np.ones(lc.shape, 'd')
    for i in set(lc.ravel()):
        new_lc[lc == i] = i
    new_lc[lc == 255] = np.nan
    return new_lc
    
def reprojectLandCover(dem_profile, lc_filename, clobber=False):
    """Reprojects land cover from the NSSI default CRS into (ugh) Lat/Long the CRS of the bands"""

    if clobber or not os.path.isfile(lc_filename):
        rio_profile = dem_profile.copy()
        rio_profile.pop('blockxsize')
        rio_profile.pop('blockysize')
        rio_profile.pop('tiled')
        rio_profile['nodata'] = 255
        rio_profile['driver'] = 'GTiff'
    
        with rasterio.open('../data/6979-north-slope-landcover-map/nssi_landcover_final_2013oct1_albers83.img', 'r') as fin:
            with rasterio.open(lc_filename, 'w', **rio_profile) as fout:
                rasterio.warp.reproject(
                    source=rasterio.band(fin, 1),
                    destination=rasterio.band(fout, 1),
                    src_transform=fin.transform,
                    src_crs=fin.crs,
                    dst_transform=rio_profile['transform'],
                    dst_crs=rio_profile['crs'],
                    resampling=rasterio.warp.Resampling.nearest)
      

def soilStructure(dz, land_cover):
    """Given a cell thickness profile dz, determine the label of each cell
    
    1000 = mineral_soil
    1001 = acrotelm
    1002 = catotelm
    """
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

