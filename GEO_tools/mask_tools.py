########################################
# Tools for creating geogrphic masks
#


############################################################################################
# Function: get_geometries
# Purpose:  return cartopy geometries for specified countries
# 
# Inputs: lat,lon - country_names - list of strings of country names
#                                Must corresponde to cartopy naming conventions
#
# Outputs: geoms - list of geometries
#              
############################################################################################    
def get_geometries(country_names):
    import cartopy.crs as ccrs
    from cartopy.io.shapereader import natural_earth, Reader

    """
    Get an iterable of Shapely geometries corrresponding to given countries.

    """
    # Using the Natural Earth feature interface provided by cartopy.
    # You could use a different source, all you need is the geometries.
    shape_records = Reader(natural_earth(resolution='110m',
                                         category='cultural',
                                         name='admin_0_s')).records()
    geoms = []
    for country in shape_records:
        if country.attributes['name_long'] in country_names:
            print  country.attributes['name_long']
            try:
                geoms += country.geometry
            except TypeError:
                geoms.append(country.geometry)
    transform = ccrs.PlateCarree()._as_mpl_transform

    return geoms, transform



############################################################################################
# Function: country_mask
# Purpose:  return a country mask for a given lat lon grid
# 
# Inputs: lat,lon - latitude and longitude grid. 
#                     Must have identical shapes
#
# Outputs: Mask - A boolean array of same shape of lat,lon.
#                    1 is where the pixel is "in" the specified country
#
# Optional Inputs: Country='land' - Country to return mask.#
#                                   Default is all land points.
#                  Fraction=0.5 - Fraction of pixel to be masked as "in" specified country.
#                  
############################################################################################
def country_mask(lat,lon, Country='land',Fraction=0.5):
    
    shape_records = Reader(natural_earth(resolution='110m',
                                         category='cultural',
                                         name='admin_0_countries')).records()
    geoms = []
    for country in shape_records:
        if country.attributes['name_long'] in country_names:
            try:
                geoms += country.geometry
            except TypeError:
                geoms.append(country.geometry)
