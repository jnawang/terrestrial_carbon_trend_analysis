import ee
import datetime

ee.Authenticate()
ee.Initialize(project='carbon-trend')

# -------- User-defined parameters --------
model = 'MIROC6'
scenario = 'ssp245'
variable = 'rsds'
member_id = 'r1i1p1f1'
start_year = 2050
end_year = 2050
year = 2005


dataset = ee.ImageCollection('NASA/GDDP-CMIP6') \
    .filterDate('2014-07-01', '2014-07-30') \
    .filter(ee.Filter.eq('model', 'ACCESS-CM2'))\
    .sort('system:time_start')

# Select the 'tasmin' band (minimum air temperature)
collection = dataset.select('rsds')

# Check number of images
print('Collection size:', collection.size().getInfo())

    # Convert daily images to multiband image
def convert_to_band(img):
    date_str = ee.Date(img.get('system:time_start')).format('YYYYMMdd')
    return img.rename(date_str)
        
multiband = collection.map(convert_to_band).toBands()

task = ee.batch.Export.image.toDrive(
    image=multiband,
    description=f'{variable}_{model}_{scenario}_{year}_daily_global',
    folder='EarthEngine',
    fileNamePrefix=f'{variable}_{model}_{scenario}_{year}_daily',
    scale=25000,
    region=ee.Geometry.BBox(-180, -90, 180, 90),
    crs='EPSG:4326',
    maxPixels=1e13
)
task.start()
print(f'Started export for {year}')




# -------- Loop over years --------
##for year in range(start_year, end_year + 1):
##    start_date = f'{year}-01-01'
##    end_date = f'{year}-12-31'
##
##    # Filter collection for the year
##    collection = ee.ImageCollection('NASA/GDDP-CMIP6') \
##        .filter(ee.Filter.eq('model', model)) \
##        .filter(ee.Filter.eq('scenario', scenario)) \
##        .filterDate("2050-01-01", "2050-12-31") 
##
##
##
##    print("Size of collection:", collection.size().getInfo())



##    # Export to Drive
##    task = ee.batch.Export.image.toDrive(
##        image=multiband,
##        description=f'{variable}_{model}_{scenario}_{year}_daily_global',
##        folder='EarthEngine',
##        fileNamePrefix=f'{variable}_{model}_{scenario}_{year}_daily',
##        scale=25000,
##        region=ee.Geometry.BBox(-180, -90, 180, 90),
##        crs='EPSG:4326',
##        maxPixels=1e13
##    )
##    task.start()
##    print(f'Started export for {year}')

##        .filter(ee.Filter.eq('member_id', member_id)) \
