# HREA

Replication code to produce High Resolution Electricity Access indicators.


## About

The High Resolution Electricity Access (HREA) project led by University of Michigan Professor Brian Min aims to provide open access time series indicators of electricity access and reliability across the world.  Electrification metrics are computed via statistical methods utilizing multiple geospatial data sets and nighttime satellite imagery.  For details on the methodology, please read our [paper](https://authors.elsevier.com/a/1j8j9925JERZpQ) and the accompanying [supplemental information](https://www.cell.com/cms/10.1016/j.joule.2024.05.001/attachment/609c1be6-e561-4462-947f-eca092193c4b/mmc1.pdf) published in Joule.  For the most up-to-date information, please see the [HREA project website](https://hrea.isr.umich.edu/).


## Description

This GitHub repository contains code to recreate HREA output from publicly available inputs.  HREA leverages nighttime satellite imagery publicly available on AWS through the Light Every Night (LEN), a partnership between the National Oceanic and Atmospheric Administration (NOAA), World Bank, and University of Michigan.  One can learn more about the LEN repository [here](https://registry.opendata.aws/wb-light-every-night/).  The HREA data introduced in "Lost in the Dark: A Survey of Energy Poverty from Space" utilizes orbital strip data from the Suomi National Polar-orbiting Partnership (SNPP) Visible Infrared Imaging Radiometer Suite (VIIRS) Day/Night Band (DNB) and associated metadata.

HREA data were generated using [R](https://www.r-project.org/).  All core functions appear in the order in which they must be run in [R/HREA_functions.R](../blob/master/R/HREA_functions.R).  Comments in this file describe each function, the requisite inputs, and the output each generates.


## Preparation Steps

To generate HREA, one must install the following packages like so:

```
install.packages(c('data.table', 'terra', 'geodata', 'MODISTools', 'stackoverflow', 'lme4'))
```

The geospatial packages may require additional installation and setup to work.  In addition, the multicore functionality of `data.table` may require additional setup for macOS users.

The code available here depends on auxiliary data.  In particular, one must download the following:

* [HREA_vflag_ints.rds](https://globalnightlight.s3.amazonaws.com/HREA_aux_data/HREA_vflag_ints.rds): Vector of VIIRS flags indicating high-quality data.  Specify local filepath with `good_vflag_path`.
* [HREA_place_lookup.csv](https://globalnightlight.s3.amazonaws.com/HREA_aux_data/HREA_place_lookup.csv): Table of places with HREA data containing information about the country (*country\_iso3*), COG identifier (*cog\_id*), and settlement layer (*set\_layer*).  Specify local filepath with `cog_lookup_path`.
* [land_cover_best_alt_type.rds](https://globalnightlight.s3.amazonaws.com/HREA_aux_data/land_cover_best_alt_type.rds): Table of land cover types matched to others ranked by albedo similarity.  Specify local filepath with `LC_sub_path`.

In addition, researchers must download the appropriate settlement layers. The lookup table (HREA\_place\_lookup.csv) describes whether Data for Good at Meta's High Resolution Settlement Layer (HRSL) or the European Commission's Global Human Settlement Layer (GHSL) was used for a particular place.  HRSL data can be obtained from the Humanitarian Data Exchange (HDX) at the [dataset page](https://data.humdata.org/dataset/?dataseries_name=Data+for+Good+at+Meta+-+High+Resolution+Population+Density+Maps+and+Demographic+Estimates), but note that some of these layers have been updated over time and may not exactly match the settlement layers used in Min et al. (2024).  GHS-POP can be downloaded from the [project website](https://human-settlement.emergency.copernicus.eu/download.php?ds=pop).  The GHSL layer employed in Min et al. (2024) is GHS\_POP\_E2015\_GLOBE\_R2019A\_54009\_250\_V1\_0, but this version has been classified as obsolete and does not appear to be available anymore.  One may need to alter the code to work with other versions.  Note too that Meta's high resolution population density maps sometimes cross multiple countries and may also have multiple tiles per country.  Brazil, for example, is provided as four tiles.  We split India into states and territories (using GADM 3.6 definitions).  Fiji was split because it crosses the antimeridian.  Splitting large countries is advisable considering computational constraints and R's vector length maximum of 2<sup>31</sup>-1.


## Processing Steps

Once the user has installed the necessary software and data, one must first prepare static *cog\_id*-level data following these steps:

1. `make_country_boundary`: Creates a cropped and simplified areal boundary polygon.
2. `make_grid_pop_rast`: Creates a settlement-layer level raster grid and fixed 15 arcsecond grid.  Extracts population data from settlement layer, as well as MODIS land cover values.
3. `find_iso_nsets`: Finds non-settlement cells that are at least 15 arcseconds away from cells with any population.

Once overlapping LEN COGs have been identified for a given location and downloaded locally, one can download and clean the data by *year* and *month* for some *cog\_id* with `get_VIIRS_dat`.

After obtaining data for a full *year* and *country_iso3*, one can proceed to generate the annual HREA model with `make_HREA_mod`.  This function reads isolated non-settlement data, excludes outliers, and randomly samples the data.  After some variable manipulation, a linear mixed model (LMM) is run predicting radiance using the [HREA regression model](https://hrea.isr.umich.edu/methods.html).

Having generated the appropriate *country_iso3*-level model, one may then generate HREA statistics for an associated *cog\_id* and *year* for settlements with `make_HREA_data`.  This function reads the settlement data, generates predicted values and standardized residuals from the regression model coefficients, and then computes HREA metrics.

Finally, one can generate HREA Cloud-Optimized GeoTIFFs like those found on the LEN repository for a given *cog_id* and *year* using `make_HREA_cogs`.


## Caveats

This code should run using the most up-to-date versions of the software and packages as of May, 2024.  However, it may need to be adapted as packages are updated or removed.

The code provided here is not parallelized due to differing operating system configurations, but some steps can be done in parallel.  Users are encouraged to tweak the code to fit their environment and goals.  Note that parallelizing many of the operations requires increasingly large amounts of RAM.

Please note that it is not recommended to attempt to recreate HREA on a personal machine.  The HREA process is best suited to run on supercomputing clusters due to the amount of disk space, memory, and CPU time required.


## Data Access

HREA data in the form of Cloud-Optimized GeoTIFFs (COGs) are available for download from Amazon Web Services (AWS) as part of the Light Every Night (LEN) repository.  We use the following variables to refer to the data:

+ *country\_iso3*: ISO ALPHA-3 country code
+ *cog\_id*: specific identifier for country or sub-country boundary
+ *year*: year in YYYY format

We offer four HREA products for download.  These data types are described in the following table:

| Description | File Stem | Resolution | Coverage | Value Range |
| ----------- | --------- | ---------- | -------- | ----------- |
| Mean annual natural logged visible radiance (in nW/cm<sup>2</sup>/sr) + 2.5 following EOG procedure for dropping low quality observations | *rade9lnmu* | 15 arcsec | All cells | [0, 1] |
| Mean *z*-score (standardized residuals) | *set\_zscore* | 1 or 7.5 arcsec | Settlement cells | [-29, 218] |
| Light score (transformation of mean *z*-score) | *set\_lightscore* | 1 or 7.5 arcsec | Settlement cells | [0, 1] |
| Proportion of nights a settlement is brighter than uninhabited areas | *set\_prplit* | 1 or 7.5 arcsec | Settlement cells | [0, 1] |

COGs are stored as 32-bit floats with missing values corresponding to -999999.  All COGs use the World Geodetic System unprojected longitude-latitude coordinate reference system (WGS 84, EPSG:4326).  Each COG asset has a corresponding item json, and searchable SpatioTemporal Asset Catalog (STAC) collections and catalogs organize the data following best practices.

HREA rasters and metadata are stored in sub-directories below the top-level HREAv1.1\_COGs folder, in which 115 directories exist named according to their *country\_iso3* country code. Inside each *country\_iso3* folder are directories for each cog_id belonging to that country.  In most cases, there is only one *cog\_id* which is equivalent to the *country\_iso3*; however, some countries are split into multiple parts out of necessity.  In such cases, the *cog_id* is formatted as *country_iso3-X*, where *X* may correspond to a compass direction indicator (as in Brazil and Fiji) or a two-letter state code as in India.  For each *cog\_id*, there exists a constant settlement population raster file HREAv1.1\_COGs/*country\_iso3*/*cog\_id*/*cog\_id*\_set\_pop.tif, as well as a STAC item with the same name but instead ending with .json.  In addition, for each year in 2013–2020, there exist the following GeoTIFFs, which also have corresponding STAC items with names that substitute .json for the .tif extension:

+ HREAv1.1\_COGs/*country\_iso3*/*country\_iso3*/*cog\_id*\_rade9lnmu\_*year*.tif: Nighttime light annual composite for all cells in a locale.
+ HREAv1.1\_COGs/*country\_iso3*/*country\_iso3*/*cog\_id*\_set\_zscore\_*year*.tif: Statistically estimated brightness levels for settlements.  Higher levels indicate more robust usage of outdoor lighting, which is correlated with overall energy consumption.
+ HREAv1.1\_COGs/*country\_iso3*/*country\_iso3*/*cog\_id*\_set\_lightscore\_*year*.tif: Predicted likelihood that a settlement is electrified (0 to 1).
+ HREAv1.1\_COGs/*country\_iso3*/*country\_iso3*/*cog\_id*\_set\_prplit\_*year*.tif: Proportion of nights a settlement is statistically brighter than matched uninhabited areas.

One can obtain HREA COGs from AWS either by direct download using the object URLs or through the AWS command line interface (CLI).  Using the CLI, one can list the contents of a directory with `aws s3 ls s3://globalnightlight/path/to/dir`.  With the AWS CLI, one must use the `--no-sign-request` flag.  To list the contents of the top-level HREA directory, one would run:

```
aws s3 ls s3://globalnightlight/HREAv1.1_COGs --no-sign-request
```

The root HREA STAC catalog, which lists links to all *country\_iso3* catalogs, can be downloaded with a tools like `wget`, `curl`, or `aws-cli`:

```
wget https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/catalog.json

curl -O https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/catalog.json

aws s3 cp s3://globalnightlight/HREAv1.1_COGs/catalog.json . --no-sign-request
```

One can find the *cog\_id*s belonging to a specific *country\_iso3* by checking the [lookup table](https://globalnightlight.s3.amazonaws.com/HREA_aux_data/HREA_place_lookup.csv) or listing the contents in the folder with:

```
aws s3 ls s3://globalnightlight/HREAv1.1_COGs/${country_iso3} --no-sign-request
```

One could read a specific *country\_iso3*'s catalog to determine linked catalogs and collections:

```
curl https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/${country_iso3}/catalog.json
```

Individual assets can be downloaded using the path structure described above.  As an example, one could obtain the mean standardized residuals (*z*-scores) for Ghanaian settlements in 2020 with:

```
aws s3 cp s3://globalnightlight/HREAv1.1_COGs/GHA/GHA/set_zscore/GHA_set_zscore_2020.tif . --no-sign-request
```

Alternatively, one can download the contents of an entire directory for some *country\_iso3* to some directory `/path/to/dir/` like so:

```
aws s3 cp s3://globalnightlight/HREAv1.1_COGs/${country_iso3}/ /path/to/dir/ --recursive --no-sign-request
```

For an example of how one could traverse the HREA STAC and download items in R, please see [R/download_HREA_COGs.R](../blob/master/R/download_HREA_COGs.R).


## How to Cite

If you use the HREA data or code, please cite our paper:

Brian Min, Zachary P. O'Keeffe, Babatunde Abidoye, Kwawu Mensan Gaba, Trevor Monroe, Benjamin P. Stewart, Kim Baugh, Bruno Sánchez-Andrade Nuño
“Lost in the Dark: A Survey of Energy Poverty from Space,” Joule (2024), [https://doi.org/10.1016/j.joule.2024.05.001](https://doi.org/10.1016/j.joule.2024.05.001)


## Acknowledgements

The HREA project was made possible through the support and partnership of the World Bank, National Oceanic Atmospheric Administration, Microsoft AI for Earth, United Nations Development Programme, National Science Foundation, and the University of Michigan's Advanced Research Computing division.


## Author/maintainer

Zachary O'Keeffe (zokeeffe@umich.edu)


## Contact

For more information, please contact Brian Min (brianmin@umich.edu)


## Links

* [HREA project website](https://hrea.isr.umich.edu/)
* ["Lost in the Dark: A Survey of Energy Poverty from Space"](https://authors.elsevier.com/a/1j8j9925JERZpQ)
* [Light Every Night repository](https://registry.opendata.aws/wb-light-every-night/)
* [Light Every Night tutorial](https://worldbank.github.io/OpenNightLights/wb-light-every-night-readme.html)
* [Achieving Universal Electricity Access](https://data.undp.org/achieving-universal-electricity-access/)
* [UNDP GeoHub](https://geohub.data.undp.org/)
* [Professor Min's personal webpage](https://websites.umich.edu/~brianmin/)


## Notes

We use the term "country" loosely, following the conventions established by ISO 3166-1 alpha-3.  We use the [Database of Global Administrative Areas (GADM)](https://gadm.org/) to define geographic boundaries across all areas of interest.  Both decisions were made to enhance accessibility and usability so that researchers may more easily obtain the data and merge them with other data sources.  We make no claims about rightful territorial ownership or political sovereignty.


## License

HREA code and data are licensed under the GNU GPLv3 license.
