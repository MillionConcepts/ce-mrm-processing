# Introduction

This bundle contains Planetary Data System 4 (PDS4) versions of calibrated
(L2C) data products from the Chinese Lunar Exploration Program
(CLEP)'s Chang'e-1 (CE-1) and Chang'e-2 (CE-2) Microwave Radiometer
(MRM) instruments, along with new products derived from those data. We
produced these products with NASA support under grant #80NSSC20K1430. Note
that NASA is generally disallowed from funding collaborations with Chinese
entities[^1], so we conducted this effort without the participation of the
China National Space Administration (CNSA). We thank NASA for this financial
support. We also thank CNSA, the CE-1 and CE-2 archival teams, China’s Lunar
and Planetary Data Release System (CLPDS), and the National Astronomical
Observatory of China (NAOC) for making these data publicly available to the
global scientific community.

We chose to archive these products in the PDS for two major reasons:

* Almost all scientists who use these data in time-series form want to analyze
  data taken across an entire mission. The NAOC-held L2C products, although
  written in a clean and consistent format, require laborious preprocessing
  before they can be used in these sorts of analyses. We have created tables
  that aggregate the MRM L2C data and add precomputed values for analytically
  crucial quantities like local time. We believe these products will improve
  time-to-value for scientists who wish to incorporate the MRM data into new
  investigations.
* Many scientists prefer to work with data of this type in map-projected
  rather than time-series form, and the map-projected products we include in
  this bundle are significantly higher-quality, better-resolution, and more
  complete than the limited set of map-projected products made available by
  the NAOC.

# Context

## Missions and Hosts

The CE-1 and CE-2 missions comprised the first phase of CLEP’s Chang’e
Project. They were intended to test navigation and communication systems,
produce high-resolution maps of the lunar surface, conduct assorted
geophysics and heliophysics investigations, and more generally enable surface
operations planned for later phases of the Project. They each hosted a
similar suite of surface-observing instruments: a visible-band stereo camera
(CCD), a laser altimeter (LAM), X-ray and gamma-ray spectrometers (XRS and
GRS), and a 4-channel microwave radiometer (MRM). They also each hosted a
heliophysics package: a high-energy particle detector (HPD) and a solar wind
ion detector (SWID). The CE-2 CCD’s design differed significantly from the
CE-1 CCD’s, but all other instruments were identical in construction, subject
to few if any modifications[^2]. In addition to these instruments, CE-1 also
hosted a 32-band NUV-NIR imaging interferometer (IIM). (Huixian et al., 2005;
Zuo et al., 2014).  

CE-1 entered lunar orbit on 2007-11-05. It ceased operation on 2009-03-01 upon
planned impact with the lunar surface. (China National Space Administration,
2009; Wang et al., 2009).  CE-2 entered lunar orbit on 2010-10-06. On
determining that CE-2 would have much more fuel remaining by the end of its
primary mission than initially expected, CNSA designed an extended mission
that re-tasked CE-2 for deep space communications testing and small body
observations. CE-2 exited lunar orbit on 2011-06-08 and subsequently
performed several deep space operations. (NSSDCA, 2022; Zuo et al., 2014).
These operations included a flyby of the asteroid Toutatis; however, Toutatis
was only briefly within the effective range of the MRM, not long enough to
collect clear Toutatis data  (Li, 2013). CE-2 left communications range in
2014, and CNSA estimates that it may reenter communications range as early as
2027 (Jones, 2021).  

## Geophysical Applications

The MRM channels fall primarily into the super high frequency
(SHF) or "centimeter" band (the 37 GHZ channel technically edges into the
extremely high frequency (EHF) band).  A handful of other planetary science
investigations have taken SHF data, including Juno MWR, Rosetta MIRO, and
Cassini RADAR. However, only the CE-1 and CE-2 MRMs have systematically
imaged a large body's surface in this band (not counting coarse observations
by ground-based radiotelescopes). 

SHF technologies are robust, cheap, and extremely mature; terrestrial
applications range from aviation radar to short-range high-bandwidth
communications (like WiFi) to ground-based weather radar to kitchen
appliances. It is also extensively used in radioastronomy. Because water
absorbs SHF very effectively, SHF radiometry has limited use in terrestrial
geology. However, in the absence of water, passive SHF radiometers can
essentially “see” into the top few meters of a body’s surface. Emissions in
this band are principally due to blackbody radiation from variations in
physical temperature, and because the depth of this “vision” is
frequency-dependent (lower frequencies see deeper into the surface), a
multiband radiometer can effectively measure subsurface temperature
gradients. Effective measurement depth is affected by dielectric material
properties as well as frequency, so regional variations in heating and
cooling over a diurnal cycle can also give insight into compositional
properties, particularly metallicity. See Siegler et al., 2020 and Siegler et
al., 2023 for further discussion and examples of application. 

## Antenna Characteristics

The MRM has 4 channels, centered respectively at 3, 7.8, 19.35, and 37 GHz,
with respective bandwidths of 100, 200, 500, and 500 MHz. Nominal radiometric
resolution is 0.5 K for all channels. Angular full width at half maximum
(FWHM) is ~13 degrees for channel 1 and ~10 degrees for the higher channels,
although the main beams are platykurtic and slightly asymmetrical. The first
side lobes are larger and highly asymmetric. Outer side lobes are not
characterized in the literature. Please refer to Wang et al., 2010b for a
more complete description of antenna patterns.

Orbital height varies significantly within each data set. CE-1 mean orbital
height was ~185 km with a standard deviation of ~32 km; CE-2’s was ~100 km
with a standard deviation of ~10 km. This means that surface footprint size
is not constant within either data set. However, at mean CE-2 orbital height
(~100 km), the diameter of the surface footprint of the main beam at 50% of
maximum response (i.e. FWHM) is ~23 km for channel 1 and ~17 km for the
higher channels; at 10%, diameter is ~40 km for channel 1 and ~25 km for the
higher channels. At mean CE-1 orbital height (~185 km), FWHM is ~45 km for
channel 1 and ~30 km for the higher channels; diameter at 10% is ~70 km for
channel 1 and ~45 km for the higher channels. (All of these figures assume
perfect nadir pointing.)

Various values for the MRM’s “spatial resolution” are quoted in the literature
without derivation; they are usually close to FWHM at mean orbital height.
Note, however, that it is not really meaningful to assign a single “spatial
resolution” value to the MRM data. The effective resolution of the data set
is much higher than the “resolution” of a sample considered in isolation, and
is determined as much by sample density as by footprint size. Furthermore,
both mean sample density and mean footprint size vary greatly between
regions.

## Temporospatial Scope

The CE-1 MRM L2C data set contains samples taken between 2007-11-27 and
2009-01-24. Most of these samples exhibit a regular observational cadence: 5
samples separated by 1.6-second intervals followed by 1 at a 3.6-second
interval[^3]. The record has two gaps of ~3 months, several other gaps of 1-5
days, and 92 gaps of ~15 minutes. The data set contains ~6.3 million samples
in total. 

The CE-2 MRM L2C data set contains samples taken between 2010-10-15 and
2011-05-20. They exhibit the same observational cadence as the CE-1 data. The
record has several time gaps ranging in size from ~15 minutes to ~7 days.
There is also a ~2-week interval during a major orbital pattern change in
February 2011 during which the mission flagged all data as off-nominal.
(The calibration for the mission-flagged data appears to us as possibly
compromised, perhaps along with geopositioning.) The data set contains ~8.7
million samples in total, and ~7.5 million excluding the 2-week period of
mission-flagged data.

Both MRM data sets include samples of the entire lunar surface, and each MRM
sampled most regions multiple times in adjacent orbital passes and in
separate sets of passes at other times of day. Note, however, that because
CE-1 held a roughly noon-midnight orbit throughout the observational phase of
the mission, the CE-1 MRM data include few samples taken during the
lunar “morning” or “evening,” particularly in lower latitudes. Conversely,
because CE-2 changed its orbital pattern more dramatically over the course of
its mission, moving from a noon-midnight orbit to a terminator orbit back to
a noon-midnight orbit, the CE-2 MRM data include samples taken across the
entire diurnal cycle (although sample density is greatest near midnight and
noon).

## Calibration

MRM calibration relied on constants determined during ground testing along
with regular warm and cold calibration measurements taken in flight (see Wang
et al., 2010b for details). The CE-2 calibration system was identical in
basic architecture, but ground calibration was significantly more rigorous
and the space-pointing cold horn was tilted up approximately 15 degrees to
address concerns about possible contamination (Feng et al., 2013; also
compare Wang et al., 2010a, 2010b).

Both the CE-1 and CE-2 MRM data exhibit significant offsets in temperature
between orbits (even adjacent orbits), and sometimes between groups of
samples within the same orbit. The data are also quite “noisy”, especially
the CE-1 data; what we interpret as "noise" includes sample-to-sample
variations, a number of obviously “bad” orbits, and many outlier points
within otherwise “good” orbits. Even with best-effort data cleaning, many
uncertainties remain, and we recommend treating relative temperature values
as more reliable in general than absolute values.

Errors in calibration appear to be partly responsible for these offsets,
noise, and inconsistencies. First, both MRM cold horns appear to have
experienced contamination by emission from the lunar surface. Various
approaches to addressing this phenomenon have been proposed (see for instance
Hu et al., 2022, 2017; Tsang et al., 2014), but no complete “corrected”
version of either data set has been published.

Second, the CE-2 MRM exhibits changes in absolute calibration across the
course of the mission (Siegler and Feng, 2017; St. Clair et al., 2022).
Physical modeling suggests that, while both MRMs were subject to surface
contamination effects, this CE-2 specific phenomenon likely had a separate
physical cause (Hu and Keihm, 2021). Feng et al., 2020 postulated solar
contamination of the cold horn (either via direct radiation or indirect
physical heating) as the most plausible explanation. This is supported by a
strong relationship between measured temperatures and the angle between the
CE-2 orbital plane and the Sun-Moon axis (“orbit angle”). This relationship
is clearest (and perhaps only present) for absolute orbit angles below ~45
degrees. This is to say that as CE-2 shifted away from a noon-midnight orbit
towards a terminator orbit, its temperature measurements dropped, and they
rose again as it shifted back towards a noon-midnight orbit. The peak effect
size in channel 1 is ~13K (towards the middle of the large range of
mission-flagged data in February 2011). This is far too large to be explained
by diurnal effects. It also cannot be explained by sampling bias, as it
appears across many different regions of the surface. 

This orbit angle-dependent effect has high statistical error. This is due in
part to the confounding effects of latitudinal cooling, regional temperature
variability, and the diurnal cycle; however, it also appears to have a
complex relationship with orbiter attitude and Sun-relative orientation. For
instance, it is asymmetric with respect to the sign of orbit angle (larger as
the orbiter approached the terminator); i.e., the mean offset is larger
at -15 degrees orbit angle than +15 degrees orbit angle. It is also worth
noting that measurements from the final phase of mapping appear to be
slightly *warmer* on average than measurements taken during early orbits with
similar absolute orbit angles. Due to these complexities, we consider it
unlikely that this calibration offset can be seriously investigated or
addressed without data on orbiter attitude, which CNSA has never released. 

## Georeferencing

The MRM L2C tables include latitude, longitude, and orbital height values for
each sample. Neither documentation or published literature specify a frame of
reference for the lat/lon values, but comparison to other orbital data
products suggests that they are referenced to MOON_ME or IAU_MOON. Because
the millidegree offsets between MOON_ME and IAU_MOON are below even the
best-case spatial resolution of the data, we have arbitrarily chosen to treat
the coordinates as MOON_ME. In the published literature, the lat/lon values
are treated both as (1) boresight intercept point coordinates and(2) directly
nadir to the spacecraft. In other words, analyses of these data assume that
the MRM was always perfectly nadir-facing, such that its boresight intercept
point is also the orbiter’s subsurface point.

However, because, both orbiters performed somewhat complex maneuvers, the MRM
was almost certainly *not* directly nadir-facing for most measurements.
Published literature and documentation do not describe how–or if–these
georeferences were corrected for orbiter attitude, but we have observed
quasi-periodic oscillations in the time-series data that are reminiscent of
typical patterns in gyroscopic attitude determination and control systems
(ADCS). These oscillations typically exhibit much higher amplitude past +/-
75 degrees latitude. This suggests the possibility that the orbiter performed
regular corrective maneuvers as it passed the poles, and that attitude
changes from these maneuvers propagated into measured temperatures, perhaps
via changes in the MRM surface footprint or offsets in calibration or both.
They are also much stronger in some orbits than others. They may be
responsible in part for the higher “noisiness” of the CE-1 data, as CE-1’s
greater average orbital height would have tended to increase the effects of
pointing error. Unfortunately, CNSA has not released orbiter attitude data,
so we are unable to investigate these possibilities in detail. However, we
recommend exercising caution when using data taken past +/- 75 degrees
latitude.

We believe that the data striping and “jitter” visible in some portions of the
map products may be related either to geopositioning/pointing errors or the
calibration issues described above. There are also two effects that we
believe may be due solely to geopositioning/pointing errors. First, there are
small apparent spatial misalignments between some surface features in some
time bins. Second, the effect of topography on observed temperatures appears
to become much stronger at higher latitudes, much more so than can be
explained by changes in effective local time due to topographic dependence of
hour angle (especially in the lower channels).

# Bundle Directory Structure

## /miscellaneous

Root directory of the bundle's miscellaneous collection. This collection 
includes products that were used to help derive some products in the data 
collection, but that were not themselves derived from CE MRM data. The root 
directory contains two products used to help derive the products in this 
collection's "tbmod" subdirectory.

## /miscellaneous/tbmod

Thermal brightness model maps.

## /browse

Bundle's browse collection.

## /browse/[image type]

Root directories per image type for browse images generated from MRM
observational data, where "image type" is "datminus", "latshift", "stdev",
or "temp".

## /browse/[image type]/[orbiter]

Per-orbiter root directories for browse images of a particular type,
where "orbiter" is "ce1" or "ce2".

## /browse/[image type]/[orbiter]/[channel]

Per-channel directories containing browse images of a particular type
generated using data from a particular orbiter, where "channel"
is "t1", "t2", "t3", or "t4", corresponding respectively to the MRM channels
centered at 3, 7.8, 19.35, and 37 GHz.

## /browse/tbmod

Root directory for browse products generated from "tbmod" data products. As
the "tbmod" products were generated from non-MRM observational data, this
directory does not have per-orbiter subdirectories.

## /browse/tbmod/[channel]

Per-channel directories containing browse images of modeled brightness
temperatures, where "channel" is "t1", "t2", "t3", or "t4", corresponding
respectively to the MRM channels centered at 3, 7.8, 19.35, and 37 GHz.

## /data

Root directory for the bundle's primary data collection. Also contains two
observational data products: ce1_mrm.fits/xml and ce2_mrm.fits/xml. These are
binary tables which contain concatenated versions of, respectively, the CE-1
and CE-2 L2C MRM data, along with derived quantities necessary for standard
corpus-level analyses of the data.

## /data/maps

Map-projected CE MRM data products.

## /data_source

Root directory for the bundle's data_source collection.

## /data_source/[orbiter]

Per-orbiter root directories for the L2C source data products, where "orbiter"
is "ce1" or "ce2".

## /data_source/[orbiter]/[YYYY][MM]

Directories containing PDS4-labeled versions of the NAOC-archived L2C data
files, split by 4-digit year and 0-padded 2-digit month.

## /document

Root directory for the bundle's document collection. This collection has no
subdirectories.

# Map-projected Products

## Overview

The bundle's data collection includes 4 types of map-projected
products: "temp", "latshift", "datminus", and "tbmod". This section describes
their contents and organization. More details on the process used to derive
them are given in the "Deconvolved Map Product Processing" section below.
Note that the "tbmod" maps are in the bundle's miscellaneous collection,
as they do not use CE MRM observational data as an input. The others are
in the bundle's data collection.

### Projection

These products contain maps of the lunar surface in equirectangular projection
at 32 pixels per degree, centered on 0 degrees latitude / 0 degrees
longitude. All maps cover -180 to 180 degrees longitude. "temp"
and "latshift" products cover -75 to 75 degrees latitude; "datminus"
and "tbmod" products cover -70 to 70 degrees latitude. All maps are
coregistered with one another aside from the aforementioned 10-degree
latitude crop. All products also include longitude and latitude arrays
corresponding to their x and y axes. Latitude and longitude values are given
in the lunar mean Earth/polar axis reference frame (SPICE MOON_ME).

### Excluded Data

We excluded some samples from the maps, including samples flagged by the
mission, samples from a handful of orbits with errors in temperature or
geopositioning, samples from a handful of extremely “noisy” orbits, samples
with physically implausible temperatures or implied temperature gradients,
and a few groups of samples with duplicate timestamps. Also, the CE-2 set has
a large contiguous range of mission-flagged orbits that appear to have
persistent calibration errors. We extended this range on both sides to
exclude orbits in which those calibration errors appear to persist
(although at lower intensity). See the “Calibration” section above for
further discussion of this phenomenon. 

The concatenated data tables include quality flags that describe the specific
rationales for sample exclusion. See the “Columns” subsection of
the “Concatenated Data Tables” section below for a full description of these
flags.

### Physical Product Structure

Each map product consists of a FITS file and a PDS4 label file. Each product's
individual maps, as well as its latitude and longitude values, are stored as
binary arrays in separate HDUs of its FITS file. We chose this physical
layout because storing multiple named array objects per file simplifies
logical product structure and physical file access, and comes with almost no
performance tradeoffs due to the affordances of the FITS format.

### Usage Notes

* Users who want to use these maps with dedicated GIS software can convert them 
to GeoTIFF using the `gdal_translate` program included with 
[GDAL](https://gdal.org/en/latest/). The GDAL PDS4 driver will include the map 
projection specified in the .xml label in the output GeoTIFF. 
[Detailed instructions for use of the GDAL PDS4 driver are here.](https://gdal.org/en/latest/drivers/raster/pds4.html) 
  * For example, `gdal_translate PDS4:ce1_t4_temp_32ppd.xml:1:3 ce2_t4_temp_4_6.tif`
    will produce a GeoTIFF from the third raster array in ce2_t4_temp_32ppd.fits.
  * For a list of which PDS4 object names correspond to which GDAL 'subdataset'
    number (the '3' in the '1:3' part of the `gdal_translate` command above)
    in a particular file, run `gdalinfo` on the .xml label file.
* FITS normally uses bottom-to-top array order, but we wrote these maps 
top-to-bottom for compatibility with simple raster display tools like PDS4 Viewer 
and `pdr`. However, this means that users of FITS-specific display tools like 
`fv` or DS9 may see the maps flipped bottom-to-top (i.e. South up). If this 
occurs, simply mirror the maps about the x axis / equator.
  * Note that this has no effect on GDAL, because its PDS4 driver respects the 
  `vertical_display_direction` keyword and will interpret the maps North up as
  intended.

### Data Objects

All maps are stored as 2-D arrays of 16-bit MSB integers. Scaling factors and
offsets are provided to convert them to physical units; standard FITS tools
will automatically convert them to 32-bit floats in physical units on load.
All map array values are in Kelvin, although their physical references vary;
they are described in the "Map Product Types" section below.

Latitude and longitude are stored as 1-D arrays of IEEE 32-bit MSB floats.
Units are degrees. We chose to store them as 1-D arrays because, since the
maps are in equirectangular projection, every grid cell in a particular
column has the same longitude value, and every grid cell in a particular row
has the same latitude value. This means that storing them as 2-D arrays
increases file size but provides no additional information. Users who wish to
produce traditional "backplanes" may simply repeat the latitude array across
the x axis and the longitude array along the y axis.

Each of these array objects also has an associated FITS header object. Because
we included all metadata necessary for interpretation in the PDS4 labels, we
left these headers fairly minimal. They primarily contain values that FITS
software needs to correctly load the arrays and values required for
compliance with the FITS standard. They also contain a handful of parameters
used by our processing software. Because some of these values are not
required for interpretation, we did not include them all in the PDS4 labels.
However, they may be of interest to users who wish to run the software
themselves.

### Time Binning

Most scientific uses of the MRM data require knowledge of local time in order
to account for or analyze brightness temperature variation over the diurnal
cycle. For this reason, we defined an array of 2-hour local time bins and
derived separate maps from data taken within each of these bins.
(Maps in "tbmod" products contain brightness temperatures modeled at the
center of the bin.) The array of local time bins starts at midnight. Time
bins are described in HDU and PDS4 object names by their left and right
edges (see the "HDU / PDS4 Object Names" section below).

Note that while all map products define the same array of bins, products
derived from CE-1 data do not actually include maps for the 0400-0600,
0600-0800, 1600-1800, or 1800-2000 bins. This is because CE-1 maintained a
noon-midnight orbital pattern throughout the entire mapping phase of the
mission and took too few observations in these intervals to permit production
of meaningful maps.

### HDU / PDS4 object names

#### Map HDUs

`[maptype]_[binstart]_[binstop]`

"maptype" may be "TEMP", "STDEV", "LATSHIFT", "DATMINUS", or "TBMOD". These
 are described in more detail below.

"binstart" and "binend" are the left and right edges of the map's time bin.
 They may be given as either 1 or 2 digits.

Examples: 
* An HDU named `TEMP_0_2` contains a TEMP map derived from data taken between
  midnight and 2 AM local time. 
* An HDU named `LATSHIFT_12_14` contains a LATSHIFT map derived from data
  taken between noon and 2 PM local time.

#### Latitude/Longitude HDUs

HDUs containing latitude and longitude arrays are simply named "LATITUDE"
and "LONGITUDE".

#### PRIMARY HDUs

The FITS standard requires all FITS files to begin with an HDU
named "PRIMARY". Because we wanted to give every PDS4 data object a
meaningful name, and wanted these names to correspond directly to the HDU
names (i.e. EXTNAME keyword card values) of the FITS files, we followed the
common practice of placing a "stub" PRIMARY HDU at the beginning of each of
these FITS files. These HDUs contain no data. They serve only to allow
software tools to recognize the files as FITS.

## Map Product Filenames

Filenames of "temp", "latshift", and "datminus" products follow the pattern:

`[orbiter]_[channel]_[product type]_32ppd.[extension]`

* "orbiter" may be "ce1" or "ce2". It indicates which orbiter's MRM data were
   used to produce the product.
* "channel" may be "t1", "t2", "t3", or "t4", corresponding respectively to
   the MRM channels centered at 3, 7.8, 19.35, and 37 GHz. Literally, it is
   the name of the column that the deconvolution program(`run_deconv`) loaded
   from the corresponding binary source data table (
   [orbiter]_mrm.fits) described in the "Source Data Products" section
   below..
* "product type" may be "temp", "latshift", or "datminus".
* "extension" may be "xml" (the product's PDS4 label) or "fits" (the product's
   FITS file).

Because "tbmod" products are derived from non-MRM observational data, their
filenames do not include an orbiter code, and instead follow the pattern:

`[channel]_tbmod_32ppd`

"channel" here has the same basic meaning as in the MRM-derived map products.
 Literally, it indicates the frequency value given to the brightness
 temperature modeling program (`map_brightness`).

*Note: the "32ppd" in these filenames does not imply that this bundle contains
 similar map products at other resolutions. We included this field partly for
 consistency with certain model inputs and partly because we may, in the
 future, choose to produce map products at other resolutions.*

## Relationships Between Derived Map Products 

Each orbiter/channel combination has three associated map products: a “temp”
product, a “latshift” product, and a “datminus” product. Each channel also
has an associated “tbmod” product. 

“temp” products are direct outputs of the deconvolution pipeline. They are
 derived from the time-series values in the concatenated data table of their
 corresponding orbiter. “latshift” products are derived from the “temp”
 product of the same orbiter and channel. “datminus” products are derived
 from the “temp” product of the same orbiter and channel and the “tbmod”
 product of the same channel.

Although the source products for “temp”, “latshift”, and “datminus” FITS files
are explicitly stated in their PDS4 labels, they can also be determined
simply from their filenames:

* `[orbiter]_[channel]_temp.fits` is derived from `[orbiter]_mrm.fits`.
* `[orbiter]_[channel]_latshift.fits` is derived from `[orbiter]_
  [channel]_temp.fits`.
* `[orbiter]_[channel]_datminus.fits` is derived from `[orbiter]_
  [channel]_temp.fits` and `[channel]_tbmod.fits`.

The time suffixes of individual HDUs can then be used to match HDUs between
these products. Examples:

* The `LATSHIFT_0_2 HDU` of `ce1_t2_latshift_32ppd.fits` is derived from the
  `TEMP_0_2` HDU of `ce1_t2_temp_32ppd.fits`.
* The `DATMINUS_12_14` HDU of `ce2_t1_datminus_32ppd.fits` is derived from the
  `TEMP_12_14` HDU of `ce2_t1_temp_32ppd.fits` and the `TBMOD_12_14` HDU of
  `t1_tbmod_32ppd.fits`.

“tbmod” products are derived from non-MRM observational data. Specific source
 products are described in their PDS4 labels and this document’s reference
 section.

## Map Product Types

### temp

Unlike the other map product types, "temp" products contain two distinct map
types. TEMP maps contain lunar brightness temperatures derived by
deconvolving MRM antenna temperatures. STDEV maps contain biased estimates of
weighted standard deviation for the antenna temperature values used to derive
the values in their corresponding TEMP maps (these values may be useful for
both physical interpretation and data quality assessment).

Each "temp" product contains either 19 (CE-1) or 27 (CE-2) HDUs:

* PRIMARY (FITS stub)
* one TEMP HDU for each 2-hour time bin with valid data (8 or 12 HDUs total)
* one STDEV HDU for each 2-hour time-bin with valid data (8 or 12 HDUs total)
* LATITUDE
* LONGITUDE

### datminus

The DATMINUS maps in these products contain model-subtracted derived
brightness temperatures. Each is created by subtracting the corresponding
TBMOD map of its matching “tbmod” product from the corresponding TEMP map of
its matching “temp” product (after cropping that TEMP map to the latitude
range of the TBMOD map). Their values can be interpreted as differences
between measured temperatures at each grid cell and expected temperatures at
each grid cell, where this expectation is set by a physical model that takes
into account channel frequency, local time, latitude, thermal conductivity
gradient, and albedo. Because this model does not account for dielectric
loss, much of the regional and diurnal variability in these maps
(especially at channels 1 and 2) is generally dependent on concentrations of
the titanium-containing mineral ilmenite.

Each "datminus" product contains either 11 (CE-1) or 15 (CE-2) HDUs:

* PRIMARY (stub)
* one DATMINUS HDU for each 2-hour time bin with valid data (8 or 12 HDUs
  total)
* LATITUDE
* LONGITUDE

### latshift

The LATSHIFT maps in these products contain latitudinally adjusted derived
brightness temperatures. Each is created by fitting a simple latitudinal
model of temperature (`a + cos(lat) ** b`, with `a` and `b` as free
variables) to the corresponding TEMP map of its matching “temp” product, then
subtracting the fitted curve from the TEMP map. Their values are differences
between measured temperatures at each grid cell and expected temperatures at
each grid cell, where this expectation is set *only* by the latitudinal trend
in data taken by this particular orbiter in this particular channel during
this particular range of local times and incorporated into this particular
TEMP map. 

We provided these maps in addition to the DATMINUS maps because there are many
confounders in these data, including simple sample selection bias produced by
offsets between different orbits, that make it difficult to model the joint
dependence of temperature on LTST and latitude. By restricting the population
of samples in the model, the LATSHIFT maps provide a stricter,
although *strictly relative*, view of temperature variation. This means that
some features may be apparent in these maps that are not apparent in the
DATMINUS maps.

Each "latshift" product contains either 11 (CE-1) or 15 (CE-2) HDUs:

* PRIMARY
* a LATSHIFT HDU for each 2-hour time bin with valid data (8 or 12 HDUs
  total)
* LATITUDE
* LONGITUDE

### tbmod

The TBMOD maps in these products contain modeled brightness temperatures
derived by interpolating Kaguya SELENE and LRO Diviner, LOLA and LROC data
within the parameter space of a physical heatflow model. This process is
described in more detail in the “Thermal Brightness Modeling” section below.
Because they do not use CE MRM data as an input, there is only one “tbmod”
product for each channel.

Each “tbmod” product contains 15 HDUs:

* PRIMARY
* one TBMOD HDU for each 2-hour time bin (12 HDUs total)
* LATITUDE
* LONGITUDE

# Source Data Products

## Overview

All products in this bundle’s data collection are based on the contents of the
NAOC-held datasets `CE1_MRM_2C_20071127161816_20090114192053_B` and
`CE2_MRM_2C_20101015085002_20110520125150_A`. These datasets present the MRM
L2C data as fixed-width ASCII tables, one table per orbit in which the MRM
took observational data, each table in a separate file. These files have
attached labels written in valid Parameter Value Language (PVL). Considered
as data products, they are generally compliant with version 3.8 of the PDS
Standards. The CE-1 dataset contains 1691 such files; the CE-2 dataset
contains 2402. The Chang’e Project’s processing level definitions are similar
to NASA EOSDIS definitions. Like EOSDIS L2 products — and most other L2
products across the Chang’e Project data corpus — the MRM L2C products
consist of georeferenced time-series samples in physical units. 

We expect that most users interested in working with the MRM data at
time-series level will prefer to use the concatenated tables in this bundle’s
data collection. However, to help ensure the reproducibility of our work
(as well as a large body of existing scientific literature), we have included
PDS4-labeled versions of these source products in this bundle’s data_source
collection. We have omitted only a handful of “empty” products from the CE-1
dataset that contain a PVL label but no table. As of writing, we are not
aware of any other copy of these data in a PDS-equivalent archive.

Because the fixed-width tables are consistently formatted and the attached PVL
labels are clearly separated from the tables, these files require no physical
modification to be described as valid PDS4 Product_Observational objects with
associated Header and Table_Character objects. We have therefore chosen to
interfere with them as little as possible. Although no checksums are
available, all data files in the data_source collection are, to the best of
our knowledge, byte-level equivalent copies of the L2C files held by the
NAOC. We have also retained their original filenames, including the
unconventional “.2C” extension. (We have, however, organized them into
separate volumes by month and year to make navigation easier.)

## Source Data Product Filenames

All source data filenames follow the pattern:

`[orbiter]_BMYK_MRM-L_SCI_P_[start time]_[stop time]_[orbit number]_
[letter].2C`

* “orbiter” may be “CE1” or “CE2”.
* “start time” and “stop time” provide nominal time bounds for the product in
   UTC time scale, expressed as YYYYMMDDHHMMSS. For
   instance, “20071127182559” is equivalent to “2007-11-27T18:25:59Z”. Note
   that these are not fixed to the first and last timestamps in the data
   tables (although they are usually within a few seconds), and should likely
   be understood as time bounds for the mission-designated orbit.
* “orbit number” is the mission-designated orbit number as a 0-padded 4-digit
   number.
* “letter” is always “B” for CE-1 products and “A” for CE-2 products.

Our PDS4 labels for these products share these filename stems (converted to
lowercase in accordance with PDS4 conventions).

## Table Contents and Format

CE-1 and CE-2 table formats are semantically equivalent, although column names
and boundaries differ slightly. (Please refer to their PDS4 or PVL labels for
details.) The columns contain:

1. Time UTC, expressed in ISO 8601 format to millisecond precision 
2. Channel 1 brightness temperature in Kelvin, given to two decimal places 
3. Channel 2 brightness temperature in Kelvin, given to two decimal places 
4. Channel 3 brightness temperature in Kelvin, given to two decimal places 
5. Channel 4 brightness temperature in Kelvin, given to two decimal places 
6. Signed solar incidence angle in degrees, given to four decimal places 
7. Solar azimuth angle in degrees, given to four decimal places 
8. Latitude in degrees, given to four decimal places 
9. Longitude in degrees, given to four decimal places
10. Orbital height in kilometers, given to six decimal places 
11. Quality state, given as an encoded string (possibly a malformatted hexadecimal
integer) in CE-1 tables and a 0-padded 2-digit decimal integer in CE-2
tables. The exact meanings of specific quality state values are not included
in documentation, but any value of this field other than “0X000000”(for CE-1)
or “00” (for CE-2) appears to indicate severely off-nominal data.

There are short missing time spans in some tables and large time gaps between
some tables, but most of the samples are separated by only 1.6 or 3.6
seconds (mean of ~2 seconds).  Most tables have about 4000 samples / rows.

# Concatenated Data Tables

## Rationale

Each NAOC-released L2C file contains data from a single orbit, but most
applications for the time-series MRM data require analysis of data from many
orbits. This means that scientists who wish to work with these data must
write their own preprocessing code to load millions of rows from thousands of
separate text files and transform them into a format suitable for their
specific application. Most applications also require scientists to write
additional custom code to derive various geometric quantities. These coding
tasks present a barrier to use and run the risk of propagating unique
preprocessing errors into downstream code. We have spent a lot of effort on
these chores in support of our own investigations. To remove this burden from
future researchers, we have provided two preprocessed tables, one per
mission, in a tabular format suitable for most anticipated uses of the
time-series data.

## Concatenated Table Filenames

There are exactly two products of this type: ce1_mrm.fits/xml, which contains
data from the CE-1 MRM, and ce2_mrm.fits/xml, which contains data from the
CE-2 MRM.

## Concatenated Table Contents

Each concatenated table product consists of a PDS4-labeled FITS file. The
TABLE HDU of each FITS file contains a binary table of annotated MRM data.
Each row of these tables corresponds to one row from an L2C ASCII table
(described in the “Source Data Products” section above). Equivalently, each
row corresponds to a single MRM sample. The CE-2 table has ~8.7 million rows;
the CE-1 table has ~6.3 million. They exclude only a handful of rows from the
source data tables that contain special constants or out-of-bounds values
(about 50 in total).

Note that these tables discard the solar azimuth and incidence fields from the
source tables due to the presence of some implausible values. They replace
them with SPICE-computed solar position vectors and local true solar time
(LTST) at the boresight intercept point. They also do not directly copy the
quality flag columns, as there is no public documentation on the meaning of
the different flags; they instead propagate the presence of any nonzero
quality flag value as an element of the FLAG column bitmask.

In addition to LTST and solar position vectors, the tables contain several
other quantities that are not directly copied from the source files: surface
normal vectors at boresight intercept point, time expressed as seconds since
J2000 (for compatibility with SPICE), and a quality flag bitmask.

A complete list of columns follows. Note that all integer and float values are
stored in MSB order.

### Columns

1. ORBIT (2-byte unsigned integer, unitless): Mission-specified orbit number,
taken from names of source files. 
2. UTC (23-byte UTF-8 string, unitless) ISO 8601-formatted time string, taken 
from source tables. 
3. ET (8-byte float, seconds):  Ephemeris Time (offset in seconds from J2000,
commonly used in the SPICE ecosystem). Derived from UTC and NAIF ephemerides. 
4. LTST (4-byte float, unitless): Local true solar time, normalized to 0-1, 
derived from ET, LON, and NAIF ephemerides. 
5. T1-T4: (4-byte float, K, 4 columns total): antenna temperature at channels 
1-4, taken from source tables. 
6. LAT (4-byte float, degrees): Selenodetic latitude of boresight position in
MOON_ME, taken from source tables. 
7. LON (4-byte float, degrees): Selenodetic longitude of boresight position in
MOON_ME, taken from source tables and converted from 0-360 to -180-180 
representation. 
8. D (4-byte float, km): Distance from orbiter to boresight intercept point, 
taken from source tables. 
9. BX, BY, BZ (8-byte float, km, 3 columns total): X, Y, and Z components of 
boresight intercept point relative to Sun body center in ecliptic coordinates 
(ECLIPJ2000), derived from LAT, LON, ET, and NAIF ephemerides. 
10. CX, CY, CZ (8-byte float, km, 3 columns total): like BX/BY/BZ, but for 
orbiter position rather than boresight intercept point. 
11. MX, MY, MZ (8-byte float, km, 3 columns total): like BX/BY/BZ, but for Moon
body center rather than boresight intercept point. 
12. NX, NY, NZ (4-byte float, unitless, 3 columns total): X, Y, and Z components 
of unit normal vector at boresight intercept point in body-fixed
(MOON_ME) coordinates. 
13. FLAG (2-byte unsigned integer, unitless): Quality flags expressed as a 
bitmask. These are our own flags and are not directly propagated from the source 
data. For most analytic purposes, we recommend discarding all samples with a 
nonzero value in this field, with the possible exception of 0b1000000. Individual
(summable) values denote:
   * 0b1: Sample flagged as bad by data providers. Note that this value is
    not used in the CE-1 table. This is because there are very few flagged rows
    in the CE-1 L2C source files (~30, as opposed to ~1.5 million in CE-2), and
    all of these rows contain invalid special constant values, so we did not
    include them in the concatenated table.
   * 0b10: Physically implausible temperatures (< 34K) at one or more channels.
   * 0b100: Physically implausible offset (> 75 K) between channel
     temperatures.
   * 0b1000: Sample is from an orbit we have manually flagged due to large
     ranges of implausible temperature or geopositioning data.
   * 0b10000: Sample is from an orbit almost completely flagged by the original
     data providers or with very few provided samples (the surviving samples in
     these orbits are uniformly suspicious).
   * 0b100000: Sample's timestamp duplicates the timestamp of another sample in
     the dataset (the data may or may not be duplicated, but it is unclear what
     the orbiter was actually doing at that point).
   * 0b1000000: Sample is from an orbit range we have flagged on suspicion of
     secular offsets in calibration that make its temperature values poorly
     comparable with other portions of the data set. Not used in CE-1 table. See
     the “Calibration” section above for details on the CE-2 intermission
     calibration offsets that led us to apply this flag.
   * 0b10000000: Sample is from an orbit with implausibly low or high latitudinal
     variation in channel 1 temperature (after excluding data with prior
     flags).
   * 0b100000000: Sample is from an orbit with suspiciously large numbers of
     individual outliers in channel 1 temperature across a variety of
     latitudinal ranges (after excluding data with prior flags).

## PRIMARY HDUs

The concatenated table FITS files each begin with a “stub” PRIMARY HDU for
compliance with the FITS standard. These HDUs contain no data and no metadata
relevant to interpretation. they are present only to make the files valid
FITS.

# Deconvolved Map Product Processing

## Rationale

There is no single, widely-used map-projected version of the MRM data. The
NAOC holds a set of equirectangular maps of the CE-1 data at 1 pixel /
degree. (The dataset title is CE1-MRM-SDP-V1.0. See Zheng et al., 2012 for a
discussion of these products (the dataset itself is cited as GSDSA 2008 in
this document’s bibliography). However, no similar maps are available for
CE-2, and, in practice, most scientific investigations that approach the MRM
data in map-projected form begin by creating ad hoc “basemaps” directly from
the L2C data. These are generally equirectangular projections produced by
creating a lat/lon grid, selecting a temperature channel, and defining the
brightness temperature of each grid cell as the mean of the temperature
values at that channel for all samples whose nominal boresight intercept
coordinates fall within that grid cell (potentially after selecting a subset
of and/or applying an inline correction to the samples). 

## On Binning and Averaging

Although this "binning and averaging" (BAA) technique is very common across
remote sensing disciplines, it can introduce unpredictable distortions, and
there are several reasons that BAA is a questionable mapmaking technique for
the CE MRM data in particular. (These observations are also applicable, at
least in part, to any other spatially-referenced time-series data in which
individual samples have spatial extent on the scale of the spatial frequency
of the data, the spatial gradient of the measured quantity, or the size of
the grid cells.)

BAA makes the following implicit assumptions:

1. A sample provides information about the grid cell in which its georeference
falls, and no other cells. Equivalently, no sample whose georeference falls
outside a grid cell provides information about that cell. 2. Every sample
whose georeference falls within a particular grid cell provides the same
amount of information about that grid cell.

These are not good assumptions for the CE MRM unless you are
making *extremely* coarse maps. This is because the MRM antenna patterns are
large enough that radiative emissions from a large swath of the lunar surface
contribute to each antenna temperature measurement. This means that every
sample provides information about a large swath of the lunar surface. Because
lunar brightness temperatures can vary significantly on this spatial scale,
it also means that looking at a single antenna temperature measurement cannot
tell you whether the MRM “saw” a wide or narrow range of temperatures within
that swath.

Applying BAA to these data discards much of the evidence they provide and
overweights the rest. Samples whose footprints cover a grid cell but whose
georeferences do not fall within that cell do not contribute to the derived
value for that cell. This sharply reduces the *apparent* resolution and
coverage of the data: a particular map may have many “holes” in extremely
well-measured areas simply because no sample’s georeference fell within those
cells. This means that insignificant differences in orbital swath
patterns–down to tens of meters in some cases–can make dozens of kilometers
of the surface appear unmeasured. It also reduces the *effective* resolution
of the data, because the presence of a feature within a particular grid cell
may be apparent only from its aggregate contributions to multiple samples
whose georeferences fall outside the cell (which BAA specifically cannot
relate to that particular cell).

The fact that BAA discards the contributions of samples to the derived values
of adjacent cells also means that it overweights the contributions of samples
to their assigned cells. This tends to amplify noise, because it includes
outlier sample values in only one cell instead of distributing them across
multiple cells. It can also arbitrarily shift features into adjacent cells–or
even fragment them into multiple adjacent cells in different directions–when
they were measured principally “from” those cells.

These effects are exacerbated by the fact that they vary unpredictably, in
both strength and location, depending on the chosen grid. They tend to be
worse at higher grid resolutions, but because they depend on essentially
arbitrary assignment of samples to cells, even a very small change in
resolution can create artifacts in different places, and some lower
resolutions might have worse distortion than some higher resolutions. In
equirectangular projection, these effects become worse at higher latitudes,
because each grid cell represents a smaller total area, meaning that BAA
discards more information from each sample. And finally, they have complex
mission-phase-dependent variability due to the fact that the size of the MRM
footprint varies with the orbiter’s distance from the surface.
 
Although we have used BAA previously, we believe it is an inferior approach.
It is something like making a photomosaic by taking a series of orbital
images, arbitrarily cropping some of them, and then blurring remaining
overlapping regions together rather than attempting to coregister them. The
CE MRM data have more “real” resolution than is apparent under BAA. Data
anomalies make it uncertain, but based on coregistration with known lunar
features, we believe meaningful resolution down to 32 pixels/degree is
available in densely-sampled regions of the surface, especially in CE-2
channels 2, 3, and 4. (By "meaningful resolution", we refer to feature
centroids and gradients; many features are still "blurred" relative to data
from instruments with higher per-sample resolution.) Also, although this
resolution is not consistently present across the MRM data sets, the method
we present here is, at worst, a strictly-correct upsampling technique
(if performed to adequate orbiter-referenced angular resolution). It
introduces fewer distortions into comparisons between the MRM data and
high-resolution data from other instruments than generic resampling methods.

## Deconvolution Method

Consider each antenna temperature value as a surface integral over the product
of two functions: one function represents the brightness temperature of the
Moon, and the other represents the MRM antenna pattern with its center
shifted to the boresight intercept point. In other words, each antenna
temperature value is a convolution of lunar brightness temperature and
antenna response. We know the antenna response curve, so we can reconstruct
the brightness temperature of the Moon by deconvolving georeferenced antenna
response from antenna temperature. Because the brightness temperature
function cannot be written in closed form, there is no way to solve this with
clever symbolic manipulation. However, it *is* possible to back out
brightness temperature by numerically approximating the integrals. To do so,
we perform the following steps (ignoring many implementation details):

1. Define an equirectangular latitude / longitude grid on the lunar surface.
2. Define three arrays referenced to this grid: W (weight), WT
   (weighted temperature) and WS (weighted squared temperature). Set all
   elements of these arrays to 0. 
3. Define the antenna pickup pattern in instrument-relative elevation / 
   azimuth coordinates, with high resolution relative to orbital height and 
   latitude / longitude resolution. 
4. Then, at each sample: 
   1. using boresight intercept coordinates and orbital height as input 
   variables, apply a reverse orthographic projection to the antenna
   pickup pattern grid to transform it into latitude/longitude coordinates.
   2. Compute a binned sum of the projected antenna pattern elements on the
   latitude / longitude grid. Call that array A. Call the temperature value of
   the sample V. Then:
      * W = W + A
      * WT = WT + A * V
      * WS = WS + A * V^2 
5. When you are done with all the samples, define the output quantities like this:
   * WT / W is the derived brightness temperature at each grid cell.(These are
     the TEMP HDUs in our “temp” deconvolved map products.)
   * sqrt(WS * W - WT^2) is the standard deviation of measured temperature at
     each grid cell. (These are the STDEV HDUs in our “temp” deconvolved map
     products.)

### Notes

* We chose an equirectangular map projection because we wanted to make
  comparisons to existing equirectangular maps based on data from other
  instruments, but this basic method is independent of the choice of map
  projection. To modify the algorithm to work with other map projections, you
  could either change the coordinate system transformation function in step
  2, add an additional step between 2 and 3 to convert lat/lon into the
  coordinates of your map projection, or simply deal with the inconveniences
  of counting on an irregular grid in step 3.
* We considered only the antenna main beams, not their sidelobes. This is
  because the MRM sidelobes are highly asymmetrical and no spacecraft
  attitude data are available, so we believed attempting to model the
  sidelobes would simply cause spatial ‘smearing’ of unknown and variable
  direction and magnitude.  
* The method we describe here approximates a deconvolution as a massive
  weighted average. As such, it is perfectly reasonable to think of BAA as a
  special case of this method that treats each measured value as an integral
  over the product of the measured quantity and a Dirac delta function
  shifted to the center of the nearest grid cell. This allows you to skip the
  complicated parts by setting almost everything to 0 or 1.
* This is formally similar to "super-res" techniques in visible-band imaging.
  In this case, however, we do not consider it an "image enhancement"
  technique, but rather a literally correct--even brute-force--approach to
  the problem. Each measurement always already provides information about an
  extended area.
* We are not aware of prior work with the CE MRM data that takes this
  approach. We suspect that other investigators may have considered it but
  found it computationally intractable. Although the procedure is
  conceptually simple, the MRM data sets include millions of data points, and
  approximating the surface integrals to useful precision without introducing
  discontinuities requires performing expensive operations on multiple arrays
  with millions of elements at each of those millions of points. We discarded
  our first, straightforward implementation because processing all of the
  data with it would have broken our timeline or our budget or both. To
  create a practical implementation of this procedure, we had to perform
  intensive, multifaceted software development efforts, ranging from finicky
  optimization of application-specific calculations to developing a
  general-purpose dynamic numerical approximation library (St. Clair, 2024).
  For full implementation details, please refer to the deconvolution source
  code (St. Clair et al., 2024).

### Related Work

There are at least two distinct existing approaches to dealing with the wide
footprint of the MRM data. One is a deconvolution approach of a very
different character that applies a maximum entropy-based image reconstruction 
technique to some unspecified “original map” of the CE-2 data 
(Xing et al., 2015). We suspect our approach is more valid; 
however, as the authors did not fully specify their inputs, and made neither 
their code nor their data products available, we have not been able to 
meaningfully assess their work.

Another approach, which is applicable only to analyses that fuse the MRM data
with other observational or model data, borrows a common technique from
imaging camera data processing and convolves the comparison data with MRM
antenna patterns. However, time-series radiometer data is very little like
imaging camera data. Because the “image” is built up from many point samples,
convolving comparison data with the radiometer antenna pattern makes
comparison of feature extents more realistic, but leaves in all the noise and
loss effects associated with BAA. We have used this technique in some
previous projects, but now consider it inferior to “correct” deconvolution.


# Thermal Brightness Modeling

This section provides a basic description of how we produced the modeled
thermal brightness temperatures used in our “tbmod” and “datminus” maps. It
is a variation of the technique described in Siegler et al., 2023.

We first create a forward model of physical temperature using an iterative,
self-equilibrated heat flow model (the physical assumptions are described in
Hayne and Aharonson, 2015; Hayne et al., 2017 describes the immediate
precursor to our software implementation of the model). This model takes
latitude, albedo, and h-parameter (a metric of density / thermal conductivity
change with respect to depth) as inputs and returns physical temperature
by “depth” (really a tuple of strictly covariant depth-related quantities)
and local time.

To do so, we define sets of latitude, albedo, and h-parameter values that span
the range of our external inputs and desired map bounds, then execute one
instance of this model for each element of their Cartesian product
(~7500 total for the products in this bundle). The outputs of this set of
model executions form a discrete, bounded 5-dimensional parameter space in
which each element has an associated scalar representing expected physical
temperature.

We then map a simple brightness temperature model across the latitude, albedo,
local time, and h-parameter axes of this space. This model takes temperature
gradient with depth, material :dielectric parameters, and emission frequency
as inputs and returns brightness temperature. We execute this model four
times (once per MRM channel) for each element of these axes—the model
essentially “collapses” all elements along the depth axis at each of these
points into a scalar brightness temperature. Note that analytic difficulties
precluded us from creating dielectric property estimates that were compatible
with the deconvolution process, so we treated dielectric properties as
constant, using mean values for the Moon as a whole. (Note that this differs
from the analyses in Siegler et al. 2020 and Siegler et al. 2023, which
included fit values of the localized loss tangent to remove the effect of
varying composition. It is also why titanium content pops out so clearly in
some of our products.)

Finally, we ingest derived maps of h-parameter, albedo, slope, and azimuth
based on observational data from Kaguya SELENE and LRO Diviner, LOLA, and
LROC then coregister them to the intended scale and lat/lon bounds of our
output maps. For each element of our array of local time bins, we adjust
effective local time using the slope/azimuth data (this essentially treats
hour angle as “real” local time); then, for each of the MRM channel center
frequencies, we produce derived temperature values by interpolating latitude,
h-parameter, albedo, local time, and frequency within the parameter space
produced by the brightness temperature model. These are the TBMOD maps in
the “tbmod” products. Our “datminus” maps are produced by subtracting the
values in these maps from the TEMP maps in the “temp” products.

(The slope and azimuth maps are available from the PDS LOLA node and
are referenced in the "tbmod" product labels. The albedo and h-parameter maps 
are, respectively, wac_hapke_604nm_70N70S_64ppd.fits and hpar_global_128ppd.fits, 
and are both included in this bundle's miscellaneous collection.)

# Bibliography

China National Space Administration, 2009. China’s lunar probe Chang’e-1 
impacts moon [WWW Document].  

Feng et al., 2013. Data processing and result analysis of CE-2 MRM. Earth 
SciJ China Univ Geosci 38, 898–906.  

Feng, J., Siegler, M.A., Hayne, P.O., 2020. New Constraints on Thermal and 
Dielectric Properties of Lunar Regolith from LRO Diviner and CE‐2 Microwave 
Radiometer. J. Geophys. Res. Planets 125, e2019JE006130. 
https://doi.org/10.1029/2019JE006130  

GSDSA, 2008. CE1-MRM-SDP-V1.0. https://doi.org/10.12350/CLPDS.GRAS.CE1.MRM-03.vA  

Hayne, P.O., Aharonson, O., 2015. Thermal stability of ice on Ceres with 
rough topography. J. Geophys. Res. Planets 120, 1567–1584. 
https://doi.org/10.1002/2015JE004887  

Hayne, P.O., Bandfield, J.L., Siegler, M.A., Vasavada, A.R., Ghent, R.R., 
Williams, J., Greenhagen, B.T., Aharonson, O., Elder, C.M., Lucey, P.G., 
Paige, D.A., 2017. Global Regolith Thermophysical Properties of the Moon 
From the Diviner Lunar Radiometer Experiment. J. Geophys. Res. Planets 122, 
2371–2400. https://doi.org/10.1002/2017JE005387  

Hu, G., Keihm, S.J., Wang, Z., 2022. An In-Flight Recalibration for Chang’E-1
and E-2 Microwave Radiometer Datasets Based on Highland Thermophysical Models. 
IEEE Trans. Geosci. Remote Sens. 60, 1–21. https://doi.org/10.1109/TGRS.2021.3125714  

Hu, G.-P., Chan, K.L., Zheng, Y.-C., Tsang, K.T., Xu, A.-A., 2017. Comparison
and evaluation of the Chang’E microwave radiometer data based on theoretical
computation of brightness temperatures at the Apollo 15 and 17 sites. Icarus 
294, 72–80. https://doi.org/10.1016/j.icarus.2017.04.009  

Hu, G.P., Keihm, S.J., 2021. Effect of the Lunar Radiation on the Cold Sky
Horn Antennas of the Chang’E-1 and -2 Microwave Radiometers. IEEE Geosci.
Remote Sens. Lett. 18, 1781–1785. https://doi.org/10.1109/LGRS.2020.3007606   

Huixian, S., Shuwu, D., Jianfeng, Y., Ji, W., Jingshan, J., 2005. Scientific
objectives and payloads of Chang’E-1 lunar satellite. J. Earth Syst. Sci. 114,
789–794. https://doi.org/10.1007/BF02715964  

Jones, A., 2021. China to launch a pair of spacecraft towards the edge of the
solar system. SpaceNews.  

Li, C. & Li, H, 2013. Chang’e 2 Flyby of Toutatis (presentation). SBAG 8 (2013).
https://www.lpi.usra.edu/sbag/meetings/jan2013/presentations/sbag8_presentations/TUES_0930_CE_Toutatis.pdf

NSSDCA, 2022. Chang’e 2 [WWW Document]. NASA Space Sci. Data Coord. Arch.

Siegler, M.A., Feng, J., 2017. Microwave Remote Sensing of Lunar Subsurface
Temperatures: Reconciling Chang’e MRM and LRO Diviner 1705.  

Siegler, M.A., Feng, J., Lehman-Franco, K., Andrews-Hanna, J.C., Economos,
R.C., Clair, M.St., Million, C., Head, J.W., Glotch, T.D., White, M.N., 2023.
Remote detection of a lunar granitic batholith at Compton–Belkovich. Nature 620,
116–121. https://doi.org/10.1038/s41586-023-06183-5  

Siegler, M.A., Feng, J., Lucey, P.G., Ghent, R.R., Hayne, P.O., White, 
M.N., 2020. Lunar Titanium and Frequency‐Dependent Microwave Loss Tangent as 
Constrained by the Chang’E‐2 MRM and LRO Diviner Lunar Radiometers. J. 
Geophys. Res. Planets 125, e2020JE006405. https://doi.org/10.1029/2020JE006405  

St. Clair, M., 2024. quickseries. https://doi.org/10.5281/ZENODO.12702826

St. Clair, M., Brown, S., Million, C., Feng, J., Siegler, M., 2024. 
ce-mrm-processing. https://doi.org/10.5281/ZENODO.12709958    

St. Clair, M., Million, C.C., Ianno, A., Feng, J., Siegler, M., 2022. 
Approaches to Production of Intermediate Data Products for Characterizing 
Systemic Anomalies in the Chang’e-2 Microwave Radiometer Data, in: 53rd Lunar 
and Planetary Science Conference. Presented at the Lunar and Planetary Science 
Conference.  

Tsang, K.T., Fu, S., Gu, J., Zhou, M., Chan, K.L., Zheng, Y.C., 2014. A 
statistical learning approach to Chang’E microwave radiometer data 
calibration, in: 2014 11th International Conference on Fuzzy Systems and 
Knowledge Discovery (FSKD). Presented at the 2014 11th International 
Conference on Fuzzy Systems and Knowledge Discovery (FSKD), pp. 374–378. 
https://doi.org/10.1109/FSKD.2014.6980863  

Wang, M., Shi, X., Jian, N., Yan, R., Ping, J., 2009. Real time monitoring 
of the Chang’E-1 lunar orbit insertion, in: 2009 15th Asia-Pacific Conference
on Communications. Presented at the 2009 15th Asia-Pacific Conference on 
Communications, pp. 442–445. https://doi.org/10.1109/APCC.2009.5375596  

Wang, Z., Li, Y., Zhang, D., Jingshan, J., 2010a. Prelaunch Calibration
of Chang’E-2 Lunar Microwave Radiometer, in: 2010 International Conference on 
Microwave and Millimeter Wave Technology. pp. 1551–1554.

Wang, Z., Li, Y., Zhang, X., JingShan, J., Xu, C., Zhang, D., Zhang, W., 
2010b. Calibration and brightness temperature algorithm of CE-1 
Lunar Microwave Sounder (CELMS). Sci. China Earth Sci. 53, 1392–1406. 
https://doi.org/10.1007/s11430-010-4008-x  

Xing et al., 2015. The deconvolution of lunar brightness temperature based 
on the maximum entropy method using Chang’e-2 microwave data. Res. Astron. 
Astrophys. 15., pp. 293-304.

Zheng, Y.C., Tsang, K.T., Chan, K.L., Zou, Y.L., Zhang, F., Ouyang, Z.Y., 2012. 
First microwave map of the Moon with Chang’E-1 data: The role of local time 
in global imaging. Icarus 219, 194–210. 
https://doi.org/10.1016/j.icarus.2012.02.017  

Zuo, W., Li, C., Zhang, Z., 2014. Scientific data and their release of 
Chang’E-1 and Chang’E-2. Chin. J. Geochem. 33, 24–44. 
https://doi.org/10.1007/s11631-014-0657-3  

# Footnotes

[^1]: "Pursuant to The Department of Defense and Full-Year Appropriation Act,
Public Law 112-10, Section 1340(a); The Consolidated and Further Continuing
Appropriation Act of 2012, Public Law 112-55, Section 539; and future-year
appropriations (hereinafter, "the Acts"), NASA is restricted from using funds
appropriated in the Acts to enter into or fund any grant or cooperative
agreement of any kind to participate, collaborate, or coordinate bilaterally
with China or any Chinese-owned company, at the prime recipient level or at
any subrecipient level, whether the bilateral involvement is funded or
performed under a no-exchange of funds arrangement."

[^2]: Most CE-2 instruments were in fact flight spares from CE-1, including
the MRM. Therefore, unless otherwise specified, all technical specifications
given in this document for the MRM apply equally to both CE-1 and CE-2.

[^3]: The MRM’s observational cadence is clearly visible in some of the
deconvolved maps, particularly in the CE-2 data at lower latitudes. Because
the spatial offset “jumps” every sixth sample, the measured temperature
gradient also “jumps” every sixth sample, creating evenly-spaced “blocks”
along each orbital track.