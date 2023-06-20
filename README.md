# decanOpy

decanOpy is an AstroPy-powered code which generates celestial coordinates of a given Egyptian decan in 4-minute intervals for a given year BCE. 

## Install

```
$ git clone https://github.com/lunazagor/decanOpy
```

## Usage

The code takes two mandatory inputs: the name(s) of the decan(s) and the year BCE. An example run might look like
```
python3 run.py --decan Sirius Saiph Rigel Hamal Rasalhague Sheratan --yearBC 1300 
```
The complete set of flags the code can take are 
```
python3 run.py --decan [names] --yearBC [year] --month ["01"] --matchStellariumJD [True] --name ["data"]
```
where the last three are the first month of the calculation (January by default), whether to match JD dates to Stellarium (these are offset from Astropy by 10-ish days between 1600 and 1100 BCE; True by default), and the name of the resulting file. 

The code will produce a .txt file with the following columns:
```
Julian Date|Local Date and Time|Sun Azimuth|Sun Altitude|[decan#0] Azimuth|[decan#0] Altitude|[decan#1] Azimuth|[decan#1] Altitude ...
```
with the last two columns repeating for the number of decans listed. The Jupyter Notebook contains some examples of visualizing the data saved in the .txt file. 

## Accompanying StoryMap

The newest version of the StoryMap discussing decanOpy, as presented at the Annual Meeting at ARCE '21:

https://storymaps.arcgis.com/stories/eea3fbc9c05b40948563ffd0ccfab59d

The code has since been updated, as will be documented further. 
