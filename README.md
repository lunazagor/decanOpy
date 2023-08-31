# decanOpy

decanOpy is an AstroPy-powered code which generates celestial coordinates of a given Egyptian decan in 4-minute intervals for a given year BCE. 

## Install

```
$ git clone https://github.com/lunazagor/decanOpy
```

## Generating Stellar Data

The code takes two mandatory inputs: the name of the decan (or list of decan names) and the year BCE. An example run might look like

```
python3 run.py --decan Sirius Saiph Rigel Hamal Rasalhague Sheratan --yearBC 1400 
```

The complete set of flags the code can take are 

```
python3 run.py --decan [names] --yearBC [year] --month ["01"] --matchStellariumJD [True] --name ["data"]
```

where the last three are the first month of the calculation (January by default), whether to match JD dates to Stellarium (these are offset from Astropy by 10-ish days between 1600 and 1100 BCE; True by default), and the name of the resulting file. 

The code will produce a .txt file in `/DecanLists` with the following columns:

```
Julian Date|Human Readable Date|<Decan name> Azimuth|<Decan name> Altitude|....|<Decan name> Altitude|Sun Azimuth|Sun Altitude
```

Alternatively, the user can crate Sirius-like mock data (currently same dec as Sirius, but equally-distributed RAs) by running

```
python3 mockrun.py --yearBC 1300 --num 48
```  

where num is the number of stars to equally distribute. WARNING: the dec and RA of Sirius are currently hard-coded to be at BCE 1300, so the yearBC parameter will actually not change anything. The unlisted optional flags may also require debugging.  

## Visualizing Data

The Jupyter notebook `decanPlotting.ipynb` contains some examples of how to visualize the data saved in the .txt file. Some of its contents may be deprecated. 

The notebook `mock_synRSCs.ipynb` allows the user to generate synthetic Ramesside star clocks, saved in the folder `/SynRSC`. These are automatically saved as `.xlsx` files, so Microsoft Excel is currently necessary to read them. 
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

The code has since been updated and is actively being updated still, as will be documented further. 

## To Do

- check stability of storymap / link to PDF if not

- update `decanPlotting.ipynb` so that it's compatible with both real and mock multi-decan data

- allow option for synRSCs in .csv 

- add more options for mock data generation, remove Sirius hardcoding, make --yearBC a useful a parameter

De