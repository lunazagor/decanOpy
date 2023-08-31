# decanOpy

decanOpy is an AstroPy-powered code which generates celestial coordinates of a given Egyptian decan in 4-minute intervals for a given year BCE. 

## Install

```
$ git clone https://github.com/lunazagor/decanOpy
```

## Generating Stellar Data

The code takes three inputs: the name of the decan, the year BCE, and the starting month. An example run might look like

```
python3 decanO.py -decan Merak -yearBC 1300 -month 01 
```

The code will produce a .txt file with the following columns:

```
Julian Date|Human Readable Date|<Decan name> Azimuth|<Decan name> Altitude|Sun Azimuth|Sun Altitude
```

The Jupyter notebook containts some examples of how to visualize the data saved in the .txt file. 

## Accompanying StoryMap

The newest version of the StoryMap discussing decanOpy, as presented at the Annual Meeting at ARCE '21:

https://storymaps.arcgis.com/stories/eea3fbc9c05b40948563ffd0ccfab59d
