# decanOpy

decanOpy is an AstroPy-powered code which generates celestial coordinates of a given Egyptian decan in 4-minute intervals for a given year BCE. 

## Install

```
$ git clone https://github.com/lunazagor/decanOpy
```

## Usage

The code takes three inputs: the name of the decan, the year BCE, and the starting month. An example run might look like

```
python3 decanO.py -decan Merak -yearBC 1300 -month 01 
```

The code will produce a .txt file with the following columns:

```
Julian Date|Human Readable Date|Merak Azimuth|Merak Altitude|Sun Azimuth|Sun Altitude
```

The Jupyter notebook (to be uplaoded soon!) containts some examples of how to visualize the data saved in the .txt file. 
