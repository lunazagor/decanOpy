# decanOpy

decanOpy is an AstroPy-powered code that generates celestial coordinates of a given Egyptian decan in 4-minute intervals for a given year BCE. 

## Install

Clone the files from this repo from the command line as

```
$ git clone https://github.com/lunazagor/decanOpy
```

or just download repo through GitHub.com. The list of required Python modules is contained in requirements.txt. 

## Generating Stellar Data

The module `decanopy.flow` contains all the machinery to generate a sky (real_sky, rand_sky, or star_like) and recreate the movement around the sky over a calendar year. This part of the code is currently being refactored and is unstable, but we have two output files to work on in the meantime. 

## Creating synthetic Ramesside Star Clocks 

Inside the module `decanopy.models.RSC` are the functions needed to generate synthetic RSCs for comparison with N&P data. An example notebook for calling the functions is stored in `/notebooks/synRSCgenerator.ipynb`, and may be updated as the relevant functions change. 

I highly recommend making a working copy of the notebook in `/notebooks` or `/dev` (which ships with the repo but its contents are kept private). This way, it won't get rewritten or cause version issues whenever a `git pull` is called!  

## Visualizing Data

TBD in `decanopy.visualization`. 

## To Do

Create a better README and proper docs, for one!


