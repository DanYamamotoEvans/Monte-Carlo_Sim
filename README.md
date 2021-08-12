# Monte-Carlo Simulation for biological experiments

This page explains how to install and execute the codes.
If you have any questions, pleaes post a question in the Discussions. 

Please make sure you have appropriate Python and pip before starting.
```sh
Python version >=3.5
pip    version >= 1.1.0
```

Dependencies :
```sh
numpy  version >=1.19 
tqdm   version >=2.2.4
```
To install these pakcages, first clone this repository by
```sh
git clone https://github.com/DanYamamotoEvans/Monte-Carlo_Sim.git
```

Next, go to the location of the Monte-Carlo_Sim folder in the terminal, and install the dependencies by
```sh
pip install .
```

Other core programs to install:
- [Jupyter-notebook](https://jupyter.org/install)
```sh
pip install jupyterlab
```
    
## Overview
This script was built to perform experimental plans. There is a single jupyter-notebooks for each experimental setup. Please add/modify your own experiment. For visualization, I use R (sorry matplot lib people).

- Monte-Carlo simulation (BFG-PCA)
- (Visualization, You will need to install [R](https://cran.r-project.org/).)

### Monte-carlo simulation of BFG screening proccess
Since BFG screenings have multiple sampling steps while handling a complex pool of strains, we suimulate the sampling process with a Monte-Carlo simulation. This notebook follows the procedures of BFG screenings, and allows the user to estimate the nessesary paramaters for sampling. 


### References
- [Yachie _et al_, 2016](https://www.embopress.org/doi/full/10.15252/msb.20156660) / Initial report of BFG. The codes here were built based on perl scripts provided from [Dr. Nozomu Yachie](http://yachie-lab.org/?nozomuyachie).
- [Evans-Yamamamto _et al_, 2021 (Preprint)](https://www.biorxiv.org/content/10.1101/2021.07.27.453987v1) / This repositry was created in part of this work to make BFG-PCA accessible. 
