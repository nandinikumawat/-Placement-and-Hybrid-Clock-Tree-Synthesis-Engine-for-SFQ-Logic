# vda2-project
## Dependencies
- [hMETIS](https://karypis.github.io/glaros/software/metis/overview.html) (Runtime dependency)
- [cgal](https://www.cgal.org/)
## Building

This repository contains two programs: an implementation of FastPlace from EE5301, and 
this [SFQ placement algorithm](https://ieeexplore.ieee.org/document/8342242).

This project uses CMake. To build:
```
mkdir build && cd build
cmake ..
make
```

## Running
### Runtime Dependencies
The SFQPlace executable requires [hMETIS](https://karypis.github.io/glaros/software/metis/overview.html) to run.
Simply make sure that the `khmetis` executable is present in the current working directory when `sfqplace` is ran.

`PA3` is also required to be in the current working directory. 

### Usage
Compiling will generate two executables: `PA3` and `sfqplace`. 
`sfqplace` accepts netlists of the [ISC format](https://davidkebo.com/wp-content/uploads/2023/10/iscas85.pdf).
To run `sfqplace` simply provide the netlist filename **without the extension .isc**.

Example: `./sfqplace c17` to run sfqplace on `c17.isc`

The program will execute and create numerous files. Most importantly are `[netlist]_spread.kiaPad` which is the initial placement of cells prior to the SFQ placement algorithm.
The final placement of supercells can be found in `supercells_spread.kiaPad`.

## Visualization

Included is a python visualization tool to plot the cell locations of a `.kiaPad` file.
Requirements of the script are listed in `requirements.txt`. 

We recommend using [virtualenv](https://virtualenv.pypa.io/en/latest/user_guide.html) to install them and run the script.

Usage:
`python placement_visualizer.py [path to .kiaPad file]`
