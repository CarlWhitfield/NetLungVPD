# NetLungVPD
A computational model of ventilation, gas/particle transport, and deposition in lung airway networks.

## To build and execute using Dockerfile

To build an executable using the Dockerfile, install Docker and in the directory Network_lung_model run the command

`docker build -t oo_lung_sim ./`

which will build a local Docker image called `oo_lung_sim`. To execute the code, use the `docker run` command

`docker run -v "/path/to/sim:/sim" -w "/sim" oo_lung_sim ARGS`

where the `-v` option maps the local directory `/path/to/sim` to the virtual directory `/sim` and the `-w` option sets the working directory in the docker image to the virtual directory `/sim`. This means that simulation input files contained in the local folder `/path/to/sim` can be passed to the code, and that simulation outputs will be outputted to this directory. For example, suppose the local directory `/path/to/sim` contains the input files `in.options`, `in.params`, `tree.branches`, `tree.nodes`, and `tree.termnodes`. These can be passed as input to the simulation as follows

`docker run -v "/path/to/sim:/sim" -w "/sim" oo_lung_sim in.options in.params tree.branches tree.nodes tree.termnodes`

## Required command line inputs

Five input files are required to run the executable from the command line: a .options file, a .params file, and three files that define the network geometry. Example options and params files are given in the directory `options_params`. The files defining a network geometry are (i) a .nodes file, which provides the location of each node, (ii) a .termnodes file, which lists the indices of all of the terminal nodes, and (iii) a .branches file, which provides the radius of each edge and the nodes it is connected to. Examples of these are published at https://doi.org/10.5281/zenodo.3709105 .

## Simulation outputs

When a particle deposition or multiple-breath washout simulation is run, the code outputs the following files: a washout.csv file, which contains the value of several key outputs at regular time points through the simulation; deposited_*.csv files, which contain the value of key outputs for each individual airway with each file corresponding to a time point as specified by values in the .params file; transport_tree and flow_tree directories, which contains .vtk files at specified time points which can be used to visualise the output in the 3D network geometry in paraview.
