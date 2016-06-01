# NEON-mosquito-data
An analysis of mosquito data from the NEON data portal

## Contents
By running the scripts in this repository, you can explore a portion of the NEON mosquito data [available on the NEON data portal](http://data.neonscience.org). Complete data is always available directly from the NEON website, but this tutorial can be applied to any amount of mosquito data downloaded from the portal. 

### Data
A folder that contains 2 zip files with mosquito data collected at 4 locations in 2014.

### Code
A folder that contains all the scripts used to process/analyze the data and generate the figures in this tutorial. By downloading current data from the NEON data portal and placing it in the 'data' folder, the same analyses can be conducted by using the scripts located in the 'code' folder.
Run each script in the following order.
1. 'getting-the-data.R' 
2. 'analysis.R'
The first script will extract all zip files and compile the data by table. The second script will run the analysis. Before either can be executed, the user must change their filepath (called 'myPathToData') so that it matches where this repository has been cloned.
