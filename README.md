The mlpack package should be installed firstly.
How to install mlpack can be found in address: https://github.com/mlpack/mlpack/.

# compile

     make

When 'make' command is fininshed, 'StructureEncode' file is generated. 

# Comands 


     ./StructureEncode Encode  DistanceAdjacencyMatrix.csv   EnergyAdjacencyMatrix.csv    OutputFile   [-config  configFile ]


# config file

     Tolerance = 0.01
     VectorSize = 128
     
   
   The Tolerance is used to determine the size of represent dictionary for distance adjacency matrix. If Tolerance is specified, threshold = 1-Tolerance; The procedure will find the maxmum r as dictionary size which satisfy: ![first equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Br%7D%5E%7B2%7D%7D%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Bn%7D%5E%7B2%7D%7D%5Cleq%20threshold)
   $$ sigma _{i} $$ 











