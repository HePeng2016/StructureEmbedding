The mlpack package should be installed firstly.
How to install mlpack can be found in address: https://github.com/mlpack/mlpack/.

# compile

     make

When 'make' command is fininshed, 'StructureEncode' file is generated. 

# Comands 


     ./StructureEncode Encode  DistanceAdjacencyMatrix.csv   EnergyAdjacencyMatrix.csv    OutputFile   [-config  configFile ]
     
  The DistanceAdjacencyMatrix is a csv file which is recorded the distance between fragments. The distance is the interfragment
 distance relative to van-der-Waals radii (-1.00 is printed if distances are not computed). 
  The size of distance adjacency matrix is  (column : K\*N , row : N), N is the number of fragments, and K is the number of conformations for a structure.  
  
                conformation.1              conformation.2               ...     conformation.K  
            frag.1  frag.2 ... frag.N   frag.1  frag.2 ... frag.N
    frag.1   0.0     0.76  ...  0.41     0.00    0.89  ...  0.42
    frag.2   0.76    0.00  ...  0.32     0.89    0.00  ...  0.33 
    ...                    ...
    frag.N   ...     ...                 ...                 ...
    

   ( Notice that the recorded value is the real value plus 1.0).  
   
   The EnergyAdjacencyMatrix is a csv file which is recorded the energy between fragments, and the energy for each fragment (all are in kcal/mol).  
   For multi-conformations, the energy of Energy Adjacency Matrix can be the average of conformations.
   
   The command will produce an Output File which record the embedding vectors for each fragment. If two fragments located in the similar  structural environment, then they will have similar vector (Cosine similar metric). 
   
   
    ./StructureEncode  FmoToM   OneBody.txt   TwoBody.txt DistanceAdjacencyMatrix.csv   EnergyAdjacencyMatrix.csv  
   
   OneBody.txt and TwoBody.txt are two parts of result of GAMESS fragment molecular orbital analysis,   

# config file

     Tolerance = 0.01
     VectorSize = 128
     
Tolerance : 

   The Tolerance is used to determine the size of represent dictionary for distance adjacency matrix. If Tolerance is specified, threshold = 1-Tolerance; The procedure will find the maxmum r as dictionary size which satisfy: ![first equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Br%7D%5E%7B2%7D%7D%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Bn%7D%5E%7B2%7D%7D%5Cleq%20threshold)
   
 { ![first equation](http://latex.codecogs.com/gif.latex?%7B%5Csigma%20%7D_j) } are  singular values for distance adjacency matrix using SVD.
 
 VectorSize :
 
   VectorSize is the size of final embedding vector for each fragement in a structure.
  
If the config file is missing, the default tolerance is 0.01 and the default VectorSize is 128.
 












