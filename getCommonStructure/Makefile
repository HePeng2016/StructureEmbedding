objects = GraphEmbedding.o
CC = g++  -std=c++11  -g   #-Wall
LIB_FLAGS = -fopenmp  -lopenblas   -llapack  -lmlpack  -larmadillo
install : $(objects)
GraphEmbedding.o : GraphEmbedding.cpp   GraphEmbedding.h
        $(CC)$(LIB_FLAGS)  -c  $<   
clean :
        rm  $(objects) 
                     
