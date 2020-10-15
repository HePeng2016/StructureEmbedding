objects =getCommonStructure/GraphEmbedding.o
CXXFLAGS =  -I    MlpackInclude
LIB_FLAGS =  -lopenblas    -llapack
CC = g++  -std=c++11  -fopenmp     -g  -I getCommonStructure/    #-Wall
all:    prepareall
        $(CC)  -o   StructureEncode    main.cpp  $(objects) $(LIB_FLAGS)   -lmlpack  -larmadillo    
prepareall:    subsystem
subsystem:
        $(MAKE) -C getCommonStructure/
clean :  cleansub
        rm   StructureEncode
cleansub :
        rm  $(objects)


