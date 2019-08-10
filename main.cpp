#include  "getCommonStructure/GraphEmbedding.h"
using namespace std;



int main( int argc , char *argv[] )
{

     GraphEncode  Test;
     Test.tolerance=0.01;
     double K =128;
     FILE * configip=NULL;
     char buffer[1024];


    for(int i =0;i<argc-1;i++)
    {
       if( strcmp(argv[i],"-config")==0)
       {
           configip = fopen(argv[i+1],"r");
       }
    }

   if( configip == NULL)
    {
        configip = fopen("config","r");
    }

   if(configip != NULL )
   {
     while (!feof(configip))
    {
        fscanf(configip,"%[\b|\t]*",buffer);
        fscanf(configip,"%s",buffer);

        if(strcmp(buffer,"Tolerance")==0)
        {
            fscanf(configip,"%[\b|\t]*",buffer);
            fscanf(configip,"%s",buffer);

            if(strcmp(buffer,"=")==0)
            {
               fscanf(configip,"%lf",&Test.tolerance);
            }
               fscanf(configip,"%[\b|\t|\n]",buffer);
         continue;
       }

       if(strcmp(buffer,"VectorSize")==0)
       {
            fscanf(configip,"%[\b|\t]*",buffer);
            fscanf(configip,"%s",buffer);

            if(strcmp(buffer,"=")==0)
            {
               fscanf(configip,"%lf",&K);
            }
               fscanf(configip,"%[\b|\t|\n]",buffer);
         continue;
       }
    }
   }




for(int i =0;i<argc-1;i++)
{

  if( strcmp(argv[i],"FmoToM")==0)
  {
     FILE * ip1 = fopen(argv[i+1],"r");
     FILE * ip2 = fopen(argv[i+2],"r");
     sp_mat DistanceMatrix;
     sp_mat EnergyMatrix;
     sp_cx_mat  OutPutM;
     Test.GamessFMOMatrixRead(ip1,ip2,DistanceMatrix,EnergyMatrix);
     mat B = (mat)DistanceMatrix;
     B.save(argv[i+3],csv_ascii);
     mat C = (mat)EnergyMatrix;
     C.save(argv[i+4],csv_ascii);

  }
  if ( strcmp(argv[i],"Encode") ==0 )
  {

        mat A;
        A.load(argv[i+1],csv_ascii);
        sp_mat  DistanceMatrix = (sp_mat)A;
        mat B;
        B.load(argv[i+2],csv_ascii);
        sp_mat EnergyMatrix = (sp_mat)B;
        FILE * ip3 = fopen(argv[i+3],"w+");



        Test.minsupport = 1;
        std::vector< std::vector< std::complex<double> > > GraphEncodeList;
        Test.K_GraphGeneration(DistanceMatrix,EnergyMatrix);
        Test.Vectorization(EnergyMatrix,K,GraphEncodeList);


     for(int  i=0;i<GraphEncodeList.size();i++)
     {

          for(int j =0;j<GraphEncodeList[i].size();j++)
          {
             fprintf(ip3,"%lf	", imag(GraphEncodeList[i][j]));
          }
          fprintf(ip3,"\n");
      }
     fclose(ip3);
    }

  }
}



