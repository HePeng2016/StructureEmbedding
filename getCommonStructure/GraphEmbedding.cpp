#include "GraphEmbedding.h"




 void GraphEncode::GamessFMOMatrixRead( FILE *  OneBody,FILE * TwoBody ,sp_mat&DistanceMatrix, mat&RecEnergyMatrix)
{

   char IgnorantWord [1024];
   char AminoAcidName[1024];
   double Value;
   unsigned IndicesSize =0;
   mat EnergyMatrix;


   std::vector <std::string> AminoAcidNames;
   std::vector <double> OneBodyEnergys;

       while(!feof(OneBody))
	  {

               fscanf(OneBody,"{\b|\t|\n}*",IgnorantWord);

               strcpy(AminoAcidName,"");

               if( fscanf(OneBody,"%s",AminoAcidName) )
               {
                     fscanf(OneBody,"{\b|\t|\n}*",IgnorantWord);
                     fscanf(OneBody,"%s",IgnorantWord);

               }else
               {
                   break;
               }


               fscanf(OneBody,"{\b|\t|\n}*",IgnorantWord);



               if( fscanf(OneBody,"%lf",&Value) )
               {

               }else
               {
                   break;
               }


                while(!feof(OneBody))
               {
		          char entry = fgetc(OneBody);
                  if(entry=='\n')
		          {
		            break;
		           }
                }
              if(strcmp(AminoAcidName,"") !=0)
              {
                AminoAcidNames.resize(AminoAcidNames.size()+1);
                AminoAcidNames[AminoAcidNames.size()-1].assign(AminoAcidName);

                OneBodyEnergys.resize(OneBodyEnergys.size()+1);
                OneBodyEnergys[OneBodyEnergys.size()-1]=Value;
              }
	  }



      EnergyMatrix.resize(OneBodyEnergys.size(),OneBodyEnergys.size());
      DistanceMatrix.resize(OneBodyEnergys.size(),OneBodyEnergys.size());

      for(int i=0;i< EnergyMatrix.n_rows;i++)
      {
          EnergyMatrix(i,i)= OneBodyEnergys[i];
      }
      IndicesSize = IndicesSize+EnergyMatrix.n_rows;
      while(!feof(TwoBody))
      {
          int I;
          int J;
          int Index;
          double Energy;
          double Distance;

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);


               fscanf(TwoBody,"    ",IgnorantWord);

               if( fscanf(TwoBody,"%d",&Index) )   //1
               {
                   I=Index;
               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%d",&Index) )  //2
               {
                   J=Index;
               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%s",IgnorantWord))  //3
               {

               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%lf",&Value) )  //4
               {

               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%lf",&Value) ) //5
               {
                   Distance = Value;
                   if(  Distance > -0.000000001)
                   {
                      DistanceMatrix(I-1,J-1) = Distance+1.0;
                      DistanceMatrix(J-1,I-1) = Distance+1.0;
                   }
               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%lf",&Value) ) //6
               {

               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);

               if( fscanf(TwoBody,"%lf",&Value) ) //7
               {
                 if(Value>16384||Value<-16384)
                 {

                  while(!feof(TwoBody))
                  {
		             char entry = fgetc(TwoBody);
                     if(entry=='\n')
		             {
		               break;
                     }
                  }

                   EnergyMatrix(I-1,J-1)=Value;
                   EnergyMatrix(J-1,I-1)=Value;


                   continue;
                 }

               }else
               {  break; }

               fscanf(TwoBody,"{\b|\t|\n}*",IgnorantWord);


               while(!feof(TwoBody))
                {
		             char entry = fgetc(TwoBody);
                     if(entry=='\n')
		             {
		               break;
                     }
                }



            IndicesSize=IndicesSize+2;
            EnergyMatrix(I-1,J-1)=Value;
            EnergyMatrix(J-1,I-1)=Value;

/*            if(  Distance > -0.000000001)
            {
               DistanceMatrix(I-1,J-1) = Distance+1.0;
               DistanceMatrix(J-1,I-1) = Distance+1.0;
            } */

      }
      for(int i=0;i<DistanceMatrix.n_cols;i++)
      {
         DistanceMatrix(i,i) = 1.0;
      }

    {

        arma::umat indices;  // contains the known indices [2 x n_entries]
 //     arma::vec values;    // contains the known values [n_entries]
        indices.resize(2, EnergyMatrix.n_rows*EnergyMatrix.n_rows - IndicesSize);
//        values.set_size(IndicesSize);
        unsigned int n_Index;
        n_Index = 0;
        for(int i=0;i<EnergyMatrix.n_rows;i++)
        {

            for( int j=0;j<=i;j++)
            {
              if( !( (EnergyMatrix(i,j)>= -16384)&&(EnergyMatrix(i,j) <= 16384) ))
              {
                indices(0,n_Index)=i;
                indices(1,n_Index)=j;
             //   values(n_Index)=EnergyMatrix(i,j);
                n_Index++;
                if(i!=j)
                { indices(1,n_Index)=i;
                  indices(0,n_Index)=j;
             //   values(n_Index)=EnergyMatrix(i,j);
                  n_Index++;
                }
              }
            }
        }

        std::vector <unsigned int> subID;
        subID.resize(DistanceMatrix.n_cols);

        for(int i=0;i<indices.n_cols;i++)
        {
          int I = indices(1,i);
          int J = indices(0,i);

          if(I!=J)
          {
             i++;
          }

          int sub_Index = 0;
          int key_Index1;
          int key_Index2;
          for(int i=0;i<DistanceMatrix.n_cols;i++)
          {
              if(DistanceMatrix(I,i)<3||DistanceMatrix(J,i)<3)
              {
                 subID[sub_Index]=i;
                 if((i==I))
                 {
                    key_Index1 = sub_Index;
                 }
                 if((i==J))
                 {
                    key_Index2 = sub_Index;
                 }
                 sub_Index ++;
              }
          }

          arma::umat indices_;
          indices_.resize(sub_Index);

          for(int i=0;i<sub_Index;i++)
          {
            indices_[i] = subID[i];
          }


          mat SubEnergyMatrix = EnergyMatrix.submat(indices_,indices_);
          arma::umat SubIndices;
          arma::mat SubValues;


          SubIndices.resize(2,sub_Index*sub_Index);
          SubValues.resize(SubIndices.n_cols);
          n_Index=0;

          for(int i=0;i<SubEnergyMatrix.n_rows;i++)
          {

             for( int j=0;j<=i;j++)
             {
                if(((SubEnergyMatrix(i,j)>= -16384)&&(SubEnergyMatrix(i,j) <= 16384)))
               {
                 SubIndices(0,n_Index)=i;
                 SubIndices(1,n_Index)=j;
                 SubValues(n_Index)=SubEnergyMatrix(i,j);
                 n_Index++;
                if(i!=j)
                { SubIndices(1,n_Index)=i;
                  SubIndices(0,n_Index)=j;
                  SubValues(n_Index)=SubEnergyMatrix(i,j);
                  n_Index++;
                }
               }
            }
         }

        SubValues.resize(n_Index);
        SubIndices.resize(2,n_Index);


         {
         	 arma::mat recovered;
         	 MatrixCompletion mc( SubEnergyMatrix.n_rows,SubEnergyMatrix.n_cols,SubIndices,SubValues);
	         mc.Recover(recovered);
             EnergyMatrix(I,J)= recovered(key_Index1,key_Index2);
             EnergyMatrix(J,I)= recovered(key_Index1,key_Index2);
	     }
    }
  }

    RecEnergyMatrix = EnergyMatrix;
}




void  laplaceNormalized(sp_mat & AdjacencyM,sp_mat & LaplaceM )
{

     std:vector <double> D;

     D.resize(AdjacencyM.n_rows);
     LaplaceM.resize( AdjacencyM.n_rows, AdjacencyM.n_cols);
     for(int i=0;i<AdjacencyM.n_rows;i++)
     {
       D[i]=0;
     }

     for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
        sp_mat::row_iterator it = AdjacencyM.begin_row(i);
        sp_mat::row_iterator it_end = AdjacencyM.end_row(i);

        for(;it !=it_end;++it)
	    {
          D[i] = D[i]+(*it);
		}
	  }

      for ( int i=0;i<AdjacencyM.n_rows;i++)
	  {
	     D[i] = sqrt(D[i]);
	  }

	  for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
	    sp_mat::row_iterator   it = AdjacencyM.begin_row(i);
		sp_mat::row_iterator   it_end = AdjacencyM.end_row(i);

		for(;it !=it_end;++it)
		{
			int row;
			int col;

			row= it.row();
			col= it.col();
			double value = (*it);
            LaplaceM.at(row,col)=-value/(D[row]*D[col]);
		}
	  }
	  for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
          LaplaceM.at(i,i) = 1+LaplaceM.at(i,i);
      }
}


void  laplace(sp_mat & AdjacencyM,sp_mat & LaplaceM )
{

    std:vector <double> D;

    D.resize(AdjacencyM.n_rows);
    LaplaceM.resize( AdjacencyM.n_rows, AdjacencyM.n_cols);
    for(int i=0;i<AdjacencyM.n_rows;i++)
    {
        D[i]=0;
    }

    for ( int i=0;i<AdjacencyM.n_rows;i++)
    {

       sp_mat::row_iterator  it = AdjacencyM.begin_row(i);
       sp_mat::row_iterator  it_end = AdjacencyM.end_row(i);

        for(;it !=it_end;++it)
	    {
          D[i] = D[i]+(*it);
		}
    }
     for ( int i=0;i<AdjacencyM.n_rows;i++)
     {
	    sp_mat::row_iterator it = AdjacencyM.begin_row(i);
		sp_mat::row_iterator it_end = AdjacencyM.end_row(i);

		for(;it !=it_end;++it)
		{
			int row;
			int col;

			row= it.row();
			col= it.col();
			double value = (*it);
            LaplaceM.at(row,col) = -value;
		}
	  }
     for ( int i=0;i<AdjacencyM.n_rows;i++ )
     {
         LaplaceM.at(i,i) = D[i]+LaplaceM.at(i,i);
     }
}


void ConnectionMatrixConstruction(sp_mat & ConnectionMatrix, double NFactor,sp_cx_mat & OutPutM)
{
    std::complex<double> I(0,1);
    OutPutM.resize(ConnectionMatrix.n_rows,ConnectionMatrix.n_cols);


       for ( int i=0;i<ConnectionMatrix.n_rows;i++)
       {
        for  ( int j=0;j<ConnectionMatrix.n_cols;j++)
         {
            if (ConnectionMatrix(i,j)!=0)
            {
              OutPutM(i,j)=exp((ConnectionMatrix(i,j)/(NFactor*2) )*I);
			};
         }
       }
}



std::vector <int> GraphEncode::Common_Set(std::vector <int>a,std::vector <int>b)
{

        std::vector <int> CoVector;

        for(int i=0,j=0;i<a.size()&&j<b.size();)
        {
            if(a[i]==b[j])
            {

               CoVector.resize(CoVector.size()+1);
               CoVector[CoVector.size()-1]= a[i];
               i++;
               j++;
               continue;
            }
            if(a[i]>b[j])
            {
                j++;
            }else
            {
                i++;
            }

        }
        return  CoVector;
}



bool GraphEncode::rankFilter(int Rank, std::vector <int> ReportPath,std::vector<int> parentSet)
{

       if((FilterSet.find(parentSet) == FilterSet.end()))
       {
           if(FilterSetID.size()!=0)
           {
              for(int i=0;i<parentSet.size();i++)
             {
                if( FilterSetID.find(parentSet[i])!= FilterSetID.end() )
                {
                   return false;
                }
             }
                  return true;
           }

            return false;
       }else
         {
            return true;
         }
}





void GraphEncode::K_GraphGeneration(sp_mat & DistanceMatrix,sp_mat & ConnectionMatrix)
{
     mat MatrixX;
     vec s;
     mat V;
     mat U;
     std::vector< std::vector<int> >ResultIDArray;
     double value1;
     ResultIDArray.resize(DistanceMatrix.n_rows);
     MatrixX.resize( DistanceMatrix.n_rows, DistanceMatrix.n_cols);


     if(DistanceMatrix.n_cols!=DistanceMatrix.n_rows)
     {
        printf("The Distance Matrix should be a square matrix\n");
     }


   {

        std::vector <unsigned int> subID;
        subID.resize(DistanceMatrix.n_cols);
        unsigned int Stack[128];
        int stackI = 0;
        int MaxStackSize=5;
        int count__ = 0;
        int K = 1;
	double CutOFF;
	CutOFF = 3.8;


        for(int I=0;I<DistanceMatrix.n_cols;I++)
        {

          int sub_Index = 0;
          stackI = 0;

          for(int i=I+1;i<DistanceMatrix.n_cols;i++)
          {
              if(abs(ConnectionMatrix(I,i))>CutOFF)
              {
                 subID[sub_Index]=i;
                 sub_Index ++;
              }
          }

          if(sub_Index==0)
            continue;

          Stack[stackI]=sub_Index;
          MaxStackSize =K-1;
         if(MaxStackSize>=0)
		 {
		   while(stackI>=0)
           {

              bool IsConnected =true;
              Stack[stackI]=Stack[stackI]-1;


              {
                 std::vector <unsigned int>key;
                 key.resize(stackI+2);
                 for(int i=0;i<stackI;i++)
                 {
                   if( abs(ConnectionMatrix(subID[Stack[i]],subID[Stack[stackI]])) > CutOFF)
                   {
                     IsConnected = false;
                     break;
                   }
                 }
                 if(IsConnected)
                 {
                    for(int i=0;i<=stackI;i++)
                    key[i] = subID[Stack[i]];
                    key[stackI+1] = I;
                    ConnectedGraphs.push_back(key);
                  }

              }


                if((stackI<MaxStackSize)&&(Stack[stackI]!=0)&&(IsConnected))
                {
                   stackI=stackI+1;
                   Stack[stackI]=Stack[stackI-1];
                }else if(Stack[stackI]==0)
                {
                   stackI=stackI-1;
                }
             };
		}
		}
      printf("AtomicComponents:%d\n",ConnectedGraphs.size());
       }

        for(int i=0;i<ConnectionMatrix.n_cols;i++)
		{
		   std::vector <unsigned int>key;
           key.resize(1);
		   key[0]=i;
           ConnectedGraphs.push_back(key);
		}

        //short range
/*       ConnectedGraphs.resize(0);

       {
	   std::vector <unsigned int>key;

           for(int I=0;I<DistanceMatrix.n_cols;I++)
           {
                for(int i=0;i<I;i++)
               {
                   if((DistanceMatrix(I,i)>=3))
                   {
                      key.resize(2);
                      key[0] = I;
                      key[1] = i;
                      ConnectedGraphs.push_back(key);
                   }
               }
           }
       }*/
   //long range
     printf("AllComponents:%d\n",ConnectedGraphs.size());
     return ;
}



void GraphEncode::Vectorization(sp_mat & EnergyMatrix, int N, std::vector< std::vector< std::complex<double> > > &GraphEncodeList )
{


    double Energy =0;
    double AllEnergy=0;
    std::complex<double> Image(0,1);
    double t_step= PI*2/(N);
    int PN;
    double Shift = round(N*ZeroPosition);
    int K; 
    size_t n_Index;         // size of unknown matrix
    double SingleAminoAcidEnerge=-12;
    double Ratio;
    std::vector < std::vector <int> > subID; 
    subID.resize(EnergyMatrix.n_rows);
    K = 64;
    std::vector<int> Flags;
	double CutOFF=0.1;


    for(int i=0;i<EnergyMatrix.n_rows;i++)
	{
       for(int j=0;j<EnergyMatrix.n_cols;j++)
       {
          if(abs( EnergyMatrix(i,j))>CutOFF&&(i!=j))
           {
			   subID[i].resize(subID[i].size()+1);	
			   subID[i][subID[i].size()-1]=j;
		  }
       }
      std:sort(subID[i].begin(),subID[i].end());
	}



    if((Shift< Diffusion)||(Shift> N-Diffusion))
    {
      printf("Error: The ZeroPosition setting is improper\n");
      return;
    }
    double UpperEnerge = ((1-ZeroPosition)*N-Diffusion)*precision;
    double LowerEnerge = ((-ZeroPosition)*N+Diffusion)*precision;


    GraphEncodeList.resize(EnergyMatrix.n_rows);

    for(PN=0;PN<N;PN++)
    {
      if( exp(-pow((PN*Diffusion/N),2)/2)<0.00001)
      break;
    }




    for(int i=0;i<(EnergyMatrix.n_rows);i++)
    {
        GraphEncodeList[i].resize(PN);

        for(int j =0;j<PN;j++)
        {
          GraphEncodeList[i][j]=0.0;
        }

    }


   // printf("%d\n",clusterResult.size()); 

	Flags.resize(64);


    for(int i=0;i<ConnectedGraphs.size();i++)
    {

      {

         AllEnergy=0;
         Energy = 0;

         for(int j =0;j<ConnectedGraphs[i].size();j++)
         {
              //for(int j2 =0;j2<=j;j2++)
             {
                AllEnergy =  AllEnergy+EnergyMatrix(ConnectedGraphs[i][j],ConnectedGraphs[i][j]);
             }
              Flags[j]=0;
          }
     }

	     //
         // double Ratio =  AllEnergy/Energy-1;


          if((Energy<LowerEnerge)||(Energy>UpperEnerge))
          {
              printf("Warning: The energy %lf is discarded due to beyond abound, larger VectorSize is necessary or the zero position is needed to be reset.\n",Energy);
	      continue;
          }

                  int MAX = -1;
		  bool BREAK;
                  bool GlobleBreak = false; 
		  bool skip;

		  while(!GlobleBreak)
		  {
	                   skip = false; 
			   MAX = MAX+1;

			   do{

				        BREAK = true;

				        for(int j =0;j<ConnectedGraphs[i].size();j++)
					{
		  
					   if( subID[ConnectedGraphs[i][j]].size() <= Flags[j])
					   {
						 GlobleBreak = true;
						 BREAK = true;
						 break;
					    }

					  if( subID[ConnectedGraphs[i][j]][Flags[j]] < MAX )
					  {
						   Flags[j]= Flags[j]+1;
						   BREAK = false;
					   } else if(MAX != subID[ConnectedGraphs[i][j]][Flags[j]] )
					   {
						   MAX = subID[ConnectedGraphs[i][j]][Flags[j]];
						   BREAK = false;
					   }
				        }

				}while(!BREAK);


           
				if( GlobleBreak == true)
				{
				    continue;
				}


				{
					double  SubEnergy = 0;

					for(int j2=0;j2<ConnectedGraphs[i].size();j2++)
					{
						
					   SubEnergy = SubEnergy + EnergyMatrix(MAX,ConnectedGraphs[i][j2]);
					
					   if(j2==MAX)
					   {
					     skip = true; 
					   }
					}

				   Ratio = SubEnergy/(AllEnergy+EnergyMatrix(MAX,MAX));

				if( Ratio<0.01||skip==true)
				continue;

				Energy = AllEnergy + EnergyMatrix(MAX,MAX) - (ConnectedGraphs[i].size()+1)*SingleAminoAcidEnerge;

			   for( int n=0;n<PN;n++)
                           GraphEncodeList[MAX][n] = GraphEncodeList[MAX][n]+(Ratio)*exp((-pow((n*Diffusion/N),2)/2))*exp(-(round(Energy/precision)+Shift)*Image*(t_step*n));
              
			    printf("Energy/AllEnergy	%lf  ",(Ratio));
                            printf(" %d  ",ConnectedGraphs[i].size());
	                    printf(" %lf \n",Energy);

				}

		  }
		   
      }

 }




