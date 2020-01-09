#include "GraphEmbedding.h"


 void GraphEncode::GamessFMOMatrixRead( FILE *  OneBody,FILE * TwoBody ,sp_mat&DistanceMatrix,sp_mat&EnergyMatrix)
{

   char IgnorantWord [1024];
   char AminoAcidName[1024];
   double Value;



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
          EnergyMatrix(i,i)= OneBodyEnergys[i];
      }

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
                 if(Value>16384||Value==0||Value<-16384)
                 {

                  while(!feof(TwoBody))
                  {
		             char entry = fgetc(TwoBody);
                     if(entry=='\n')
		             {
		               break;
                     }
                  }
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

            EnergyMatrix(I-1,J-1)=Value;
            EnergyMatrix(J-1,I-1)=Value;

            if( Value >0.000000001|| Value < -0.000000001)
            {
               DistanceMatrix(I-1,J-1) = Distance+1.0;
               DistanceMatrix(J-1,I-1) = Distance+1.0;
            }

      }
      for(int i=0;i<DistanceMatrix.n_cols;i++)
      {
         DistanceMatrix(i,i) = 1.0;
      }



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
     for(int j=0;j< DistanceMatrix.n_cols;j++)
        for(int i=0;i< DistanceMatrix.n_rows;i++)
        {

             double x = DistanceMatrix(i,j);

             x=x-1;
             if( x>=0 )
             {
                x=x*x;
                MatrixX(i,j)= exp(-x);
             }else
             {
                MatrixX(i,j)=0;
             }

        }
      for(int j=0;j< DistanceMatrix.n_cols;j++)
      {
         double sum = -1;

         for(int i=0;i< DistanceMatrix.n_rows;i++)
         {
               sum = sum + MatrixX(i,j);
         };
         MatrixX(j,j)=sum;
      }

      for(int j=0;j< DistanceMatrix.n_cols;j++)
      {
         for(int i=0;i< DistanceMatrix.n_rows;i++)
         {
            if(i!=j)
            MatrixX(i,j)= MatrixX(i,j);
         };
       }



      double PC =0;// principal components
      double AC =0;// all  components
      int Rank =0;
      //MatrixX.save("DebugMatrix",csv_ascii);
      eig_sym(s,V,MatrixX);
      for(int i=0;i<V.n_cols;i++)
      MatrixX.col(i)= V.col(i)*sqrt(s[i]);



      svd_econ(U,s,V,MatrixX);
      for(int i=0;i<s.n_elem;i++)
      {
         AC = AC+s[i]*s[i];
      }

      for(int i=0;i<s.n_elem;i++)
      {
         PC = PC+s[i]*s[i];
         Rank ++;
         if(PC>((1.0-precision)*AC))
         {
           break;
         }
      }



     //MatrixX.save("DebugMatrix_2",csv_ascii);

     V.resize(MatrixX.n_cols,Rank);
     SparseCoding  sc(0,0,0,86,0.1,1e-1);
     sc.Lambda1() = 0.18;
     sc.Lambda2() = 0;
     sc.Atoms() = Rank;
     sc.MaxIterations()=1;
     sc.Dictionary() = V;
     //sc.Dictionary().save("disctionary1",csv_ascii);
     mat D;
     //printf("%d\n",sc.MaxIterations());
     sc.Encode(MatrixX.t(),D);
     ResultIDArray.resize(2*Rank);
     //D.save("debug",csv_ascii);
     //sc.Dictionary().save("disctionary2",csv_ascii);



     for(int i=0;i<Rank;i++)
     {
        for(int j=0;j<D.n_cols;j++)
        {
              if(D(i,j)>0.01)
              {
               ResultIDArray[i].resize(ResultIDArray[i].size()+1);
               ResultIDArray[i][ResultIDArray[i].size()-1] = j;
              }

             if(D(i,j)<-0.01)
             {
               ResultIDArray[i+Rank].resize(ResultIDArray[i+Rank].size()+1);
               ResultIDArray[i+Rank][ResultIDArray[i+Rank].size()-1] = j;
             }
        }
     }

    DFSCpath.resize(0);


    for(int i=0;i<ResultIDArray.size();i++)
    {
            DFSCpath.push(ResultIDArray[i],i);
            DFSCpath[DFSCpath.size()-1].Rank  = 1;
            DFSCpath[DFSCpath.size()-1].depth = 0;
    }

    {
         std::vector <int>reportPath;
         reportPath.resize(0);
         int Rank=1;

         while(DFSCpath.size())
        {
               DFSC * DFSset = &DFSCpath[DFSCpath.size()-1];
               std::vector<int> parentSet;
               std::vector<int> CoVector;


               int order = DFSset->index;
               int old_rank = DFSset->Rank;
               int old_depth = DFSset->depth;

               parentSet.resize(DFSset->Projected.size());
               std::copy(DFSset->Projected.begin(),DFSset->Projected.end(),parentSet.begin());
               reportPath.resize(DFSset->depth+1);
               reportPath[DFSset->depth]  = DFSset->index;
               DFSCpath.pop();
               /*Filter Code*/


               if(rankFilter(old_rank,reportPath,parentSet))
                {

                   continue;
                }

                if(reportPath.size()> MaxDepth)
                {
                    continue;
                }
               int MinRank = old_rank;
               bool Isminimum = true;
               int max_support = 0;
               int j;



               for(int i=0;i<order;i++)
               {


                  CoVector = Common_Set(ResultIDArray[i],parentSet);
                  if(rankSkip)
                  {
                     if(CoVector.size() < minsupport)
                     {
                         Rank = -1;

                     }else
                     {
                         Rank = 1;
                       }
                  }

                  if(max_support<CoVector.size())
                  {
                     max_support=CoVector.size();
                  }


                  if(Rank>0)
                  {
                        DFSCpath.push(CoVector,i);
                        DFSCpath[DFSCpath.size()-1].depth = old_depth+1;
                        DFSCpath[DFSCpath.size()-1].Rank = Rank;
                        if(MinRank >= Rank)
                        {
                           Isminimum = false;
                        }
                  }
               }


                if((max_support<parentSet.size()||order==0||Isminimum)&(!(rankFilter(old_rank,reportPath,parentSet))))
                {

                   FilterSet.insert(parentSet);

                   if(parentSet.size() >= minsupport)
                  {

                      bool Trim = false;

                      for(int i=0;i<parentSet.size()&&!Trim;i++)
                      {
                         for(int j=0;j<i;j++)
                         {

                              double  x = DistanceMatrix(parentSet[i],parentSet[j]) -1;

                               if (x<0)
                                x = 1024;

                            if( x > CutOff )
                            {
                               Trim = true;
                               break;
                            }
                         }
                      }
                    if(!Trim)
                    {
                       clusterEntry  Temp;
                       Temp.IDArray = reportPath;
                       Temp.Rank = old_rank;
                       Temp.PatternIDs = parentSet;
                       clusterResult.push_back(Temp);
                    }
                  }
               }
        }
   }

}



void GraphEncode::Vectorization(sp_mat & ConnectionMatrix, int N, std::vector< std::vector< std::complex<double> > > &GraphEncodeList )
{


    double Energy =0;
    std::complex<double> Image(0,1);
    double t_step= PI*2/(N);
    int PN;
    double Shift = round(N*ZeroPosition);

    if((Shift< Diffusion)||(Shift> N-Diffusion))
    {
      printf("Error: The ZeroPosition setting is improper\n");
      return;
    }
    double UpperEnerge = ((1-ZeroPosition)*N-Diffusion)*precision;
    double LowerEnerge = ((-ZeroPosition)*N+Diffusion)*precision;


    GraphEncodeList.resize(ConnectionMatrix.n_rows);

    for(PN=0;PN<N;PN++)
    {
      if( exp(-pow((PN*Diffusion/N),2)/2)<0.00001)
      break;
    }




    for(int i=0;i<(ConnectionMatrix.n_rows);i++)
    {
        GraphEncodeList[i].resize(PN);

        for(int j =0;j<PN;j++)
        {
          GraphEncodeList[i][j]=0.0;
        }

    }


    printf("%d/n",clusterResult.size());


    for(int i=0;i<clusterResult.size();i++)
    {

      int  EnergyEncodeType = 1;

      if( EnergyEncodeType == 1)
     {
     	   
     	   Energy=0; 
     	   
         for(int j =0;j<clusterResult[i].PatternIDs.size();j++)
        {
            for(int j2 =0;j2<=j;j2++)
          {
             Energy =  Energy+ConnectionMatrix(clusterResult[i].PatternIDs[j],clusterResult[i].PatternIDs[j2]);
          }
        }
          // printf("E:%lf\n",Energy);
          // printf("\n ");
          if((Energy<LowerEnerge)||(Energy>UpperEnerge))
          {
             printf("Warning: The energy %lf is discarded due to beyond abound, larger VectorSize is necessary or the zero position is needed to be reset.\n",Energy);
             continue; 
          }



         for(int j =0;j<clusterResult[i].PatternIDs.size();j++)
         {

            printf("%d ",clusterResult[i].PatternIDs[j]);
            for( int n=0;n<PN;n++)
              GraphEncodeList[clusterResult[i].PatternIDs[j]][n] = GraphEncodeList[clusterResult[i].PatternIDs[j]][n]+exp((-pow((n*Diffusion/N),2)/2))*exp((round(Energy/precision)+Shift)*Image*(t_step*n));

         }

            
      }

      if( EnergyEncodeType == 2)
      {

            for(int j =0;j<clusterResult[i].PatternIDs.size();j++)
           {

                Energy=0;


                for(int j2 =0;j2!=clusterResult[i].PatternIDs.size();j2++)
                {
                  Energy =  Energy+ConnectionMatrix(clusterResult[i].PatternIDs[j],clusterResult[i].PatternIDs[j2]);
                }

                 Energy =  Energy+ConnectionMatrix(clusterResult[i].PatternIDs[j],clusterResult[i].PatternIDs[j]);
               //  printf("E:%lf\n",Energy);
               // printf("\n ");
                 if((Energy<LowerEnerge)||(Energy>UpperEnerge))
                 {
                   printf("Warning: The energy %lf is discarded due to beyond abound, larger VectorSize is necessary or the zero position is neaded to be reset.\n",Energy);
                   continue;
                 }

                for( int n=0;n<PN;n++)
                 GraphEncodeList[clusterResult[i].PatternIDs[j]][n] = GraphEncodeList[clusterResult[i].PatternIDs[j]][n]+exp((-pow((n*Diffusion/N),2)/2))*exp((round(Energy/precision)+Shift)*Image*(t_step*n));
            }
      }

      double  GaussianBias;

      for(int i=0;i<(ConnectionMatrix.n_rows);i++)
      {
         GaussianBias = 0.0;
         for( int n=0;n<PN;n++)
         {
            GaussianBias = GaussianBias+real(GraphEncodeList[i][n]);
         }

         GraphEncodeList[i][0]= GraphEncodeList[i][0]-GaussianBias;
      }



 }
}














