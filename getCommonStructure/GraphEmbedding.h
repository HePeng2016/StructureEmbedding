#include <armadillo>
#include <map>
#include <vector>
#include <list>
#include <set>

#include "mlpack/methods/sparse_coding/random_initializer.hpp"
#include "mlpack/methods/sparse_coding/sparse_coding.hpp"
using namespace mlpack::sparse_coding;





using namespace arma;
using namespace std;


typedef struct DFSC{
    unsigned int Rank;
    unsigned int index;
    unsigned int depth;
    std::vector <int> Projected;
}DFSC;




struct Compare
{
    bool operator()(const std::vector<int> &a, const std::vector<int> &b) const
    {

      if(a.size()<b.size())
      {
          return false;
      }
        if(a.size()>b.size())
      {
          return true;
      }
    for(int i=0;i<a.size();i++)
    {
              if (a[i] > b[i])
              {
                 return  true;
              }
              if(a[i] < b[i])
              {
                 return  false;
              }


    }
          return   false;
    }
};


struct DFSCPath: public std::vector <DFSC> {
public:
	void push (std::vector <int>&Projected,int index)
	{
		resize (size() + 1);
        DFSC &d = (*this)[size()-1];
        d.Projected = Projected;
		d.index = index;
	}
	void pop () { resize (size()-1); }
};//The  type  of deep first search


typedef struct clusterEntry
{
   int Rank;
   std::vector <int> IDArray;
   std::vector <int> PatternIDs;
}clusterEntry;


class GraphEncode{

     DFSCPath DFSCpath;
     std::set<int> FilterSetID;
     std::set < std::vector<int>,Compare>FilterSet;
     public:
     std::vector<clusterEntry> clusterResult;


     int MaxDepth = 1024;
     bool rankSkip = true;
     int minsupport=1;
     double tolerance;
     std::vector <int> Common_Set(std::vector <int>a,std::vector <int>b);
     void GamessFMOMatrixRead( FILE *  OneBody,FILE * TwoBody ,sp_mat&DistanceMatrix,sp_mat&EnergyMatrix);
     bool rankFilter(int Rank,std::vector <int> ReportPath,std::vector<int> parentSet);
     void K_GraphGeneration(sp_mat & DistanceMatrix,sp_mat & ConnectionMatrix);
     void Vectorization(sp_mat & ConnectionMatrix,int N, std::vector< std::vector< std::complex<double> > > &GraphEncodeList );


};
