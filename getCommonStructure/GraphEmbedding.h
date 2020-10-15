//#include <armadillo>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <mlpack/core.hpp>
#include <mlpack/methods/sparse_coding/random_initializer.hpp>
#include <mlpack/methods/sparse_coding/sparse_coding.hpp>
#include <mlpack/methods/matrix_completion/matrix_completion.hpp>
#include <complex>
#include <cmath>
#define  PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706
using namespace mlpack;
using namespace mlpack::sparse_coding;
using namespace mlpack::matrix_completion;

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
     std::vector<std::vector<unsigned int > > ConnectedGraphs;  

     double CutOff = 3.3;
     int MaxDepth = 1024;
     bool rankSkip = true;
     int minsupport=1;
     double tolerance=0.001;
     double precision=0.01;
     double Diffusion=64;
     double ZeroPosition = 0.75;
     std::vector <int> Common_Set(std::vector <int>a,std::vector <int>b);
     void GamessFMOMatrixRead( FILE *  OneBody,FILE * TwoBody ,sp_mat&DistanceMatrix,mat&EnergyMatrix);
     bool rankFilter(int Rank,std::vector <int> ReportPath,std::vector<int> parentSet);
     void K_GraphGeneration(sp_mat & DistanceMatrix,sp_mat & ConnectionMatrix);
     void Vectorization(sp_mat & ConnectionMatrix,int N, std::vector< std::vector< std::complex<double> > > &GraphEncodeList );


};
