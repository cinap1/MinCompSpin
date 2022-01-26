#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>       /* tgamma */
#include <map>
#include <vector>

#include "support.h"

using namespace std;

/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/

/********* of a Sub-Complete Model based on m basis operators     *************/
/******************************************************************************/
// for 1 <= m <= n . Rem C(m=1) = log(pi)
//' GeomComplexity_SubCM
//' 
//' Geometric Complexity of a Sub-Complete Model based on m basis operators.
//' for 1 <= m <= n . Rem C(m=1) = log(pi)
//' 
//' @param m The number of basis elements/operators.
//' 
//' @return GeomComplexity The complexity.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double GeomComplexity_SubCM(unsigned int m)     // Geometric complexity
{
  double pow = (double) ( 1UL << (m-1) );
  return (log(M_PI) * pow - lgamma(pow)  );   // lgamma(x) = log(gamma(x))
}

//' Parametric Complexity of a sub-complete model
//' 
//' Parametric Complexity of a Sub-Complete Model based on m basis operators.
//' for 1 <= m <= n . Rem C(m=1) = log(pi).
//' 
//' @param m The number of basis elements/operators.
//' @param N The size of the data-set.
//' 
//' @return ParamComplexity The complexity.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double ParamComplexity_SubCM(unsigned int m, unsigned int N)  // Parameter Complexity
{
  uint32_t K = (1UL << m) - 1;  // number of interactions
  return K * log(((double) N)/2./M_PI) / 2.;
}

/********* of the Minimally Complex Model (MCM) defined by "Partition"   ******/
/******************************************************************************/
// Compute separately: -- the first order complexity    --> stored in C_param
//                     -- and the geometric complexity  --> stored in C_geom

double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, double *C_param, double *C_geom)
{
  *C_param = 0;   *C_geom = 0;
  uint32_t m_i = 0;  // number of elements in Ai

  for (map<uint32_t, uint32_t>::iterator Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    m_i = countSetBits((*Part).second);
    (*C_param) += ParamComplexity_SubCM(m_i, N);
    (*C_geom) += GeomComplexity_SubCM(m_i);
  }  

  return (*C_param) + (*C_geom);
}

//' Complexity of an MCM based on the geometric and parametric complexity.
//' 
//' Complexity of and MCM defined by the addition of the total geometric complexity and
//' total (parametric) first order complexity. 
//' Compute separately: -- the first order complexity    --> stored in C_param
//'                     -- and the geometric complexity  --> stored in C_geom 
//' 
//' @param Partition Partition of MCM.
//' @param N The size of the dataset
//' @param n The amount of spin variables.
//' @param C_param The (parametric) first order complexity.
//' @param C_geom The geometric complexity.
//' 
//' @return Complexity A double containing the complexity of the MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double Complexity_MCM(std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, Rcpp::NumericVector C_param, Rcpp::NumericVector C_geom)
{
  vector<uint32_t> PartitionValue = Partition.back();

  C_param[0] = 0;   C_geom[0] = 0;
  uint32_t m_i = 0;  // number of elements in Ai
  
  for ( const auto &Part : PartitionValue )
  {
    m_i = countSetBits(Part);
    C_param[0] += ParamComplexity_SubCM(m_i, N);
    C_geom[0] += GeomComplexity_SubCM(m_i);
  }
  
  return C_param[0] + C_geom[0];
}
