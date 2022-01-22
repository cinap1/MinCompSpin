#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>       /* tgamma */
#include <map>
#include <iostream>
#include <vector>
#include <algorithm>

#include "support.h"

using namespace std;

/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/
double GeomComplexity_SubCM(unsigned int m);

/******************************************************************************/
/**************** Log-Evidence (LogE) of a sub-complete model  ****************/
/******************************************************************************/
// Compute the log-evidence of a sub-complete model based on m basis elements
// ! Kset must have been previously reduced to these m basis elements !
// This function is mainly used for call by `LogE_PartMCM`,
// but can also be used to compute the log-likelihood of a complete model
//
double LogE_SubC_forMCM(map<uint32_t, unsigned int > Kset, uint32_t m, unsigned int N)
{
  double LogE = 0;
  
  map<uint32_t, unsigned int >::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) { Rcpp::Rcout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogE += lgamma(Ks + 0.5);
  }
  if (Ncontrol != N) {  Rcpp::Rcout << "Error Likelihood function: Ncontrol != N" << endl;  }
  
  return LogE - GeomComplexity_SubCM(m) - lgamma( (double)( N + (1UL << (m-1)) ) );
}

double LogE_SubC_forMCM(std::vector<std::vector<int>> Kset, uint32_t m, unsigned int N)
{
  vector<int> KsetKey = Kset.front();
  vector<int> KsetValue = Kset.back();
  
  double LogE = 0;
  
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  
  for(vector<int>::size_type i = 0; i != KsetValue.size(); i++) 
  {
    Ks = (uint32_t) KsetValue[i]; Ncontrol += Ks;
    if (Ks == 0) { Rcpp::Rcout << "problem Ks = 0 for mu_m = " << (KsetKey[i]) << endl; }
    LogE += lgamma(Ks + 0.5);
  }
  if (Ncontrol != N) {  Rcpp::Rcout << "Error Likelihood function: Ncontrol != N" << endl;  }
  
  return LogE - GeomComplexity_SubCM(m) - lgamma( (double)( N + (1UL << (m-1)) ) );
}

/******************************************************************************/
/*********  Log-Evidence (LogE) of a sub-complete part of a MCM   *************/
/******************************************************************************/
// Compute the log-evidence of the sub-complete part (of an MCM) defined by Ai.
// This function could be also used directly by the user
// to compute the log-likelihood of a sub-complete model

double LogE_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false)
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, unsigned int > Kset_new;
  
  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset
  
  if (print_bool)  { 
     Rcpp::Rcout << endl << "--->> Build Kset for SC Model based on "  << Ai << " = " << int_to_bstring(Ai, n) << " for MCM.." << endl;
  }
  //Build Kset:
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    s = it->first;      // initial state s 
    ks = it->second;    // # of times s appears in the data set
    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << " \t" ;  }
    
    s &= Ai;   // troncated state: take only the bits indicated by Ai
    //    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << endl; }
    
    Kset_new[s] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  if (print_bool)  {   Rcpp::Rcout << endl;  }
  
  return LogE_SubC_forMCM(Kset_new, countSetBits(Ai), N);
}


//' Compute the log-evidence of the sub-complete part (of an MCM).
//' 
//' Compute the log-evidence of the sub-complete part (of an MCM) defined by Ai.
//' This function could be also used directly by the user
//' to compute the log-likelihood of a sub-complete model
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param Ai The binary representation of the sub-complete part (of an MCM).
//' @param N The size of the data-set.
//' @param n The amount of binary (spin) variables.
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return LogE The log-evidence of the sub-complete part (of an MCM) as a double.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double LogE_SubCM(std::vector<std::vector<int>> Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false)
{
  vector<int> KsetKey = Kset.front();
  vector<int> KsetValue = Kset.back();
  
  // Create Key, Value vectors for the new Kset!
  vector<int> key;
  vector<int> value;
  
  vector<int>::iterator itr;
  std::vector<std::vector<int>> Kset_new_vec;
  
  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset
  
  if (print_bool)  { 
     Rcpp::Rcout << endl << "--->> Build Kset for SC Model based on "  << Ai << " = " << int_to_bstring(Ai, n) << " for MCM.." << endl;
  }
  
  // Build Kset_new vector based
  for(vector<int>::size_type i = 0; i != KsetKey.size(); i++) 
  {
    s = (uint32_t) KsetKey[i];       // state s
    ks = (uint32_t) KsetValue[i];    // # of times s appears in the data set
    
    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << " \t" ;  }
    
    s &= Ai;   // troncated state: take only the bits indicated by Ai
    
    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << endl; }
    
    itr = std::find(key.begin(), key.end(), s); // The key is s.
    if (itr != key.cend()) // Key was present.
    { 
      int idx = std::distance(key.begin(), itr); 
      value.at(idx) += ks; // Update value
    }
    else 
    { 
      key.push_back(s);
      value.push_back(ks);
    }
  }
  
  Kset_new_vec.push_back(key);
  Kset_new_vec.push_back(value);
  
  if (print_bool)  {   Rcpp::Rcout << endl;  }
  
  return LogE_SubC_forMCM(Kset_new_vec, countSetBits(Ai), N);
}

/******************************************************************************/
/****************************   LogE of a MCM   *******************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
bool check_partition(map<uint32_t, uint32_t> Partition);

double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false)
{
  //if (!check_partition(Partition)) { Rcpp::Rcout << "Error, the argument is not a partition." << endl; return 0;  }
  
  //else
  //{
  double LogE = 0; 
  unsigned int rank = 0;
  map<uint32_t, uint32_t>::iterator Part;
  
  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    LogE += LogE_SubCM(Kset, (*Part).second, N, n);
    rank += countSetBits((*Part).second);
  }  
  return LogE - ((double) (N * (n-rank))) * log(2.);
  //}
}

//' Compute the log-evidence of of an MCM.
//' 
//' Compute the log-evidence of of an MCM.
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param Partition Partition of MCM.
//' @param N The size of the data-set.
//' @param n The amount of binary (spin) variables.
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return LogE The log-evidence of an MCM as a double.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double LogE_MCM(std::vector<std::vector<int>> Kset, std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, bool print_bool = false)
{
  vector<uint32_t> PartitionValue = Partition.back();
  
  double LogE = 0; 
  unsigned int rank = 0;
  uint32_t Ai;
  
  for ( const auto &Part : PartitionValue )
  {
    Ai = Part;
    LogE += LogE_SubCM(Kset, Ai, N, n);
    rank += countSetBits(Ai);
  }
  
  return LogE - ((double) (N * (n-rank))) * log(2.);
}
