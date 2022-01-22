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
/**************** Log-likelihood (LogL) of a complete model  ******************/
/******************************************************************************/
// Compute the log-likelihood of a complete model on Kset:
// This function is mainly used for call by `LogL_SC_PartMCM`,
// but can also be used to compute the log-likelihood of a complete model
//
double LogL_CM(map<uint32_t, unsigned int > Kset, unsigned int N)
{
  double LogL = 0;

  map<uint32_t, unsigned int >::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  double Nd = N;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) { Rcpp::Rcout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogL += (Ks * log((double) Ks / Nd) );
  }
  if (Ncontrol != N) {  Rcpp::Rcout << "Error in function 'LogLikelihood_SCforMCM': Ncontrol != N" << endl;  }

  return LogL;
}

//' Compute the log-likelihood of a complete model on Kset.
//' 
//' Compute the log-likelihood of a complete model on Kset:
//' This function is mainly used for call by `LogL_SC_PartMCM`,
//' but can also be used to compute the log-likelihood of a complete model
//' 
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' 
//' @return LogL The log-likelihood of a complete model on Kset as a double.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double LogL_CM(std::vector<std::vector<int>> Kset, unsigned int N)
{
  vector<int> KsetKey = Kset.front();
  vector<int> KsetValue = Kset.back();
  
  double LogL = 0;

  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  double Nd = N;
  
  for(vector<int>::size_type i = 0; i != KsetValue.size(); i++) 
  {
    Ks = (uint32_t) KsetValue[i]; Ncontrol += Ks;
    if (Ks == 0) { Rcpp::Rcout << "problem Ks = 0 for mu_m = " << (KsetKey[i]) << endl; }
    LogL += (Ks * log((double) Ks / Nd) );
  }
  if (Ncontrol != N) {  Rcpp::Rcout << "Error in function 'LogLikelihood_SCforMCM': Ncontrol != N" << endl;  }

  return LogL;
}

/******************************************************************************/
/***************************** Log-likelihood (LogL) **************************/
/***********************   of a sub-complete part of a MCM   ******************/
/******************************************************************************/
// Compute the log-likelihood of the sub-complete part (of an MCM) defined by Ai.
// This function could be also used directly by the user
// to compute the log-likelihood of a sub-complete model

double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false)
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, unsigned int > Kset_new;

  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  if (print_bool)  { 
     Rcpp::Rcout << endl << "--->> Build Kset for SC Model based on "  << Ai << " = " << int_to_bstring(Ai, n) << " for MCM.." << endl;
  }
  //Build Kset_new:
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

  return LogL_CM(Kset_new, N);
}

//' Compute log-likelihood of the sub-complete part (SC-part) of an MCM.
//' 
//' Compute the log-likelihood of the sub-complete part (of an MCM) defined by Ai.
//' This function could be also used directly by the user to compute the 
//' log-likelihood of a sub-complete model
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param Ai The binary representation of the SC-part.
//' @param N The size of the data-set.
//' @param n Number of binary (spin) variables of the data set.
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return LogL The log-likelihood of a sub-complete part (of an MCM) defined by Ai.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double LogL_SubCM(std::vector<std::vector<int>> Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false)
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
    
    // Build new Kset vec.
    itr = std::find(key.begin(), key.end(), s); // Create iterator for vector in order to get the index.
    
    if (itr != key.cend()) { // Key was present so use index to adjust the value vector at the right index.
      int idx = std::distance(key.begin(), itr); // Determine index.
      value.at(idx) += ks; // Update value
    }
    else { // Key was not present so add new Key Value pair to their corresponding vectors.
      key.push_back(s);
      value.push_back(ks);
    }
  }
  
  Kset_new_vec.push_back(key);
  Kset_new_vec.push_back(value);
  
  if (print_bool)  {   Rcpp::Rcout << endl;  }
  
  return LogL_CM(Kset_new_vec, N);
}

/******************************************************************************/
/******************** Log-likelihood (LogL) of a MCM  *************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
bool check_partition(map<uint32_t, uint32_t> Partition, unsigned int n);
bool check_partition(std::vector<std::vector<uint32_t>> Partition, unsigned int n);

double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false)
{
  //if (!check_partition(Partition, n)) { Rcpp::Rcout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogL = 0; 
    unsigned int rank = 0;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogL += LogL_SubCM(Kset, (*Part).second, N, n);
      rank += countSetBits((*Part).second);
    }  
    return LogL - ((double) (N * (n-rank))) * log(2.);
  //}
}

//' Compute log-likelihood of an MCM.
//' 
//' Compute the log-likelihood of of an MCM.
//' This function could be also used directly by the user
//' to compute the log-likelihood of a model.
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param Partition Partition of MCM.
//' @param N The size of the data-set.
//' @param n Number of binary (spin) variables of the data set.
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return LogL The log-likelihood of a complete model on Kset as a double.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
double LogL_MCM(std::vector<std::vector<int>> Kset, std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, bool print_bool = false)
{
  //if (!check_partition(Partition, n)) { Rcpp::Rcout << "Error, the argument is not a partition." << endl; return 0;  }
  
  //else
  //{
  vector<uint32_t> PartitionValue = Partition.back();
  
  double LogL = 0; 
  unsigned int rank = 0;
  uint32_t Ai;
  
  for ( const auto &Part : PartitionValue )
  {
    Ai = Part;
    LogL += LogL_SubCM(Kset, Ai, N, n);
    rank += countSetBits(Ai);
  }
  
  return LogL - ((double) (N * (n-rank))) * log(2.);
  //}
}
