#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <vector>
#include <algorithm>

#include "support.h"

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/

map<uint32_t, unsigned int> read_datafile(unsigned int *N, unsigned int n, string filename)    // O(N)  where N = data set size
{
  string line, line2;    
  (*N) = 0;            // N = dataset size
   Rcpp::Rcout << endl << "--->> Read \"" << filename << "\",\t Build Nset...";
  
  // ***** data are store in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set
  
  ifstream myfile (filename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);
      Nset[stoi(line2, 0, 2)] += 1;
      // Rcpp::Rcout << nb.to_ulong() << endl;   // Rcpp::Rcout << nb << " :  " << bitset<n>(nb) << endl;
      (*N)++;
    }
    myfile.close();
  }
  else  Rcpp::Rcout << "Unable to open file"; 
  
   Rcpp::Rcout << "\t\t data size N = " << (*N) << endl;
  
  return Nset;
}

//' Read data and store in Nset.
//' 
//' Reads the binary data set (with location and name provided in argument filename) 
//' as binary integers and returns a mapping of each observed state to the number of times 
//' they occur in the data set. Note that each state of the system is encoded as an n-bit integer.
//' 
//' @param N Integer that will contain the data set size.
//' @param n Number of binary (spin) variables of the data set.
//' @param filename The string containing the path to the data file and the name.
//' 
//' @return Nset The first row contains the observed state and the second row contains their occurrence (matching column-wise).
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<int>> read_datafile(Rcpp::IntegerVector N, unsigned int n, std::string filename)    // O(N)  where N = data set size
{
  string line, line2;     uint32_t nb = 0;
  N[0] = 0;            // N = dataset size
  Rcpp::Rcout << endl << "--->> Read \"" << filename << "\",\t Build Nset...";
  
  // ***** data are store in Nset:  *******************************
  vector<vector<int>> NsetVec; // = #of time state mu appears in the data set
  
  ifstream myfile (filename.c_str());
  if (myfile.is_open())
  {
    // Create Key, Value vectors!
    vector<int> key;
    vector<int> value;
    
    vector<int>::iterator itr;
    
    while (getline (myfile,line))
    { 
      line2 = line.substr (0,n); //take the n first characters of line
      nb = stoi(line2, 0, 2); //convert string line2 into a binary integer
        
      itr = std::find(key.begin(), key.end(), nb); // nb is the key.
      
      if (itr != key.cend()) // Key was present.
      { 
        int idx = std::distance(key.begin(), itr);
        value.at(idx) += 1; // Update value
      }
      else // New key-value pair.
      { 
        key.push_back(nb);
        value.push_back(1);
      }
      
      N[0]++;
    }
    
    // Create final vector!
    NsetVec.push_back(key);
    NsetVec.push_back(value);
    
    myfile.close();
  }
  else  Rcpp::Rcout << "Unable to open file";
  
   Rcpp::Rcout << "\t\t data size N = " << N[0] << endl;
  
  return NsetVec;
}

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a model (defined by the m basis vector) --> return the new m-state (state in the new m-basis)
// Rem: must have m <= n 
uint32_t transform_mu_basis(uint32_t mu, unsigned int n, list<uint32_t> basis)
{
  uint32_t bit_i = 1;
  uint32_t final_mu = 0;

  list<uint32_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    if ( (countSetBits((*phi_i) & mu) % 2) == 1) // odd number of 1, i.e. sig_i = 1
      {
        final_mu += bit_i;
      }
    bit_i = (bit_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/

// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:
map<uint32_t, unsigned int> build_Kset(map<uint32_t, unsigned int> Nset, list<uint32_t> Basis, unsigned int n, bool print_bool=false)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, unsigned int > Kset;

  uint32_t s;        // initial state
  uint32_t sig_m;    // transformed state and to the m first spins

  unsigned int ks=0; // number of time state s appear in the dataset

   Rcpp::Rcout << endl << "--->> Build Kset..." << endl;

  //Build Kset:
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    s = it->first;       // state s
    ks = it->second;    // # of times s appears in the data set
    sig_m = transform_mu_basis(s, n, Basis);
    //    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    
    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << " \t" << sig_m << ": \t" << int_to_bstring(sig_m, n) << endl; }

    Kset[sig_m] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
   Rcpp::Rcout << endl;

  return Kset;
}

//' Build K_set
//' 
//' Build Kset for the states written in the basis of the m-chosen independent 
//' operator on which the SC model is based. Changes the basis of the data-set from its original basis
//' to the basis provided as an argument in basis. Print Kset to terminal by setting print_bool to TRUE.
//' 
//' @param Nset A list containing 2 vectors which correspond to a index based mapping of the states in the data set to their occurrence.
//' @param Basis The basis on which the model is based. 
//' @param n Number of binary (spin) variables of the data set.
//' @param print_bool print relevant information. Default: FALSE.
//' 
//' @return Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<int>> build_Kset(std::vector<std::vector<int>> Nset, std::list<uint32_t> Basis, unsigned int n, bool print_bool=false) 
{
  vector<int> NsetKey = Nset.front();
  vector<int> NsetValue = Nset.back();

  std::vector<std::vector<int>> Kset_new_vec; // map<uint32_t, unsigned int > Kset;
  
  vector<int> KsetKey; // Create Key, Value vectors for the new Kset!
  vector<int> KsetValue;
  vector<int>::iterator itr;
  
  uint32_t s;        // initial state
  uint32_t sig_m;    // transformed state and to the m first spins
  
  unsigned int ks=0; // number of time state s appear in the dataset
  
   Rcpp::Rcout << endl << "--->> Build Kset..." << endl;
  
  //Build Kset:
  for(vector<int>::size_type i = 0; i != NsetKey.size(); i++) 
  {
    s = (uint32_t) NsetKey[i];       // state s
    ks = (uint32_t) NsetValue[i];    // # of times s appears in the data set
    sig_m = transform_mu_basis(s, n, Basis);

    if (print_bool)  {   Rcpp::Rcout << s << ": \t" << int_to_bstring(s, n) << " \t" << sig_m << ": \t" << int_to_bstring(sig_m, n) << endl; }
    
    // Build Kset vec.
    itr = std::find(KsetKey.begin(), KsetKey.end(), sig_m); // sig_m is the key.
    if (itr != KsetKey.cend()) // Key was present.
    { 
      int idx = std::distance(KsetKey.begin(), itr); 
      KsetValue.at(idx) += ks; // Update value.
    }
    else // New key-value pair.
    { 
      KsetKey.push_back(sig_m);
      KsetValue.push_back(ks);
    }
  }
  Kset_new_vec.push_back(KsetKey);
  Kset_new_vec.push_back(KsetValue);
  
   Rcpp::Rcout << endl;
  
  return Kset_new_vec; 
}
