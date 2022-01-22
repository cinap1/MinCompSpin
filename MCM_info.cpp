#include <Rcpp.h>
using namespace Rcpp;

#include <map>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

#include "support.h"

using namespace std;

/******************************************************************************/
/**************** Log-likelihood (LogL), Geometric Complexity *****************/
/*************************  and Log-evidence (LogE) ***************************/
/******************************************************************************/
double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false);
double LogL_MCM(std::vector<std::vector<int>> Kset, std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, bool print_bool = false);
  
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false);
double LogE_MCM(std::vector<std::vector<int>> Kset, std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, bool print_bool = false);

double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, double *C_param, double *C_geom);
double Complexity_MCM(std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, Rcpp::NumericVector C_param, Rcpp::NumericVector C_geom);

double LogE_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);
double LogE_SubCM(std::vector<std::vector<int>> Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);

double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);
double LogL_SubCM(std::vector<std::vector<int>> Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);

double GeomComplexity_SubCM(unsigned int m);
double ParamComplexity_SubCM(unsigned int m, unsigned int N);


/******************************************************************************/
/***************************    Define an MCM   *******************************/
/******************************************************************************/
map<uint32_t, uint32_t> Create_MCM(list<uint32_t> MCM_table, bool print_bool)
{
  map<uint32_t, uint32_t> MCM_partition;
  uint32_t integer = 0;
  
  // for (int i=0; i<k; i++)
  for (auto const& i : MCM_table) {
    // MCM_partition[integer]=MCM_table[i];
    MCM_partition[integer]=i;
    integer++;
  }
  
  return MCM_partition;
}

//' Define an MCM
//' 
//' Define an MCM. Input is an integer representation of all MCM parts. For example, 
//' MCM_Choice = c(384, 64, 32, 16, 8, 4, 2, 1) defines an MCM with 8 independent parts, 
//' based on n=9 spins. You can check that the model (i.e., list of parts) that you have 
//' provided properly defines an MCM by calling the function check_partition(Partition, n).
//' 
//' @param MCM_table List containing the integer representation of the MCM parts. 
//' 
//' @return MCM_partition The defined MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<uint32_t>> Create_MCM(std::list<uint32_t> MCM_table) {
  // Create Key, Value vectors for the new MCM_partition!
  vector<uint32_t> key;
  vector<uint32_t> value;
  
  std::vector<std::vector<uint32_t>> MCM_partition;
  
  uint32_t integer = 0;
  
  // for (int i=0; i<k; i++)
  for (auto const& i : MCM_table) {
    // MCM_partition[integer]=MCM_table[i];
    key.push_back(integer);
    value.push_back(i);
    integer++;
  }
  
  MCM_partition.push_back(key);
  MCM_partition.push_back(value);
  
  return MCM_partition;
}

/******************************************************************************/
/*** VERSION a) Operators are written as the binary          ******************/
/****           representation of the interactions           ******************/
/******************************************************************************/
map<uint32_t, uint32_t> Read_MCMParts_BinaryRepresentation(string MCM_binary_filename, unsigned int n, bool print_bool)
{
  map<uint32_t, uint32_t> MCM_partition;
  uint32_t integer = 0;
  
  ifstream myfile (MCM_binary_filename.c_str());
  string line, line2;     
  
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line
      MCM_partition[integer]=stoi(line2, 0, 2);   //convert string line2 into a binary integer
      integer++;
    }
    myfile.close();
  }
  return MCM_partition;
}

//' Read_MCMParts_BinaryRepresentation
//' 
//' Define MCM from a file in which the operators are written as the binary 
//' representation of the interactions.
//' 
//' @param MCM_binary_filename The file name (incl. the path to it from the current directory).
//' @param n The amount of binary (spin) variables.
//' 
//' @return MCM_partition The defined MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<uint32_t>> Read_MCMParts_BinaryRepresentation(std::string MCM_binary_filename, unsigned int n)
{
  // Create Key, Value vectors for the new MCM_partition!
  vector<uint32_t> key;
  vector<uint32_t> value;
  
  std::vector<std::vector<uint32_t>> MCM_partition;
  uint32_t integer = 0;
  
  ifstream myfile (MCM_binary_filename.c_str());
  string line, line2;     
  
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      key.push_back(integer);
      
      line2 = line.substr (0,n);          //take the n first characters of line
      // MCM_partition[integer]=stoi(line2, 0, 2);   //convert string line2 into a binary integer
      value.push_back(stoi(line2, 0, 2));   //convert string line2 into a binary integer
      
      integer++;
    }
    
    MCM_partition.push_back(key);
    MCM_partition.push_back(value);
    
    myfile.close();
  }
  
  return MCM_partition;
}

/********************************************************************/
/*************    CHECK if "Partition" IS A PARTITION   *************/
/********************************************************************/
//check if *Partition* is an actual partition of r basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

bool check_partition(map<uint32_t, uint32_t> Partition, unsigned int n)
{
  map<uint32_t, uint32_t>::iterator Part;
  uint32_t sum = 0;
  uint32_t rank = 0; 
  
  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += countSetBits((*Part).second);
    // Rcpp::Rcout << bitset<n>( (*Part).second ) << " \t";
  }
  // Rcpp::Rcout << bitset<n>(sum) << endl;
  
  return (countSetBits(sum) == rank);
}

//' check_partition
//' 
//' Check if *Partition* is an actual partition of r basis elements,
//' i.e., that no basis element appears in more than 1 part of the partition.
//' i.e., that each basis element only appears in a single part of the partition.
//' 
//' @param Partition Partition of MCM.
//' @param n The amount of binary (spin) variables.
//' 
//' @return Bool_results Returns FALSE if there is an overlap between the parts.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
bool check_partition(std::vector<std::vector<uint32_t>> Partition, unsigned int n)
{
  vector<uint32_t> PartitionValue = Partition.back();

  uint32_t sum = 0;
  uint32_t rank = 0; 
  
  for ( const auto &Part : PartitionValue )
  {
    sum |= Part;
    rank += countSetBits(Part);
  }
  
  return (countSetBits(sum) == rank);
}

/********************************************************************/
/*******    PRINT INFO on each PART of an MCM (= a partition)   *****/
/********************************************************************/
void PrintTerminal_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, unsigned int n, map<uint32_t, uint32_t> MCM_Partition)
{
  uint32_t Part = 0, m=0;
  double C_param=0, C_geom=0;
  Complexity_MCM(MCM_Partition, N, n, &C_param, &C_geom);
  double LogL = LogL_MCM(Kset, MCM_Partition, N, n);
  
   Rcpp::Rcout << "********** General Information about the MCM: **********" << endl; 
   Rcpp::Rcout << "Best MCM has " << MCM_Partition.size() << " partitions and the following properties:" << endl;
   Rcpp::Rcout << "\t LogL = " << LogL << endl;
   Rcpp::Rcout << " \t C_param = " << C_param << " \t \t C_geom = " << C_geom << endl;
   Rcpp::Rcout << " \t Total complexity = " << C_param + C_geom << endl;
   Rcpp::Rcout << " \t MDL = " << LogL - C_param - C_geom << endl;
   Rcpp::Rcout << "  \t LogE = " << LogE_MCM(Kset, MCM_Partition, N, n) << endl;
  
   Rcpp::Rcout << endl << "********** Information about each part of the MCM: **********";
   Rcpp::Rcout << endl << "\t (the total LogE of the model is the sum of the values for each part)";
   Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
   Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
   Rcpp::Rcout << "## 1:Part_int \t 2:Part_binary \t 3:LogL \t 4:C_param \t 5:C_geom \t 6:C_tot \t 7:LogE" << endl;
  
  for (map<uint32_t, uint32_t>::iterator i = MCM_Partition.begin(); i != MCM_Partition.end(); i++)
  {    
    Part = (*i).second;
    m = countSetBits(Part);  // rank of the part (i.e. rank of the SCM)
    C_param = ParamComplexity_SubCM(m, N);
    C_geom = GeomComplexity_SubCM(m);
    
     Rcpp::Rcout << " \t " << Part << " \t " << int_to_bstring(Part, n) << " \t";
     Rcpp::Rcout << LogL_SubCM(Kset, Part, N, n) << " \t";
     Rcpp::Rcout << C_param << " \t " << C_geom << " \t" << C_param + C_geom << " \t";
     Rcpp::Rcout << LogE_SubCM(Kset, Part, N, n) << endl;
  }
   Rcpp::Rcout << endl;
}

//' PRINT INFO on each PART of an MCM (= a partition)
//' 
//' The function prints the number of Independent Components of the MCM 
//' (i.e., the number of partitions -- or communities -- of the MCM),
//' as well as the log-likelihood (LogL), the parameter complexity (C_param), 
//' the geometric complexity (C_geom), the total complexity (C_tot), 
//' the Minimum Description Length (MDL) and the log-evidence (LogE). 
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' @param n The amount of binary (spin) variables.
//' @param MCM_Partition Partition of MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
void PrintTerminal_MCM_Info(std::vector<std::vector<int>> Kset, unsigned int N, unsigned int n, std::vector<std::vector<uint32_t>> MCM_Partition)
{
  uint32_t m=0;
  // double C_param=0, C_geom=0;
  Rcpp::NumericVector C_param, C_geom;
  C_param[0] = 0;
  C_geom[0] = 0;
  
  Complexity_MCM(MCM_Partition, N, n, C_param, C_geom);
  double LogL = LogL_MCM(Kset, MCM_Partition, N, n);

   Rcpp::Rcout << "********** General Information about the MCM: **********" << endl; 
   Rcpp::Rcout << "Best MCM has " << MCM_Partition[1].size() << " partitions and the following properties:" << endl;
   Rcpp::Rcout << "\t LogL = " << LogL << endl;
   Rcpp::Rcout << " \t C_param = " << C_param[0] << " \t \t C_geom = " << C_geom[0] << endl;
   Rcpp::Rcout << " \t Total complexity = " << C_param[0] + C_geom[0] << endl;
   Rcpp::Rcout << " \t MDL = " << LogL - C_param[0] - C_geom[0] << endl;
   Rcpp::Rcout << "  \t LogE = " << LogE_MCM(Kset, MCM_Partition, N, n) << endl;
  
   Rcpp::Rcout << endl << "********** Information about each part of the MCM: **********";
   Rcpp::Rcout << endl << "\t (the total LogE of the model is the sum of the values for each part)";
   Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
   Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
   Rcpp::Rcout << "## 1:Part_int \t 2:Part_binary \t 3:LogL \t 4:C_param \t 5:C_geom \t 6:C_tot \t 7:LogE" << endl;
  
  
  vector<uint32_t> MCM_PartitionValue = MCM_Partition.back();
  for ( const auto &Part : MCM_PartitionValue )
  {    
    // Part = (*i).second;
    m = countSetBits(Part);  // rank of the part (i.e. rank of the SCM)
    C_param[0] = ParamComplexity_SubCM(m, N);
    C_geom[0] = GeomComplexity_SubCM(m);
    
     Rcpp::Rcout << " \t " << Part << " \t " << int_to_bstring(Part, n) << " \t";
     Rcpp::Rcout << LogL_SubCM(Kset, Part, N, n) << " \t";
     Rcpp::Rcout << C_param[0] << " \t " << C_geom[0] << " \t" << C_param[0] + C_geom[0] << " \t";
     Rcpp::Rcout << LogE_SubCM(Kset, Part, N, n) << endl;
  }
   Rcpp::Rcout << endl;
}

/********************************************************************/
/**************************    PRINT INFO    ************************/
/******    ON SUCCESSIVE INDEPENDENT MODELS IN THE NEW BASIS   ******/
/********************************************************************/
void PrintInfo_All_Indep_Models(map<uint32_t, unsigned int> Kset, unsigned int N, unsigned int n)
{
  map<uint32_t, uint32_t> Partition_Indep;  uint32_t Op = 1;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_Indep[i] = Op;
     Rcpp::Rcout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_Indep.clear();
}

//' PRINT INFO ON SUCCESSIVE INDEPENDENT MODELS IN THE NEW BASIS.
//' 
//' The function prints the integer and binary representation of the part 
//' , as well as the the log-likelihood (LogL) and the log-evidence (LogE).
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' @param n The amount of binary (spin) variables.
//'  
//' @export
// [[Rcpp::export(rng=false)]]
void PrintInfo_All_Indep_Models(std::vector<std::vector<int>> Kset, unsigned int N, unsigned int n)
{
  vector<vector<uint32_t>> Partition_Indep;  uint32_t Op = 1;
  
  // Create Key, Value vectors.
  vector<uint32_t> key;
  vector<uint32_t> value;
  
  Partition_Indep.push_back(key);
  Partition_Indep.push_back(value);
  
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_Indep[0].push_back(i);
    Partition_Indep[1].push_back(Op);
    // Partition_Indep[i] = Op;
    
    Rcpp::Rcout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_Indep.clear();
}

/********************************************************************/
/**************************    PRINT INFO    ************************/
/******    ON SUCCESSIVE SUB_COMPLETE MODELS IN THE NEW BASIS   *****/
/********************************************************************/
void PrintInfo_All_SubComplete_Models(map<uint32_t, unsigned int> Kset, unsigned int N, unsigned int n)
{
  map<uint32_t, uint32_t> Partition_SC;  uint32_t Op = 1;
  Partition_SC[0] = 0;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_SC[0] += Op;
    Rcpp::Rcout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_SC, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_SC, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_SC.clear();
}

//' PrintInfo_All_Indep_Models
//' 
//' PRINT INFO ON SUCCESSIVE SUB_COMPLETE MODELS IN THE NEW BASIS.
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' @param n The amount of binary (spin) variables.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
void PrintInfo_All_SubComplete_Models(std::vector<std::vector<int>> Kset, unsigned int N, unsigned int n)
{
  vector<vector<uint32_t>> Partition_SC;  uint32_t Op = 1;
  // Create Key, Value vectors.
  vector<uint32_t> key = {0};
  vector<uint32_t> value = {0};
  
  Partition_SC.push_back(key);
  Partition_SC.push_back(value);
  
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_SC[1][0] += Op;
    Rcpp::Rcout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_SC, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_SC, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_SC.clear();
}
