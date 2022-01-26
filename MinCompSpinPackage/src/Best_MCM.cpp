#include <Rcpp.h>
using namespace Rcpp;

#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "support.h"

using namespace std;

/******************************************************************************/
/**************** Log-likelihood (LogL), Geometric Complexity *****************/
/*************************  and Log-evidence (LogE) ***************************/
/******************************************************************************/
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false);
double LogE_MCM(std::vector<std::vector<int>> Kset, std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, bool print_bool = false);

double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, double *C_param, double *C_geom);
double Complexity_MCM(std::vector<std::vector<uint32_t>> Partition, unsigned int N, unsigned int n, Rcpp::NumericVector C_param, Rcpp::NumericVector C_geom);

/********************************************************************/
/*************    CHECK if "Partition" IS A PARTITION   *************/
/********************************************************************/
//check if *Partition* is an actual partition of r basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

bool check_partition(map<uint32_t, uint32_t> Partition);
bool check_partition(std::vector<std::vector<uint32_t>> Partition, unsigned int n);

/********************************************************************/
/**********************    PRINT PARTITION   ************************/
/********************************************************************/
void Print_Partition(uint32_t *a, unsigned int n)
{
  for (int i=0; i<n; i++)
  {     Rcpp::Rcout << a[i];  }
}

void Print_Partition_Converted(map<uint32_t, uint32_t>  partition, unsigned int n)
{
  for (map<uint32_t, uint32_t>::iterator i = partition.begin(); i != partition.end(); i++)
  {     Rcpp::Rcout << (*i).second << " = " << int_to_bstring((*i).second, n) << "\n";  }
   Rcpp::Rcout << endl;
}

/********************************************************************/
/****************    CONVERSION  of a partition    ******************/
/**********************   SPECIFIC TO MCM    ************************/
/********************************************************************/
// *** map<uint32_t, uint32_t>   --> .first = i = index    --> .second = a[i] = number of element in the part
map<uint32_t, uint32_t> Convert_Partition_forMCM(uint32_t *a, unsigned int r)
{
  map<uint32_t, uint32_t> Partition;
  uint32_t element = 1;
  
  for (int i=r-1; i>=0; i--)  // read element from last to first
  {    
    Partition[(a[i])] += element;
    element = element << 1;      // Rcpp::Rcout << a[i] << "\t :";   Print_Partition_Converted(Partition);
  }
  
  //  Rcpp::Rcout << "Convert, " << Partition.size() << " parts: \t " ;
  //  bool ok = check_partition(Partition);
  //  Rcpp::Rcout << " \t --> Partition Ok? " << ok << endl << endl;
  
  return Partition;
}

// The vector based version of the above.
vector<vector<uint32_t>> Convert_Partition_forMCM_vector(uint32_t *a, unsigned int r)
{
  // map<uint32_t, uint32_t> Partition;
  vector<uint32_t> key;
  vector<uint32_t> value;
  
  vector<vector<uint32_t>> Partition;
  Partition.push_back(key);
  Partition.push_back(value);
  
  vector<uint32_t>::iterator itr;
  uint32_t a_key;
  
  uint32_t element = 1;
  
  for (int i=r-1; i>=0; i--)  // read element from last to first
  {    
    a_key = a[i];
    
    itr = std::find(Partition[0].begin(), Partition[0].end(), a_key); 
    
    if (itr != Partition[0].cend()) // Key was present so use index to adjust the value vector at the right index.
    { 
      Partition[1].at(std::distance(Partition[0].begin(), itr)) += element; // Update value
    }
    else // Key was not present so add new Key Value pair to their corresponding vectors.
    { 
      Partition[0].push_back(a_key);
      Partition[1].push_back(element);
    }
    
    element = element << 1;
  }
  
  return Partition;
}

// *** map<uint32_t, uint32_t> --> .first = i = index of the partition    --> .second = a[i] = number of element in the part
map<uint32_t, uint32_t> Convert_Partition_forMCM_withSubPart(uint32_t *a, bool *keep_SubPartition, unsigned int r)
{
  map<uint32_t, uint32_t> Partition;
  
  uint32_t element = 1;
  bool switch_ = false;
  *keep_SubPartition = true;
  
  for (int i=r-1; i>=0; i--)  // read element from last to first
  {    
    Partition[(a[i])] += element;  //  Rcpp::Rcout << a[i] << "\t ";
    element = element << 1;
    if(switch_ == true && a[i] != 0)  { *keep_SubPartition = false;  }
    else if(a[i] == 0) { switch_ = true;  }
  }
  
  //  Rcpp::Rcout << "Convert, " << Partition.size() << " parts: \t " ;
  //  bool ok = check_partition(Partition);
  //  Rcpp::Rcout << " \t --> Partition Ok? " << ok << endl << endl;
  
  return Partition;
}

// Vector based version of the above.
vector<vector<uint32_t>> Convert_Partition_forMCM_withSubPart_vector(uint32_t *a, bool *keep_SubPartition, unsigned int r)
{
  vector<uint32_t> key;
  vector<uint32_t> value;
  
  vector<vector<uint32_t>> Partition;
  Partition.push_back(key);
  Partition.push_back(value);
  
  vector<uint32_t>::iterator itr;
  uint32_t a_key;
  
  uint32_t element = 1;
  bool switch_ = false;
  *keep_SubPartition = true;
  
  for (int i=r-1; i>=0; i--)  // read element from last to first
  { 
    // Partition[(a[i])] += element;  //  Rcpp::Rcout << a[i] << "\t ";
    a_key = a[i];

    itr = std::find(Partition[0].begin(), Partition[0].end(), a_key);
    
    if (itr != Partition[0].cend()) // Key was present.
    { 
      Partition[1].at(std::distance(Partition[0].begin(), itr)) += element; // Update value
    }
    else // Add new Key Value pair to their corresponding vectors.
    { 
      Partition[0].push_back(a_key);
      Partition[1].push_back(element);
    }
    
    element = element << 1;
    if(switch_ == true && a[i] != 0)  { *keep_SubPartition = false;  }
    else if(a[i] == 0) { switch_ = true;  }
  }
  return Partition;
}

/******************************************************************************/
/*********************  Compute all Partitions of a set   *********************/
/***************************   with Algorithm H   *****************************/
/******************************************************************************/
// *** find the first index j (from the right) such that a[j] != b[j]
int find_j(uint32_t *a, uint32_t *b, unsigned int r)
{
  int j = r-2;
  while (a[j] == b[j])  {   j--;  }
  return j;
}

/******************************************************************************/
// *** Version 1: 
// ***            Compare all the MCM of rank r, 
// ***            based on the r first elements of the basis used to build Kset:
/******************************************************************************/
map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r, unsigned int n, string OUTPUT_directory, bool print_bool)
{
   Rcpp::Rcout << "--->> Search for the best MCM.." << endl << endl;
  
  int counter = 0, i = 0;
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM((OUTPUT_directory + "/BestMCM_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file all MCMs:
  fstream file_MCM_Rank_r((OUTPUT_directory + "/AllMCMs_Rank_r" + to_string(r) + ".dat").c_str(), ios::out);
  if(print_bool)
  {
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    file_MCM_Rank_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_MCM_Rank_r << "To activate the prints for all the MCMs of rank r="<< r << ","<< endl;
    file_MCM_Rank_r << " specify `print_bool=true` in the last argument of the function MCM_GivenRank_r();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition;
  
  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(n*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  
  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N, n);
  
  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_MCM_Rank_r << counter << ": \t";
    Partition = Convert_Partition_forMCM(a, r);
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    
    // *** Print in file:
    if(print_bool)
    {
      file_MCM_Rank_r << xx_st;
      for (i=0; i<r; i++)   {    file_MCM_Rank_r << a[i];  }     //Print_Partition(a);
      
      Complexity_MCM(Partition, N, n, &C_param, &C_geom);    //Complexity
      file_MCM_Rank_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (*LogE_best) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[r-1] up to reaching b[r-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[r-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_MCM_Rank_r.close();
  
  Rcpp::Rcout << "--> Number of MCModels (of rank r=" << r << ") that were compared: " << counter << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<r; i++) {   Rcpp::Rcout << aBest[i];  }
  Rcpp::Rcout << "\t \t LogE = " << (*LogE_best) << endl << endl;
  
  Partition = Convert_Partition_forMCM(aBest, r);
  free(a); free(b); free(aBest);
  
  return Partition;
}

// VECTOR BASED VER. 1
//' MCM_GivenRank_r
//' 
//' Compare all the MCM of rank r,
//' based on the r first elements of the basis used to build Kset:
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' @param LogE_best The double which will contain the maximized Log-evidence.
//' @param r The rank r which is the r first elements of the basis.
//' @param n The amount of binary (spin) variables.
//' @param OUTPUT_directory The name of the output directory (and path to it) in which to save the output data. 
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return Partition The best MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<uint32_t>> MCM_GivenRank_r(std::vector<std::vector<int>> Kset, unsigned int N, Rcpp::NumericVector LogE_best, unsigned int r, unsigned int n, std::string OUTPUT_directory, bool print_bool)
{
   Rcpp::Rcout << "--->> Search for the best MCM.." << endl << endl;
  
  int counter = 0, i = 0;
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM((OUTPUT_directory + "/BestMCM_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file all MCMs:
  fstream file_MCM_Rank_r((OUTPUT_directory + "/AllMCMs_Rank_r" + to_string(r) + ".dat").c_str(), ios::out);
  if(print_bool)
  {
    Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
    Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    file_MCM_Rank_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_MCM_Rank_r << "To activate the prints for all the MCMs of rank r="<< r << ","<< endl;
    file_MCM_Rank_r << " specify `print_bool=true` in the last argument of the function MCM_GivenRank_r();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  // double C_param = 0, C_geom = 0;
  Rcpp::NumericVector C_param, C_geom;
  C_param[0] = 0;
  C_geom[0] = 0;
  
  // map<uint32_t, uint32_t> Partition; 
  vector<vector<uint32_t>> Partition;
  
  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(n*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  
  LogE_best[0] = LogE_MCM(Kset, Convert_Partition_forMCM_vector(a, r), N, n);
  
  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_MCM_Rank_r << counter << ": \t";
    Partition = Convert_Partition_forMCM_vector(a, r);
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    
    // *** Print in file:
    if(print_bool)
    {
      file_MCM_Rank_r << xx_st;
      for (i=0; i<r; i++)   {    file_MCM_Rank_r << a[i];  }     //Print_Partition(a);
      
      Complexity_MCM(Partition, N, n, C_param, C_geom);    //Complexity
      file_MCM_Rank_r << " \t" << LogE << " \t" << C_param[0] << " \t" << C_geom[0] << " \t" << (C_param[0] + C_geom[0]) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (LogE_best[0])) 
    { 
      LogE_best[0] = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (LogE_best[0]) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[r-1] up to reaching b[r-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[r-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_MCM_Rank_r.close();
  
  Rcpp::Rcout << "--> Number of MCModels (of rank r=" << r << ") that were compared: " << counter << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<r; i++) {   Rcpp::Rcout << aBest[i];  }
  Rcpp::Rcout << "\t \t LogE = " << (LogE_best[0]) << endl << endl;
  
  Partition = Convert_Partition_forMCM_vector(aBest, r);
  free(a); free(b); free(aBest);
  
  return Partition;
}

/******************************************************************************/
// *** Version 2:  
// ***            Compare all the MCM 
// ***            based on the k first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size()
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_Ordered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r, unsigned int n, string OUTPUT_directory, bool print_bool)
{
  int counter = 0, i = 0;
  int counter_subMCM = 0;
  
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "/BestMCM_Rank_r<=" + to_string(r) + "_Ordered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file all models:
  fstream file_allMCM_r((OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat").c_str(), ios::out);
  
  if(print_bool)
  {
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;
    
    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition;
  
  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  
  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N, n);
  
  
  // *** SubPartitions (rank < n):
  bool keep_SubPartition = false;
  
  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_allMCM_r << counter << ": \t";
    
    // *** Original Partition:
    Partition = Convert_Partition_forMCM_withSubPart(a, &keep_SubPartition, r);     //Print_Partition_Converted(Partition); 
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    
    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }     //Print_Partition(a);
      
      Complexity_MCM(Partition, N, n, &C_param, &C_geom);    //Complexity
      file_allMCM_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (*LogE_best) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }
    
    // *** Sub-Partition:
    if (keep_SubPartition)
    {
      counter_subMCM++;
      
      Partition.erase(0); //Print_Partition_Converted(Partition); 
      LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
      
      // *** Print in file:
      if(print_bool)
      { 
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++)     //Print_Partition(a);
          {
          if (a[i] == 0 )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << (a[i]-1);  } 
          }
        
        Complexity_MCM(Partition, N, n, &C_param, &C_geom);    //Complexity
        file_allSubMCM << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter_subMCM << endl;
      }
      
      // *** Best MCM LogE:
      if ( LogE > (*LogE_best) )
      { 
        *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a); 
        file_BestMCM << xx_st; 
        for (i=0; i<r; i++)   
        {    
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (*LogE_best) )
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);   aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";  aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();
  
  Rcpp::Rcout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<r; i++) {  if(aBest[i] != -1)  { Rcpp::Rcout << aBest[i];}   else { Rcpp::Rcout << "x";}   }
  Rcpp::Rcout << "\t \t LogE = " << (*LogE_best) << endl << endl;
  
  Partition = Convert_Partition_forMCM(aBest, r); 
  free(a); free(b); free(aBest);
  
  return Partition;
}

// VECTOR BASED VER. 2
//' MCM_AllRank_SmallerThan_r_Ordered
//' 
//' Compare all the MCM 
//' based on the k first elements of the basis used to build Kset
//' for all k=1 to r, where r <= basis.size().
//' By default: - r=n
//'             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the data-set.
//' @param LogE_best The double which will contain the maximized Log-evidence.
//' @param r The rank r which is the r first elements of the basis.
//' @param n The amount of binary (spin) variables.
//' @param OUTPUT_directory The name of the output directory (and path to it) in which to save the output data. 
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return Partition The best MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<uint32_t>> MCM_AllRank_SmallerThan_r_Ordered(std::vector<std::vector<int>> Kset, unsigned int N, Rcpp::NumericVector LogE_best, unsigned int r, unsigned int n, std::string OUTPUT_directory, bool print_bool)
{
  int counter = 0, i = 0;
  int counter_subMCM = 0;
  
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "/BestMCM_Rank_r<=" + to_string(r) + "_Ordered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file all models:
  fstream file_allMCM_r((OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat").c_str(), ios::out);
  
  if(print_bool)
  {
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;
    
    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  
  // double C_param = 0, C_geom = 0;
  Rcpp::NumericVector C_param, C_geom;
  C_param[0] = 0;
  C_geom[0] = 0;
  
  // map<uint32_t, uint32_t> Partition;
  vector<vector<uint32_t>> Partition;
  
  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  
  LogE_best[0] = LogE_MCM(Kset, Convert_Partition_forMCM_vector(a, r), N, n);
  
  
  // *** SubPartitions (rank < n):
  bool keep_SubPartition = false;
  
  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_allMCM_r << counter << ": \t";
    
    // *** Original Partition:
    Partition = Convert_Partition_forMCM_withSubPart_vector(a, &keep_SubPartition, r);     //Print_Partition_Converted(Partition); 
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    
    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }     //Print_Partition(a);
      
      Complexity_MCM(Partition, N, n, C_param[0], C_geom[0]);    //Complexity
      file_allMCM_r << " \t" << LogE << " \t" << C_param[0] << " \t" << C_geom[0] << " \t" << (C_param[0] + C_geom[0]) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (LogE_best[0])) 
    { 
      LogE_best[0] = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (LogE_best[0]) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }
    
    // *** Sub-Partition:
    if (keep_SubPartition)
    {
      counter_subMCM++;
      
      // Partition.erase(0); //Print_Partition_Converted(Partition); 
      vector<uint32_t>::iterator itr;
      itr = std::find(Partition[0].begin(), Partition[0].end(), 0);
      
      if (itr != Partition[0].cend()) { // Key was present.
        int idx = std::distance(Partition[0].begin(), itr); // Determine index.
        // Erase key from key vector and value from value vector.
        Partition[0].erase( Partition[0].begin() + idx);  
        Partition[1].erase( Partition[1].begin() + idx); //Print_Partition_Converted(Partition); 
      }
      
      
      LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
      
      // *** Print in file:
      if(print_bool)
      { 
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++)     //Print_Partition(a);
          {
          if (a[i] == 0 )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << (a[i]-1);  } 
          }
        
        Complexity_MCM(Partition, N, n, C_param, C_geom);    //Complexity
        file_allSubMCM << " \t" << LogE << " \t" << C_param[0] << " \t" << C_geom[0] << " \t" << (C_param[0] + C_geom[0]) << " \t" << counter_subMCM << endl;
      }
      
      // *** Best MCM LogE:
      if ( LogE > (LogE_best[0]) )
      { 
        LogE_best[0] = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a); 
        file_BestMCM << xx_st; 
        for (i=0; i<r; i++)   
        {    
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (LogE_best[0]) )
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);   aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";  aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();
  
  Rcpp::Rcout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<r; i++) {  if(aBest[i] != -1)  { Rcpp::Rcout << aBest[i];}   else { Rcpp::Rcout << "x";}   }
  Rcpp::Rcout << "\t \t LogE = " << (LogE_best[0]) << endl << endl;
  
  Partition = Convert_Partition_forMCM_vector(aBest, r); 
  free(a); free(b); free(aBest);
  
  return Partition;
}

/******************************************************************************/
// *** Version 3:  
// ***            Compare all the MCMs based on any subset of k elements 
// ***            of the of r first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size() 
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_nonOrdered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r, unsigned int n, string OUTPUT_directory, bool print_bool)
{
   Rcpp::Rcout << "All MCM based on all subsets of r operators among n chosen independent operators, r<=n: " << endl;
  
  int counter = 0, i = 0;
  int counter_subMCM = 0;
  
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "/BestMCM_Rank_r<=" + to_string(r) + "_NonOrdered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file:
  fstream file_allMCM_r((OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_NonOrdered.dat").c_str(), ios::out);
  
  if(print_bool)
  {
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;
    
    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition, Partition_buffer;
  
  //  *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N, n);
  
  // *** for SubModels:
  uint32_t amax = 0, atest = 0;
  
  //SubPartitions (rank < n):
  
  //ALGO H:
  while(j != 0) // && counter < 200)
  {
    // *** H2: Visit:   ******
    counter++;
    
    // *** Partition:
    Partition = Convert_Partition_forMCM(a, r); 
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    Complexity_MCM(Partition, N, n, &C_param, &C_geom);    //Complexity
    
    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }
      file_allMCM_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE; 
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New" << endl;  
    }
    else if ( LogE == (*LogE_best)) 
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
    }
    
    // *** Find max value in a[]:  //amax=0;   for(i=0; i<n; i++)  {  if (a[i] > amax) { amax = a[i]; } }
    if ( a[r-1] == b[r-1] ) { amax = b[r-1]; } else { amax = b[r-1]-1; } 
    
    // *** Sub-Partition: ***************************** //
    for(atest=0; atest<=amax; atest++)
    {
      counter_subMCM++;
      
      // *** Partition:
      Partition_buffer = Partition;
      Partition_buffer.erase(atest);
      
      LogE = LogE_MCM(Kset, Partition_buffer, N, n);     //LogE
      Complexity_MCM(Partition_buffer, N, n, &C_param, &C_geom);    //Complexity
      
      // *** Print in file:
      if(print_bool)
      {
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++) 
        {
          if (a[i] == atest )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << a[i];  } 
        }
        file_allSubMCM << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter_subMCM << endl;
      }
      
      // *** Best MCM LogE:
      if ( LogE > (*LogE_best)) 
      { 
        *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
        
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)   
        {    
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (*LogE_best)) 
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();
  
  Rcpp::Rcout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<n; i++) {  if(aBest[i] != -1)  { Rcpp::Rcout << aBest[i];}   else { Rcpp::Rcout << "x";}   }
  Rcpp::Rcout << "\t \t LogE = " << (*LogE_best) << endl << endl;
  
  Partition = Convert_Partition_forMCM(aBest, r); 
  free(a); free(b); free(aBest);
  
  return Partition;
}

// VECTOR BASED VER. 3
//' MCM_AllRank_SmallerThan_r_Ordered
//' 
//' Compare all the MCMs based on any subset of k elements
//' of the of r first elements of the basis used to build Kset
//' for all k=1 to r, where r <= basis.size()
//' By default: - r=n
//'             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true
//' 
//' @param Kset Mapping of frequency of occurrence of each state of the data-set in the new basis.
//' @param N The size of the dataset
//' @param LogE_best The double which will contain the maximized Log-evidence.
//' @param r The rank r which is the r first elements of the basis.
//' @param n The amount of binary (spin) variables.
//' @param OUTPUT_directory The name of the output directory (and path to it) in which to save the output data. 
//' @param print_bool Boolean to print related information to terminal. Default: FALSE.
//' 
//' @return Partition The best MCM.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::vector<std::vector<uint32_t>> MCM_AllRank_SmallerThan_r_nonOrdered(std::vector<std::vector<int>> Kset, unsigned int N, Rcpp::NumericVector LogE_best, unsigned int r, unsigned int n, std::string OUTPUT_directory, bool print_bool)
{
   Rcpp::Rcout << "All MCM based on all subsets of r operators among n chosen independent operators, r<=n: " << endl;
  
  int counter = 0, i = 0;
  int counter_subMCM = 0;
  
  string xx_st = "";
  for(int i=0; i<n-r; i++)
  {  xx_st += "_";  }
  
  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "/BestMCM_Rank_r<=" + to_string(r) + "_NonOrdered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;
  
  // *** Print in file:
  fstream file_allMCM_r((OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_NonOrdered.dat").c_str(), ios::out);
  
  if(print_bool)
  {
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    
     Rcpp::Rcout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
     Rcpp::Rcout << (OUTPUT_directory + "/AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;
    
    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }
  
  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;
  
  // *** LogE and Complexity
  double LogE = 0;
  
  // double C_param = 0, C_geom = 0;
  Rcpp::NumericVector C_param, C_geom;
  C_param[0] = 0;
  C_geom[0] = 0;
  
  // map<uint32_t, uint32_t> Partition, Partition_buffer;
  vector<vector<uint32_t>> Partition, Partition_buffer;
  
  //  *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  LogE_best[0] = LogE_MCM(Kset, Convert_Partition_forMCM_vector(a, r), N, n);
  
  // *** for SubModels:
  uint32_t amax = 0, atest = 0;
  
  //SubPartitions (rank < n):
  
  //ALGO H:
  while(j != 0) // && counter < 200)
  {
    // *** H2: Visit:   ******
    counter++;
    
    // *** Partition:
    Partition = Convert_Partition_forMCM_vector(a, r); 
    LogE = LogE_MCM(Kset, Partition, N, n);     //LogE
    Complexity_MCM(Partition, N, n, C_param, C_geom);    //Complexity
    
    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }
      file_allMCM_r << " \t" << LogE << " \t" << C_param[0] << " \t" << C_geom[0] << " \t" << (C_param[0] + C_geom[0]) << " \t" << counter << endl;
    }
    
    // *** Best MCM LogE:
    if ( LogE > (LogE_best[0])) 
    { 
      LogE_best[0] = LogE; 
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New" << endl;  
    }
    else if ( LogE == (LogE_best[0])) 
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
    }
    
    // *** Find max value in a[]:  //amax=0;   for(i=0; i<n; i++)  {  if (a[i] > amax) { amax = a[i]; } }
    if ( a[r-1] == b[r-1] ) { amax = b[r-1]; } else { amax = b[r-1]-1; } 
    
    // *** Sub-Partition: ***************************** //
    vector<uint32_t>::iterator itr;
    for(atest=0; atest<=amax; atest++)
    {
      counter_subMCM++;
      
      // *** Partition:
      Partition_buffer = Partition;
      
      // Partition_buffer.erase(atest);
      itr = std::find(Partition_buffer[0].begin(), Partition_buffer[0].end(), atest);
      
      if (itr != Partition_buffer[0].cend()) { // Key was present.
        int idx = std::distance(Partition_buffer[0].begin(), itr); // Determine index.
        // Erase key from key vector and value from value vector.
        Partition_buffer[0].erase( Partition_buffer[0].begin() + idx);  
        Partition_buffer[1].erase( Partition_buffer[1].begin() + idx); //Print_Partition_Converted(Partition); 
      }
      
      LogE = LogE_MCM(Kset, Partition_buffer, N, n);     //LogE
      Complexity_MCM(Partition_buffer, N, n, C_param, C_geom);    //Complexity
      
      // *** Print in file:
      if(print_bool)
      {
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++) 
        {
          if (a[i] == atest )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << a[i];  } 
        }
        file_allSubMCM << " \t" << LogE << " \t" << C_param[0] << " \t" << C_geom[0] << " \t" << (C_param[0] + C_geom[0]) << " \t" << counter_subMCM << endl;
      }
      
      // *** Best MCM LogE:
      if ( LogE > (LogE_best[0])) 
      { 
        LogE_best[0] = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
        
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)   
        {    
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (LogE_best[0])) 
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }
    
    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }
  
  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();
  
  Rcpp::Rcout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl;
  
  Rcpp::Rcout << endl << "********** Best MCM: **********";
  Rcpp::Rcout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  Rcpp::Rcout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  
  Rcpp::Rcout << "\t >> Best Model = ";
  Rcpp::Rcout << xx_st;
  for(int i=0; i<n; i++) {  if(aBest[i] != -1)  { Rcpp::Rcout << aBest[i];}   else { Rcpp::Rcout << "x";}   }
  Rcpp::Rcout << "\t \t LogE = " << (LogE_best[0]) << endl << endl;
  
  Partition = Convert_Partition_forMCM_vector(aBest, r); 
  free(a); free(b); free(aBest);
  
  return Partition;
}
