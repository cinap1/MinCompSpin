//g++ -std=c++11 -O3 main.cpp Data_Manipulation.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp MCM_info.cpp
//time ./a.out

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <cmath>       /* tgamma */

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "library.h"

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/

//' Perform the exhaustive search for the best MCM at once.
//' 
//' Performs the community detection on the data given in datafilename with the given basis Basis_li 
//' and returns the best fitting MCM for the data. Optional: set the comparison
//' variable to TRUE to test what happens when the model is compared to 3 different cases.
//' For general use set to FAlSE. Output is stored in the (newly created)
//' OUTPUT directory in the same working directory.
//' 
//' @param n Number of binary (spin) variables of the data set.
//' @param datafilename The name of input file including path containing the binary data.
//' @param Basis_li The basis on which the model is based. 
//' @param OUTPUT_directory The name of the output directory (and path to it) in which to save the output data. 
//' @param comparison Boolean indicating whether to run 2 other versions of the exhaustive search. For comparison purposes. Default: FALSE. 
//' 
//' @export
// [[Rcpp::export(rng=false)]]
int MCM_search(unsigned int n, std::string datafilename, std::list<uint32_t> Basis_li, std::string OUTPUT_directory, bool comparison)
{
   Rcpp::Rcout << "--->> Create OUTPUT Folder: (if needed) ";
   int status = system( ("mkdir -p " + OUTPUT_directory + "/").c_str() );
   // if (status < 0) Rcpp::Rcout << "Error: " << strerror(errno) << '\n';
  // int ret = system("mkdir -p OUTPUT/");
   Rcpp::Rcout << endl;

   Rcpp::Rcout << endl << "*******************************************************************************************";
   Rcpp::Rcout << endl << "***********************************  Read the data:  **************************************";
   Rcpp::Rcout << endl << "*******************************************************************************************" << endl;
  unsigned int N=0; // will contain the number of datapoints in the dataset
  map<uint32_t, unsigned int> Nset = read_datafile(&N, n, datafilename);
  
  if (N == 0) {
    return 1;
  }

   Rcpp::Rcout << endl << "*******************************************************************************************";
   Rcpp::Rcout << endl << "******************************  Choice of the basis:  *************************************";
   Rcpp::Rcout << endl << "*******************************************************************************************" << endl;

   Rcpp::Rcout << endl << "Choice of the basis for building the Minimally Complex Model (MCM):" << endl;

  // *** Basis elements are written using the integer representation of the operator
  // *** For instance, a basis element on the last two spin variable would be written:
  // ***      -->  Op = s1 s2           Spin operator
  // ***      -->  Op = 000000011       Binary representation
  // ***      -->  Op = 3               Integer representation   ( 000000011 = 3 )
  
  // *** The basis can be specified by hand here:
  // uint32_t Basis_Choice[] =  {3, 5, 9, 48, 65, 129, 272, 81, 1};    // Ex. This is the best basis for the "Shapes" dataset
  
  // unsigned int m = sizeof(Basis_Choice) / sizeof(uint32_t);
  // list<uint32_t> Basis_li;  Basis_li.assign (Basis_Choice, Basis_Choice + m);
  
  // *** The basis can also be read from a file:
  //   list<uint32_t> Basis_li = Read_BasisOp_IntegerRepresentation();
  //   list<uint32_t> Basis_li = Read_BasisOp_BinaryRepresentation();
  
  // *** Or one can simply use the original basis of the data:
  //   list<uint32_t> Basis_li = Original_Basis();

  // *** Print info about the Basis:
  PrintTerm_Basis(Basis_li, n);

  Rcpp::Rcout << "Number of spin variables, n=" << n << endl;
  Rcpp::Rcout << "Number of basis elements, m=" << Basis_li.size() << endl;
  if (Basis_li.size() > n) {  Rcpp::Rcout << " -->  Error: the number 'm' of basis elements is larger than the size 'n' of the system." << endl;  }
  else {
     Rcpp::Rcout << " -->  m <= n :  Everything seems fine." << endl;
     Rcpp::Rcout << "Make sure that the set of basis elements provided are orthogonal to each other." << endl;
  }

  Rcpp::Rcout << endl << "*******************************************************************************************";
  Rcpp::Rcout << endl << "*********************************  Change the data basis   ********************************";
  Rcpp::Rcout << endl << "**************************************  Build Kset:  **************************************";
  Rcpp::Rcout << endl << "*******************************************************************************************" << endl;
  
  Rcpp::Rcout << endl << "Transform the data in the specified basis." << endl;
  Rcpp::Rcout << endl << "/!\\ If m<n:";
  Rcpp::Rcout << endl << "\t If the size 'm' of the basis is strictly smaller than the number 'n' of variables, ";
  Rcpp::Rcout << endl << "\t then the data will be troncated to the 'm' first basis elements." << endl;

  map<uint32_t, unsigned int> Kset = build_Kset(Nset, Basis_li, n, false);

  Rcpp::Rcout << endl << "*******************************************************************************************";
  Rcpp::Rcout << endl << "************************************  All Indep Models:  **********************************";
  Rcpp::Rcout << endl << "*******************************************************************************************" << endl << endl;

  Rcpp::Rcout << "Independent models in the new basis:" << endl;
  
  PrintInfo_All_Indep_Models(Kset, N, n);

  Rcpp::Rcout << endl << "*******************************************************************************************";
  Rcpp::Rcout << endl << "**************************  All Successive Sub-Complete Models:  **************************";
  Rcpp::Rcout << endl << "*******************************************************************************************" << endl << endl;
  
  Rcpp::Rcout << "Sub-Complete models in the new basis:" << endl;

  PrintInfo_All_SubComplete_Models(Kset, N, n);

  Rcpp::Rcout << endl << "*******************************************************************************************";
  Rcpp::Rcout << endl << "***********************************  VERSION 1  *******************************************";
  Rcpp::Rcout << endl << "*******************************************************************************************";
  Rcpp::Rcout << endl << "**********************  Compare all MCMs of a given rank 'r' ******************************";
  Rcpp::Rcout << endl << "*********************  based on the 'r' first basis Operators:  ***************************";
  Rcpp::Rcout << endl << "*******************************************************************************************" << endl << endl;
  
  Rcpp::Rcout << endl << "Search among all MCMs based on the 'r' first basis operators (i.e., the models of rank exactly equal to 'r')" << endl;
  Rcpp::Rcout << endl << "/!\\ Conditions on the value of 'r':  r <= m <= n ";
  Rcpp::Rcout << endl << "\t 'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
  Rcpp::Rcout << endl << "\t which must be smaller or equal to the number 'n' of spin variables." << endl << endl;

  int r1 = n;
  double LogE_BestMCM1 = 0;

  if (r1 <= Basis_li.size())
  {
    map<uint32_t, uint32_t> MCM_Partition1 = MCM_GivenRank_r(Kset, N, &LogE_BestMCM1, r1, n, OUTPUT_directory, false);
    // Rcpp::Rcout << "\t Best LogE = " << LogE_BestMCM1 << endl;
    PrintTerminal_MCM_Info(Kset, N, n, MCM_Partition1);
  }
  else {  Rcpp::Rcout << "The condition on the value of 'r' is not respected" << endl;  }
  
  if (comparison == true) 
  {
     Rcpp::Rcout << endl << "*******************************************************************************************";
     Rcpp::Rcout << endl << "***********************************  VERSION 2  *******************************************";
     Rcpp::Rcout << endl << "*******************************************************************************************";
     Rcpp::Rcout << endl << "****************  Compare all MCMs of rank 'k', with 1 <= k <=r' **************************";
     Rcpp::Rcout << endl << "******************  based on the 'k' first basis Operators:  ******************************";
     Rcpp::Rcout << endl << "*******************************************************************************************" << endl << endl;
  
     Rcpp::Rcout << endl << "Search among all MCMs based on the 'k' FIRST basis operators provided, for all k=1 to r" << endl;
     Rcpp::Rcout << endl << "/!\\ Conditions on the value of 'r':  r <= m <= n ";
     Rcpp::Rcout << endl << "\t 'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
     Rcpp::Rcout << endl << "\t which must be smaller or equal to the number 'n' of spin variables." << endl << endl;
  
     Rcpp::Rcout << endl << "\t Check function declaration for the default options." << endl << endl;
  
  /******************************************************************************/
  // *** By default: - r=n
  // ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true
  /******************************************************************************/
  
    int r2 = n;
    double LogE_BestMCM2 = 0;
  
    if (r2 <= Basis_li.size())
    {
      map<uint32_t, uint32_t> MCM_Partition2 = MCM_AllRank_SmallerThan_r_Ordered(Kset, N, &LogE_BestMCM2, r2, n, OUTPUT_directory, false);
      // Rcpp::Rcout << "\t Best LogE = " << LogE_BestMCM2 << endl;
      PrintTerminal_MCM_Info(Kset, N, n, MCM_Partition2);
    }
    else {  Rcpp::Rcout << "The condition on the value of 'r' is not respected" << endl;  }
  
  
     Rcpp::Rcout << endl << "*******************************************************************************************";
     Rcpp::Rcout << endl << "***********************************  VERSION 3  *******************************************";
     Rcpp::Rcout << endl << "*******************************************************************************************";
     Rcpp::Rcout << endl << "****************  Compare all MCMs of rank 'k', with 1 <= k <=r' **************************";
     Rcpp::Rcout << endl << "***********  based on any 'k' subset of the basis Operators provided  *********************";
     Rcpp::Rcout << endl << "*******************************************************************************************" << endl << endl;
  
     Rcpp::Rcout << endl << "Search among all MCMs based on ANY SUBSET of 'k' operators of the basis provided, for all k=1 to r" << endl;
     Rcpp::Rcout << endl << "/!\\ Conditions on the value of 'r':  r <= m <= n ";
     Rcpp::Rcout << endl << "\t 'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
     Rcpp::Rcout << endl << "\t which must be smaller or equal to the number 'n' of spin variables." << endl << endl;
  
     Rcpp::Rcout << endl << "\t Check function declaration for the default options." << endl << endl;
  
    /******************************************************************************/
    // *** By default: - r=n
    // ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true
    /******************************************************************************/
  
    int r3 = n;
    double LogE_BestMCM3 = 0;
  
    if (r3 <= Basis_li.size())
    {
      map<uint32_t, uint32_t> MCM_Partition3 = MCM_AllRank_SmallerThan_r_nonOrdered(Kset, N, &LogE_BestMCM3, r3, n, OUTPUT_directory, false);
      // Rcpp::Rcout << "\t Best LogE = " << LogE_BestMCM3 << endl;
      PrintTerminal_MCM_Info(Kset, N, n, MCM_Partition3);
    }
    else {  Rcpp::Rcout << "The condition on the value of 'r' is not respected" << endl;  }
    
  }

  return 0;
}
