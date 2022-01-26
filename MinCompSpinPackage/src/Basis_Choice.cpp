#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "support.h"

using namespace std;

/******************************************************************************/
/*****************   Read Basis Operators from file  **************************/
/******************************************************************************/

/******  1)  Operators should be written in one single column          ********/
/******  2)  Operators can be written in two different versions:       ********/
/***   a) as a binary representation of the spin involved in the Operator;    */
/***   b) as the integer value of that binary representation.                 */
/******************************************************************************/
/****  Ex. for a system with n=4 spin:  ***************************************/
/****      -- the operator s1 would be encoded as 0001,                      **/
/****      which has the integer representation 1  -->   0001 = 1      ********/
/****      -- the operator s1 s2 would be encoded as 0011,             ********/
/****      which has the integer representation 3  -->   0011 = 3      ********/
/******************************************************************************/


/******************************************************************************/
/*** VERSION a) Operators are written as the binary          ******************/
/****           representation of the interactions           ******************/
/******************************************************************************/
//' Read basis operators in binary representation from input file.
//' 
//' Reads the Basis operators in binary representation from the given input 
//' file. The operators should be written in one single column and as 
//' binary representation of the spin involved in the Operator. Returns a 
//' list of the basis with the integer value of that binary representation
//' 
//' @param Basis_binary_filename The name of input file including path.
//' @param n Number of binary (spin) variables of the data set.
//' 
//' @return Basis_li List of the basis with the integer value of that binary representation.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::list<uint32_t> Read_BasisOp_BinaryRepresentation(std::string Basis_binary_filename, unsigned int n)
{
  list<uint32_t> Basis_li;

  ifstream myfile (Basis_binary_filename.c_str());
  string line, line2;
  
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n); //take the n first characters of line
      Basis_li.push_back(stoi(line2, 0, 2)); //convert string line2 into a binary integer
    }
    myfile.close();
  }
  
  return Basis_li;
}

/******************************************************************************/
/*** VERSION b) Operators are written as the integer values of the binary *****/
/****           representation of the interactions           ******************/
/******************************************************************************/
//' Read basis operators in integer representation from input file.
//' 
//' Reads the Basis operators in integer representation from the given input 
//' file. The operators should be written in one single column and as the
//' integer value of that binary representation. Returns a 
//' list of the given input file.
//' 
//' @param Basis_integer_filename The name of input file including path.
//' 
//' @return Basis_li List of the given input file.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::list<uint32_t> Read_BasisOp_IntegerRepresentation(std::string Basis_integer_filename)
{
  uint32_t Op = 0;
  list<uint32_t> Basis_li;

  ifstream myfile (Basis_integer_filename.c_str());
  string line;    

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      Op = stoi(line);
      Basis_li.push_back(Op);
    }
    myfile.close();
  }

  return Basis_li;
}

/******************************************************************************/
/*************************    Original Basis     ******************************/
/******************************************************************************/
//' Original Basis
//' 
//' Return the original basis, i.e., {s1, s2, ..., sn}
//' 
//' @param n Number of binary (spin) variables of the data set.
//' 
//' @return Basis_li List of the original basis i.e. {s1, s2, ..., sn}.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
std::list<uint32_t> Original_Basis(unsigned int n)
{
  uint32_t Op = 1;
  list<uint32_t> Basis_li;

  for (unsigned int i=0; i<n; i++)
  {
    Basis_li.push_back(Op);
    Op = Op << 1;
  }

  return Basis_li;
}

/******************************************************************************/
/***************************    Print Basis     *******************************/
/******************************************************************************/
//' Print Basis
//' 
//' Print the basis info in the terminal.
//' 
//' @param Basis_li A list containing the basis.
//' @param n Number of binary (spin) variables of the data set.
//' 
//' @export
// [[Rcpp::export(rng=false)]]
void PrintTerm_Basis(std::list<uint32_t> Basis_li, unsigned int n)
{
  int i = 1;
  for (list<uint32_t>::iterator it = Basis_li.begin(); it != Basis_li.end(); it++)
  {
    Rcpp::Rcout << "##\t " << i << " \t " << (*it) << " \t " << int_to_bstring(*it, n) << endl; i++;
  } Rcpp::Rcout << "##" << endl;
}
