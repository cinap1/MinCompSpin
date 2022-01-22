/* 
 * Header for file containing support fuctions.
 */

#include <sstream>
#include <vector>
#include <map>
using namespace std;

// create a string of binary values from an int.
string int_to_bstring(unsigned x, unsigned int n);

// Counts all the set bits of an int.
unsigned int countSetBits(unsigned int n);

// Convert map to vector.
vector<std::vector<int>> mapToVec(map<uint32_t, unsigned int> old_map) ;
