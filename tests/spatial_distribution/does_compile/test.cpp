/*
  Checks that the file does compile
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "spatial_distribution.hpp"

using namespace std;

int main()
{
  // Write data to file
  ofstream f;
  f.open("res.txt");
  f << 0;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
