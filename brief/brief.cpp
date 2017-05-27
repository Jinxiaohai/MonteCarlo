#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <ctime>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <random>

using namespace std;
const unsigned int N = 30000000;
int main(int argc, char *argv[])
{
  static default_random_engine e(time(nullptr));
  /// calculate (0, 3);
  static uniform_real_distribution<double> u(0, 1.5);
  static uniform_real_distribution<double> u2(0, 1);
  double *randomNum = new double[N];
  double *randomValue = new double[N];
  unsigned int M = 0;
  for (unsigned int i = 0; i != N; ++i)
    {
      randomNum[i] = u(e);
      randomValue[i] = u2(e);
      if (pow(sin(1./randomNum[i]), 2) >= randomValue[i])
	{
	  ++M;
	}
    }
  double integrate = static_cast<double>(M)/N * (1.5-0);
  cout << integrate << endl;
  delete [] randomNum;
  delete [] randomValue;
  return 0;
}
