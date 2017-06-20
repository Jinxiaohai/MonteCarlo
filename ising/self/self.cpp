#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstddef>
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
#define cout << "\033[32m"  Prompt
#define cout << "\033[31m"  Error

int GetMax(int, int);
int GetMin(int, int);

int main(int argc, char *argv[])
{
  /**
   * Usage: self    m    n    iterations    thresh    seed
   * M, N:  the number of iterations.
   * THRESH:  the threshhold.
   * SEED:   a seed for the random number generator
   */
  int *c1 = nullptr;
  int i = 0, iterations = 0, m = 0, n = 0，seed = 0;
  double thresh = 0.;
  string plotFileName = "isingFinal.txt";
  string pngFileName = "isingFinal.png";
  string title = "Final Configuration";
  /// 概率
  const unsigned int NUMBERPROB = 5;
  double prob[NUMBERPROB] = {0.98, 0.85, 0.50, 0.15, 0.02};

  
  return 0;
}
