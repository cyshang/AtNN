#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include "Parameter.h"

#define DEBUG_MODE
#define OUTPUT_TO_SCREEN

#ifndef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
#endif

#define PI 3.141592653589793

using std::endl;
using std::cout;
using std::vector;
using std::string;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::ostringstream;

extern ofstream	debug;
extern Parameter parameter;


#endif // !GLOBAL_H
