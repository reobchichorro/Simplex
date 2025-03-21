// #include <iostream>
// #include <iomanip>
#include <string>
// #include <sstream>
#include <vector>
// #include <list>
// #include <algorithm>
// #include <unordered_set>
// #include <unordered_map>
// #include <map>
// #include <fstream>
// #include <math.h>
// #include <filesystem>
 
#include "mpsReader.h"
#include "eigen/Eigen/src/Core/Matrix.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
 
#define eps 1e-6
#define all(v) v.begin(),v.end()
typedef std::string str;
 
#ifndef __SIMPLEX_
#define __SIMPLEX_
class Simplex {
public:
    mpsReader instance;
 
    Simplex(mpsReader& instance) : instance(instance) {};

    void Standard();
    void RevisedNaive();
private:
    std::vector<int> x_b; // idx of basic variables
    std::vector<int> x_n; // idx of non-basic variables

    VectorXd c_b; // objcoef of basic variables
    VectorXd c_n; // objcoef of non-basic variables

    MatrixXd B; // constraint coefs of basic variables
    MatrixXd A_n; // constraint coefs of non-basic variables

    double Z; // obj function value
    VectorXd b; // value of basic variables

    MatrixXd B_inv; // inverse of B
    
    // MatrixXd G;
    // void CreateG();


    void SetInitialSolution();
    bool yan_cn(int& enteringVar, VectorXd& yan);
};
#endif
