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
    void RevisedNaive(); // Deprecated
    void Revised();
private:
    std::vector<int> x_b; // idx of basic variables
    std::vector<int> x_n; // idx of non-basic variables

    VectorXd c_b; // objcoef of basic variables
    VectorXd c_n; // objcoef of non-basic variables

    MatrixXd B; // constraint coefs of basic variables
    MatrixXd A_n; // constraint coefs of non-basic variables

    VectorXd lb; // lower bound of variables
    VectorXd ub; // upper bound of variables

    double Z; // obj function value
    VectorXd ans; // value of variables

    VectorXd b; // Deprecated

    MatrixXd B_inv; // inverse of B
    
    // MatrixXd G;
    // void CreateG();

    void SetInitialSolution();
    bool CheckBounds(); // Returns true if firstPhase is necessary
    int SelectEnteringVar(int& enteringVar, VectorXd& yan);
};
#endif
