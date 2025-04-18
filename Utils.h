#include <iostream>
#include <string>
#include <vector>

#include "mpsReader.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/src/Core/Matrix.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef std::string str;
#define eps 1e-6

void createDualInstance(const mpsReader& instance, mpsReader& dual_instance);
