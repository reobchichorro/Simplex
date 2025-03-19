#include <string>
#include <vector>

#include "Simplex.h"
#include "mpsReader.h"

// void Simplex::CreateG() {
//     G = MatrixXd(instance.m+1, instance.n+1);

// }

void Simplex::SetInitialSolution() {

    // x_b = VectorXd::LinSpaced(instance.m, instance.n - instance.m, instance.n - 1);
    // x_n = VectorXd::LinSpaced(instance.n - instance.m, 0, instance.n - instance.m - 1);
    x_b = std::vector<int>(instance.m);
    x_n = std::vector<int>(instance.n);
    c_b = VectorXd::Zero(instance.m);
    c_n = VectorXd::Zero(instance.n);
    B = MatrixXd::Zero(instance.m, instance.m);
    A_n = MatrixXd::Zero(instance.m, instance.n-instance.m);
    Z = 0.0;
    // auto AA = instance.A;

    for (int j=0; j<instance.m; j++) {
        x_b[j] = instance.n - instance.m + j;
        c_b(j) = instance.c(x_b[j]);
        B(Eigen::all, j) = instance.A(Eigen::all, x_b[j]);
    }
    for (int j=0; j<instance.n - instance.m; j++) {
        x_n[j] = j;
        c_n(j) = instance.c(x_n[j]);
        A_n(Eigen::all, j) = instance.A(Eigen::all, x_n[j]);
    }

    B_inv = B.inverse();
}

bool Simplex::yan_cn(int& enteringVar, VectorXd& yan){
    enteringVar = -1;
    double biggest = 0.0;
    bool ans = false;
    for (int j=0; j<c_n.size(); j++){
        if (yan(j) < c_n(j))
            ans = true;
        if (biggest < yan(j)){
            enteringVar = j;
            biggest = yan(j);
        }
    }
    return ans;
}

// void Simplex::Standard() {

// }

void Simplex::RevisedNaive() {
    // Solve the system yB = c_B
    // Choose an entering column. This may be any column a of A_N such that ya is less than the corresponding component of c_N. If there is no such column, then the current solution is optimal.
    // Solve the system Bd = a
    // Find the largest t such that x_B^* - td >= 0. If there is no such t, then the problem in unbounded; otherwise, at least one component of x_B^* - td equals zero and the corresponding variable is leaving the basis
    // Set the value of the entering variable at t and replace the values x_B^* of the basic variables by x_B^* - td. Replace the leaving column of B by the entering column and, in the basis heading, replace the leaving variable by the entering variable.

    SetInitialSolution(); // 2.0

    VectorXd y = c_b.transpose()*B_inv; // 2.1
    
    VectorXd yan = y.transpose()*A_n; // 2.2
    std::cout << yan << std::endl;

    // std::cout << (yan < c_n) << std::endl;
    // for (auto it = yan.begin(); it != yan.end(); it++)
    //     std::cout << *it << " ";
    // std::cout << std::endl;

    int enteringVar = -1;
    int exitVar = -1;

    while(yan_cn(enteringVar, yan)) { // 2.2.5
        std::cout << enteringVar << " " << yan(enteringVar) << std::endl;
        break;

        VectorXd d = B_inv*A_n(Eigen:all, enteringVar);
        VectorXd t = 
    }

    // while(2.2.5){
        //2.3
        //2.4
        //2.5
        //2.1
        //2.2
    // }
}