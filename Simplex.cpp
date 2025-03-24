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
    x_n = std::vector<int>(instance.n-instance.m);
    c_b = VectorXd::Zero(instance.m);
    c_n = VectorXd::Zero(instance.n-instance.m);
    B = MatrixXd::Zero(instance.m, instance.m);
    A_n = MatrixXd::Zero(instance.m, instance.n-instance.m);
    Z = 0.0;
    b = instance.b;
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

bool Simplex::yan_cn(int& enteringVar, VectorXd& yan) {
    enteringVar = -1;
    bool ans = false;
    VectorXd objImpr = c_n-yan;
    for (int j=0; j<objImpr.size(); j++) {
        if (objImpr(j) > .0) {
            enteringVar = j;
            return true;
        }
    }
    return ans;
}

// void Simplex::Standard() {

// }

void Simplex::RevisedNaive() {
    SetInitialSolution(); // 2.0

    VectorXd y = c_b.transpose()*B_inv; // 2.1
    
    VectorXd yan = y.transpose()*A_n; // 2.2
    std::cout << "yan\n" << yan.transpose() << std::endl;

    int enteringVar = -1;
    int exitVar = 0;
    int count = 0;
    while(yan_cn(enteringVar, yan) && count < 10) { // 2.2.5
        count++;
        std::cout << "Entering var: " << x_n[enteringVar]+1 << std::endl;
        
        VectorXd d = B_inv*A_n(Eigen::all, enteringVar); // 2.3
        std::cout << "d\n" << d.transpose() << std::endl;

        VectorXd t = (b.array() / d.array()).matrix(); // 2.4
        std::cout << "t\n" << t.transpose() << std::endl;
        exitVar = 0;
        for (int j=1; j<t.size(); j++) {
            if (t(j) < t(exitVar) && t(j) > 0)
                exitVar = j;
        }
        std::cout << "Exiting var: " << x_b[exitVar]+1 << std::endl;

        b = b - t(exitVar)*d; // 2.5
        b(exitVar) = t(exitVar);

        std::swap(c_b(exitVar), c_n(enteringVar));
        
        B.col(exitVar).swap(A_n.col(enteringVar));
        B_inv = B.inverse();

        std::swap(x_b[exitVar], x_n[enteringVar]);

        y = c_b.transpose()*B_inv; // 2.1
        yan = y.transpose()*A_n; // 2.2
        std::cout << "yan\n" << yan.transpose() << std::endl;

        std::cout << "cb\n" << c_b.transpose() << std::endl;
        std::cout << "cn\n" << c_n.transpose() << std::endl;

        std::cout << "B\n" << B << std::endl;
        std::cout << "A_n\n" << A_n << std::endl << std::endl;
    }


    std::cout << b.transpose() << std::endl;
}