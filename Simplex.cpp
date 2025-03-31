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
    lb = instance.lb; // VectorXd::Zero(instance.n);
    ub = instance.ub; // VectorXd::Zero(instance.n);
    ans = VectorXd::Zero(instance.n);
    // auto AA = instance.A;

    // Non-basic variables
    for (int j=0; j<instance.n - instance.m; j++) {
        x_n[j] = j;
        c_n(j) = instance.c(x_n[j]);
        A_n(Eigen::all, j) = instance.A(Eigen::all, x_n[j]);
        bool infLB = lb(j) == -numeric_limits<double>::infinity();
        bool infUB = ub(j) == numeric_limits<double>::infinity();
        if (!infLB)
            ans(j) = lb(j);
        else if (!infUB)
            ans(j) = ub(j);
    }
    // Basic variables
    for (int j=0; j<instance.m; j++) {
        x_b[j] = instance.n - instance.m + j;
        c_b(j) = instance.c(x_b[j]);
        B(Eigen::all, j) = instance.A(Eigen::all, x_b[j]);
        ans(x_b[j]) = (instance.A(j, Eigen::seqN(0, instance.n - instance.m)).dot(ans(Eigen::seqN(0, instance.n - instance.m))));
    }
    std::cout << ans << std::endl;

    B_inv = B.inverse();
}

bool Simplex::CheckBoundsOnInit() {
    bool firstPhase = false;
    for (int j=0; j<instance.m; j++) {
        if (ans(x_b[j]) < lb(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_n(jj) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_b(jj) = 0.0;
                firstPhase = true;
            }
            c_b(j) = 1;
        }
        else if (ans(x_b[j]) > ub(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_n(jj) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_b(jj) = 0.0;
                firstPhase = true;
            }
            c_b(j) = -1;
        }
        else if(firstPhase)
            c_b(j) = 0.0;
    }
    // std::cout << c_n << std::endl;
    // std::cout << c_b << std::endl;
    return firstPhase;
}

int Simplex::SelectEnteringVar(int& enteringVar, VectorXd& yan) {
    enteringVar = -1;
    VectorXd objImpr = c_n-yan;
    for (int j=0; j<objImpr.size(); j++) {
        if (objImpr(j) > .0 && ans(x_n[j]) < ub(x_n[j])) {
            enteringVar = j;
            return 1;
        }
        else if (objImpr(j) < .0 && ans(x_n[j]) > lb(x_n[j])) {
            enteringVar = j;
            return -1;
        }
    }
    return 0;
}

// void Simplex::Standard() {

// }

void Simplex::RevisedNaive() {
    std::cout << "Deprecated\n";
    return;

    SetInitialSolution(); // 2.0

    VectorXd y = c_b.transpose()*B_inv; // 2.1
    
    VectorXd yan = y.transpose()*A_n; // 2.2
    std::cout << "yan\n" << yan.transpose() << std::endl;

    int enteringVar = -1;
    int exitVar = 0;
    int count = 0;
    while(SelectEnteringVar(enteringVar, yan) && count < 10) { // 2.2.5
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


    std::cout << ans.transpose() << std::endl;
}

void Simplex::Revised() {
    SetInitialSolution(); // 2.0
    bool firstPhase = CheckBoundsOnInit();
    
    if (firstPhase) {
        // TODO
        
        // Non-basic variables
        for (int j=0; j<instance.n - instance.m; j++)
            c_n(j) = instance.c(x_n[j]);
        // Basic variables
        for (int j=0; j<instance.m; j++)
            c_b(j) = instance.c(x_b[j]);
    }

    VectorXd y = c_b.transpose()*B_inv; // 2.1
    
    VectorXd yan = y.transpose()*A_n; // 2.2
    std::cout << "yan\n" << yan.transpose() << std::endl;

    int enteringVar = -1;
    int exitVar = 0;
    int count = 0;
    int direction = SelectEnteringVar(enteringVar, yan);
    while(direction && count < 10) { // 2.2.5
        count++;
        std::cout << "Entering var: " << x_n[enteringVar]+1 << std::endl;
        
        VectorXd d = B_inv*A_n(Eigen::all, enteringVar); // 2.3
        std::cout << "d\n" << d.transpose() << std::endl;

        int entidx = x_n[enteringVar];

        double tlb = direction * (lb(entidx) - ans(entidx));
        double tub = direction * (ub(entidx) - ans(entidx));

        VectorXd tdlb = direction * ((ans(x_b) - lb(x_b)).array() / d.array()).matrix();
        VectorXd tdub = direction * ((ans(x_b) - ub(x_b)).array() / d.array()).matrix(); // 2.4

        std::cout << "t\n" << tlb << " " << tub << "\n" << tdlb.transpose() << "\n" << tdub.transpose() << std::endl;

        return;
        
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

        direction = SelectEnteringVar(enteringVar, yan);
    }


    std::cout << ans.transpose() << std::endl;
}
