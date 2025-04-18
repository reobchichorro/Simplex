#include "Utils.h"

void createDualInstance(const mpsReader& instance, mpsReader& dual_instance) {
    int nm = instance.n-instance.m;
    std::vector<str> varNames;
    for (int i=0; i<instance.b.size(); i++)
        varNames.push_back("Cons" + std::to_string(i));
    std::vector<bool> GtZeroVars = std::vector<bool>(instance.lb.size(), false);
    std::vector<std::pair<int, bool>> BoundVars;

    for (int j=0; j<nm; j++) {
        if (instance.lb(j) > -eps)
            GtZeroVars[j] = true;
        
        if (instance.lb(j) > -numeric_limits<double>::infinity() && abs(instance.lb(j)) > eps) {
            varNames.push_back("LBx" + std::to_string(j));
            BoundVars.push_back({j, 0});
        }
        else
            varNames.push_back("Unused");

        if (instance.ub(j) < numeric_limits<double>::infinity() && abs(instance.ub(j)) > eps) {
            varNames.push_back("UBx" + std::to_string(j));
            BoundVars.push_back({j, 1});
        }
        else
            varNames.push_back("Unused");
    }

    dual_instance.b = VectorXd::Zero(nm);
    dual_instance.c = VectorXd::Zero(BoundVars.size() + instance.n);
    dual_instance.lb = VectorXd::Zero(BoundVars.size() + instance.n);
    dual_instance.ub = VectorXd::Zero(BoundVars.size() + instance.n);
    dual_instance.A = MatrixXd::Zero(nm, BoundVars.size() + instance.n);
    for (int j=0; j<instance.m; j++) {
        for (int i=0; i<nm; i++)
            dual_instance.A(i, j) = instance.A(j, i);
        if (instance.lb(nm+j) == -numeric_limits<double>::infinity()) {
            dual_instance.c(j) = instance.ub(nm+j);
            dual_instance.lb(j) = 0;
            dual_instance.ub(j) = numeric_limits<double>::infinity();
        }
        else if (instance.ub(nm+j) == numeric_limits<double>::infinity()) {
            dual_instance.c(j) = instance.lb(nm+j);
            dual_instance.lb(j) = -numeric_limits<double>::infinity();
            dual_instance.ub(j) = 0;
        }
        else {
            dual_instance.c(j) = instance.ub(nm+j);
            dual_instance.lb(j) = -numeric_limits<double>::infinity();
            dual_instance.ub(j) = numeric_limits<double>::infinity();
        }
    }
    
    for (size_t i=0; i<BoundVars.size(); i++) {
        int mi = instance.m+i;
        std::cout << i << ": " << BoundVars[i].first << " " << BoundVars[i].second << std::endl;
        dual_instance.A(BoundVars[i].first, mi) = 1;
        dual_instance.c(mi) = BoundVars[i].second ? instance.ub(BoundVars[i].first) : instance.lb(BoundVars[i].first);
        dual_instance.lb(mi) = BoundVars[i].second ? 0 : -numeric_limits<double>::infinity();
        dual_instance.ub(mi) = BoundVars[i].second ? numeric_limits<double>::infinity() : 0;
    }
    // std::cout << "bclbubA" << std::endl;
    for (int i=0; i<nm; i++) {
        dual_instance.A(i, instance.m+BoundVars.size()+i) = -1;
        dual_instance.lb(instance.m+BoundVars.size()+i) = instance.c(i);
        dual_instance.ub(instance.m+BoundVars.size()+i) = GtZeroVars[i] ? numeric_limits<double>::infinity() : instance.c(i);
    }
    dual_instance.c = -dual_instance.c;
    dual_instance.m = dual_instance.A.rows();
    dual_instance.n = dual_instance.A.cols();

    std::cout << "A" << std::endl;
    std::cout << dual_instance.A << std::endl;
    std::cout << "b" << std::endl;
    std::cout << dual_instance.b << std::endl;
    std::cout << "lb" << std::endl;
    std::cout << dual_instance.lb.transpose() << std::endl;
    std::cout << "ub" << std::endl;
    std::cout << dual_instance.ub.transpose() << std::endl;
    std::cout << "c" << std::endl;
    std::cout << dual_instance.c.transpose() << std::endl;
}
