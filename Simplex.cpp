#include "Simplex.h"

Simplex::Simplex(mpsReader& instance) : instance(instance) {
    E_k = std::vector<std::pair<int, VectorXd> >(maxRefact);
    nll = (double *)NULL;
    A = instance.A.sparseView();
}

Simplex::~Simplex()
{
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
}

void Simplex::addEk(std::pair<int, VectorXd> E) {
    if (E_k_size >= maxRefact)
        refactor();

    E_k[E_k_size] = E;
    E_k_size++;
}

void Simplex::refactor() {
    E_k = std::vector<std::pair<int, VectorXd> >(maxRefact);
    E_k_size = 0;

    for (long unsigned int i=0; i<x_b.size(); i++) {
        B.col(i) = A.col(x_b[i]);
    }

    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);

    (void)umfpack_di_symbolic(x_b.size(), x_b.size(), B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, nll, nll);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, nll, nll);
}

VectorXd Simplex::BTRAN() {
    // std::cout << "c " << c_curr.transpose() << std::endl;
    // std::cout << "xb: ";
    // for (int i=0; i<x_b.size(); i++)
    //     std::cout << x_b[i] << " ";
    // std::cout << std::endl;
    VectorXd v = c_curr(x_b);
    double vdot = 0.0;

    // std::cout << "v " << v.transpose() << std::endl;
    for (int i=E_k_size-1; i>=0; i--){ 
        // std::cout << "BTRAN Eki " << E_k[i].first << std::endl;
        // std::cout << E_k[i].second.transpose() << std::endl;

        vdot = v.dot(E_k[i].second) - v(E_k[i].first)*E_k[i].second(E_k[i].first);
        v(E_k[i].first) -= vdot;
        v(E_k[i].first) /= E_k[i].second(E_k[i].first);
        // std::cout << "v" << i << " " << vdot << " " << v.transpose() << std::endl;
    }

    // std::cout << "v " << v.transpose() << std::endl;
    VectorXd y = VectorXd::Zero(x_b.size());

    (void)umfpack_di_solve(UMFPACK_At, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), y.data(), v.data(), Numeric, nll, nll);
    // std::cout << "y " << y.transpose() << std::endl;

    return y;
}

VectorXd Simplex::FTRAN(int entIdx) {
    VectorXd d(x_b.size());
    
    // std::cout << entIdx << " entIdx - a " << instance.A.col(entIdx).transpose() << std::endl;

    (void)umfpack_di_solve(UMFPACK_A, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), d.data(), instance.A.col(entIdx).data(), Numeric, nll, nll);
    
    // std::cout << "Mid-d " << d.transpose() << std::endl;
    for (int i=0; i<E_k_size; i++) {
        double p = d(E_k[i].first) / E_k[i].second(E_k[i].first);
        
        d -= p * E_k[i].second;
        
        d(E_k[i].first) = p;

        // std::cout << "d-calc " << i << " " << p << "\n\t" << d.transpose() << std::endl;
    }
    // std::cout << "ddddd" << std::endl;

    return d;
}

void Simplex::SetInitialSolution() {

    // x_b = VectorXd::LinSpaced(instance.m, instance.n - instance.m, instance.n - 1);
    // x_n = VectorXd::LinSpaced(instance.n - instance.m, 0, instance.n - instance.m - 1);
    x_b = std::vector<int>(instance.m);
    x_n = std::vector<int>(instance.n-instance.m);
    // c_b = VectorXd::Zero(instance.m);
    // c_n = VectorXd::Zero(instance.n-instance.m);
    c_curr = instance.c;
    // B = MatrixXd::Zero(instance.m, instance.m);
    // A_n = MatrixXd::Zero(instance.m, instance.n-instance.m);
    Z = 0.0;
    lb = instance.lb; // VectorXd::Zero(instance.n);
    ub = instance.ub; // VectorXd::Zero(instance.n);
    ans = VectorXd::Zero(instance.n);
    B = Eigen::SparseMatrix<double>(instance.m, instance.m);
    
    // Non-basic variables
    for (int j=0; j<instance.n - instance.m; j++) {
        x_n[j] = j;
        // c_n(j) = instance.c(x_n[j]);
        // A_n(Eigen::all, j) = A(Eigen::all, x_n[j]);
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
        // c_b(j) = instance.c(x_b[j]);
        B.col(j) = A.col(x_b[j]);
        ans(x_b[j]) = (instance.A(j, Eigen::seqN(0, instance.n - instance.m)).dot(ans(Eigen::seqN(0, instance.n - instance.m))));
    }

    (void)umfpack_di_symbolic(x_b.size(), x_b.size(), B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, nll, nll);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, nll, nll);
}

bool Simplex::CheckBounds() {
    bool firstPhase = false;
    lb = instance.lb;
    ub = instance.ub;
    for (int j=0; j<instance.m; j++) {
        if (ans(x_b[j]) < instance.lb(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_curr(x_n[jj]) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_curr(x_b[jj]) = 0.0;
                firstPhase = true;
                // std::cout << j << " xj: " << x_b[j] << " " << ans(x_b[j]) << " " << lb(x_b[j]) << "\n";
            }
            c_curr(x_b[j]) = 1;
            ub(x_b[j]) = lb(x_b[j]);
            lb(x_b[j]) = -numeric_limits<double>::infinity();
        }
        else if (ans(x_b[j]) > instance.ub(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_curr(x_n[jj]) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_curr(x_b[jj]) = 0.0;
                firstPhase = true;
            }
            c_curr(x_b[j]) = -1;
            lb(x_b[j]) = ub(x_b[j]);
            ub(x_b[j]) = numeric_limits<double>::infinity();
        }
        else if(firstPhase)
            c_curr(x_b[j]) = 0.0;
    }
    // std::cout << c_n << std::endl;
    // std::cout << c_b << std::endl;
    return firstPhase;
}

int Simplex::SelectEnteringVar(int& enteringVar, VectorXd& yA) {
    enteringVar = -1;
    VectorXd objImpr = c_curr-yA;

    // int tprint = 14;
    // std::cout << "entering14 " << x_n[tprint] << ", " << c_curr(x_n[tprint]) << ", " << yA(x_n[tprint]) << ", " << lb(x_n[tprint]) << ", " <<  ans(x_n[tprint]) << ", " << ub(x_n[tprint]) << "\n";

    // std::cout << c_curr.transpose() << "\n" << yA.transpose() << std::endl;
    for (long unsigned int j=0; j<x_n.size(); j++) {
        // if (j==tprint)
            // std::cout << "objImpr14 " << objImpr(x_n[tprint]) << ", " << ans(x_n[tprint]) << ", " << ub(x_n[tprint]) << "\n";
        if (objImpr(x_n[j]) > eps && ans(x_n[j]) < ub(x_n[j])) {
            enteringVar = j;
            // std::cout << "enteringidx " << j << "\n";
            return 1;
        }
        else if (objImpr(x_n[j]) < -eps && ans(x_n[j]) > lb(x_n[j])) {
            enteringVar = j;
            // std::cout << "enteringidx " << j << "\n";
            return -1;
        }
    }
    return 0;
}

void Simplex::Revised() {
    SetInitialSolution(); // 2.0
    bool firstPhase = CheckBounds();
    // std::cout << "phase "  << firstPhase << std::endl;
    // int pcount =0, qcount=0;
    // for (int i=0; i<c_curr.size(); i++) {
    //     if (abs(c_curr[i]-1) <= eps)
    //         pcount++;
    //     if (abs(c_curr[i]+1) <= eps)
    //         qcount++;
    // }
    // std::cout << pcount << " PQ " << qcount << "\n";
    // std::cout << "ans\n"  << ans.transpose() << std::endl;

    VectorXd y = BTRAN(); // 2.1
    
    VectorXd yA = y.transpose()*A; // 2.2
    // std::cout << "yan\n" << yan.transpose() << std::endl;

    int enteringVar = -1;
    int entidx = -1;
    int exitVar = -1;
    // int exitidx = -1;
    int count = 0;
    int direction = SelectEnteringVar(enteringVar, yA);
    VectorXd dd;
    
    while (direction && count < 1e300) {
        count++;
        // if (!firstPhase)
        //     std::cout << "Entering var: " << x_n[enteringVar]+1 << " direction " << direction << std::endl;
        
        entidx = x_n[enteringVar];
        VectorXd d = FTRAN(entidx); // 2.3
        
        // std::cout << "d\n" << d.transpose() << std::endl;
        
        // std::cout << "Ek" << std::endl;
        // std::cout << E_k_size << " E_k " << E_k[E_k_size-1].first << "\n" << E_k[E_k_size-1].second.transpose() << "\n";

        double tlb = direction * (lb(entidx) - ans(entidx));
        double tub = direction * (ub(entidx) - ans(entidx));
        
        if (direction == -1)
            std::swap(tlb, tub);
        
        VectorXd tdlb = direction * ((ans(x_b) - lb(x_b)).array() / d.array()).matrix();
        VectorXd tdub = direction * ((ans(x_b) - ub(x_b)).array() / d.array()).matrix(); // 2.4
        
        if (direction == 1)
            std::swap(tdlb, tdub);
        
        for (int j=0; j<d.size(); j++)
            if (d(j) < 0.0)
                std::swap(tdlb(j), tdub(j));

        for (int j=0; j<tdub.size(); j++)
            if (isnan(tdub(j)) || abs(d(j)) <= eps)
                tdub(j) = numeric_limits<double>::infinity();
        
        // std::cout << "t\n" << tlb << " " << tub << "\n" << tdlb.transpose() << "\n" << tdub.transpose() << std::endl;
        
        int lowestUBIdx = std::distance(tdub.begin(), std::min_element(tdub.begin(), tdub.end()));
        
        // int tprint = 19;
        // std::cout << "step19 " << tdub(tprint) << " " << lb(x_b[tprint]) << ", " <<  ans(x_b[tprint]) << ", " << ub(x_b[tprint]) << ", " << d(tprint) << "\n";

        // std::cout << lowestUBIdx << " lubidx\n";
        // std::cout << tdub.hasNaN() << " tlb tub " << tub << " " << lowestUBIdx << " " << tdub(lowestUBIdx) << "\n";
        // std::cout << ans(x_b[1]) << " " << ub(x_b[1]) << " " << d(1) << " " << direction << "\n";
        
        if (tdub(lowestUBIdx) == numeric_limits<double>::infinity() && tub == numeric_limits<double>::infinity()) {
            std::cout << "Unbounded\n";
            return;
        }

        if (tub <= tdub(lowestUBIdx)) {
            double t = tub;
            ans(entidx) += direction*t;
            ans(x_b) += -direction*t*d;
            if (ans.hasNaN()) {
                std::cout << count << " nan\n";
                return;
            }
            // std::cout << "ans\n" << ans.transpose() << std::endl;
        }
        else {
            double t = tdub(lowestUBIdx);
            exitVar = lowestUBIdx;
            // std::cout << "Exiting var: " << x_b[exitVar]+1 << std::endl;
            // std::cout << "td: " << direction << " " << t << " " << d(exitVar) << std::endl;
            // exitidx = x_b[exitVar];
            ans(entidx) += direction*t;
            ans(x_b) += -direction*t*d;
            if (ans.hasNaN()) {
                std::cout << count << " nan\n";
                return;
            }

            std::pair<int, VectorXd> E_ = {exitVar, d};
            addEk(E_);

            // std::cout << "ans\n"  << ans.transpose() << std::endl;

            // std::swap(c_b(exitVar), c_n(enteringVar));
            // B.col(exitVar).swap(A_n.col(enteringVar));
            // B_inv = B.inverse();
            std::swap(x_b[exitVar], x_n[enteringVar]);
        }
        
        // std::cout << "\nNext iteration\n";
        std::cout << "cost: " << instance.c.transpose() * ans << " 1a fase? " << firstPhase << std::endl;

        if (firstPhase) {
            firstPhase = CheckBounds();
            if (!firstPhase) {
                std::cout << "Second phase\n";
                // Non-basic variables
                // for (int j=0; j<instance.n - instance.m; j++)
                //     c_n(j) = instance.c(x_n[j]);
                // // Basic variables
                // for (int j=0; j<instance.m; j++)
                //     c_b(j) = instance.c(x_b[j]);
                c_curr = instance.c;
                refactor();
            }
            // else {
            //     pcount =0, qcount=0;
            //     for (int i=0; i<c_curr.size(); i++) {
            //         if (abs(c_curr[i]-1) <= eps)
            //             pcount++;
            //         if (abs(c_curr[i]+1) <= eps)
            //             qcount++;
            //     }
            //     std::cout << pcount << " PQ " << qcount << "\n";
            // }
        }

        y = BTRAN(); // 2.1
        yA = y.transpose()*A; // 2.2
        // std::cout << "yan\n" << yan.transpose() << std::endl;

        // std::cout << "cb\n" << c_b.transpose() << std::endl;
        // std::cout << "cn\n" << c_n.transpose() << std::endl;

        // std::cout << "B\n" << B << std::endl;
        // std::cout << "A_n\n" << A_n << std::endl << std::endl;

        direction = SelectEnteringVar(enteringVar, yA);
    }

    std::cout << count << std::endl;
    std::cout << ans.transpose() << std::endl;
    std::cout << instance.c.transpose() * ans << std::endl;
}

void Simplex::DualInstance() {
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
    std::cout << "bclbubA" << std::endl;
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
    std::cout << "bclbubA" << std::endl;
    for (size_t i=0; i<BoundVars.size(); i++) {
        int mi = instance.m+i;
        std::cout << i << ": " << BoundVars[i].first << " " << BoundVars[i].second << std::endl;
        dual_instance.A(BoundVars[i].first, mi) = 1;
        dual_instance.c(mi) = BoundVars[i].second ? ub(BoundVars[i].first) : lb(BoundVars[i].first);
        dual_instance.lb(mi) = BoundVars[i].second ? 0 : -numeric_limits<double>::infinity();
        dual_instance.ub(mi) = BoundVars[i].second ? numeric_limits<double>::infinity() : 0;
    }
    std::cout << "bclbubA" << std::endl;
    for (int i=0; i<nm; i++) {
        dual_instance.A(i, instance.m+BoundVars.size()+i) = -1;
        dual_instance.lb(instance.m+BoundVars.size()+i) = instance.c(i);
        dual_instance.ub(instance.m+BoundVars.size()+i) = GtZeroVars[i] ? numeric_limits<double>::infinity() : instance.c(i);
    }

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
