#include "Simplex.h"

Simplex::Simplex(mpsReader& instance) : instance(instance) {
    E_k = std::vector<std::pair<int, VectorXd> >(maxRefact);
    // EE = std::vector<VectorXd>();
    nll = (double *)NULL;
    A = instance.A.sparseView();

    (void)umfpack_di_symbolic(x_b.size(), x_b.size(), B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, nll, nll);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, nll, nll);
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

    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);

    (void)umfpack_di_symbolic(x_b.size(), x_b.size(), B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, nll, nll);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, nll, nll);
}

VectorXd Simplex::BTRAN() {
    VectorXd v = c_curr(x_b);

    for (int i=E_k_size-1; i>=0; i--)
        v(E_k[i].first) = v.dot(E_k[i].second) / E_k[i].second(E_k[i].first);

    VectorXd y = VectorXd::Zero(x_b.size());

    (void)umfpack_di_solve(UMFPACK_At, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), y.data(), v.data(), Numeric, nll, nll);

    return y;
}

VectorXd Simplex::FTRAN(int entIdx) {
    VectorXd d;
    // std::cout << entIdx << " entIdx - ABounds " << A.size() << std::endl;
    VectorXd a = A.col(entIdx);
    
    // std::cout << "umfpack" << std::endl;
    (void)umfpack_di_solve(UMFPACK_At, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), d.data(), instance.A.col(entIdx).data(), Numeric, nll, nll);    
    
    // std::cout << "Ek" << std::endl;
    for (int i=0; i<E_k_size; i++) {
        double p = d(E_k[i].first) / E_k[i].second(E_k[i].first);
        
        d -= p * E_k[i].second;
        
        d(E_k[i].first) = p;

        // for (int j = 0; j < x_b.size(); j++)
        // {
        //     if (j != E_k[i].first)
        //     v(j) -= p * E_k[i].second(j);
        // }
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
    // std::cout << ans.transpose() << std::endl;

    // B_inv = B.inverse();
}

bool Simplex::CheckBounds() {
    bool firstPhase = false;
    for (int j=0; j<instance.m; j++) {
        if (ans(x_b[j]) < lb(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_curr(x_n[jj]) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_curr(x_b[jj]) = 0.0;
                firstPhase = true;
                std::cout << j << " xj: " << x_b[j] << " " << ans(x_b[j]) << " " << lb(x_b[j]) << "\n";
            }
            c_curr(x_b[j]) = 1;
        }
        else if (ans(x_b[j]) > ub(x_b[j])) {
            if (!firstPhase) {
                for (int jj=0; jj<instance.n - instance.m; jj++)
                    c_curr(x_n[jj]) = 0.0;
                for (int jj=0; jj<j; jj++)
                    c_curr(x_b[jj]) = 0.0;
                firstPhase = true;
            }
            c_curr(x_b[j]) = -1;
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
    for (int j=0; j<x_n.size(); j++) {
        if (objImpr(x_n[j]) > eps && ans(x_n[j]) < ub(x_n[j])) {
            enteringVar = j;
            return 1;
        }
        else if (objImpr(x_n[j]) < eps && ans(x_n[j]) > lb(x_n[j])) {
            enteringVar = j;
            return -1;
        }
    }
    return 0;
}

void Simplex::Revised() {
    SetInitialSolution(); // 2.0
    bool firstPhase = CheckBounds();
    std::cout << "ans\n"  << ans.transpose() << std::endl;

    VectorXd y = BTRAN(); // 2.1
    
    VectorXd yA = y.transpose()*A; // 2.2
    // std::cout << "yan\n" << yan.transpose() << std::endl;

    int enteringVar = -1;
    int entidx = -1;
    int exitVar = -1;
    int exitidx = -1;
    int count = 0;
    int direction = SelectEnteringVar(enteringVar, yA);
    VectorXd dd;
    
    while (direction) {
        count++;
        std::cout << "Entering var: " << x_n[enteringVar]+1 << std::endl;
        
        entidx = x_n[enteringVar];
        std::cout << "FTRAN" << std::endl;
        VectorXd d = FTRAN(entidx); // 2.3
        std::cout << "FTRANOUT" << std::endl;
        
        std::cout << "d\n" << d.transpose() << std::endl;
        
        std::pair<int, VectorXd> E_ = {enteringVar, d};
        addEk(E_);
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
            if (isnan(tdub(j)) || d(j) == 0.0)
                tdub(j) = numeric_limits<double>::infinity();
        
        std::cout << "t\n" << tlb << " " << tub << "\n" << tdlb.transpose() << "\n" << tdub.transpose() << std::endl;
        
        int lowestUBIdx = std::distance(tdub.begin(), std::min_element(tdub.begin(), tdub.end()));
        
        std::cout << lowestUBIdx << " lubidx\n";
        std::cout << tdub.hasNaN() << " tlb tub " << tub << " " << lowestUBIdx << " " << tdub(lowestUBIdx) << "\n";
        std::cout << ans(x_b[1]) << " " << ub(x_b[1]) << " " << d(1) << " " << direction << "\n";
        
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
            std::cout << "ans\n" << ans.transpose() << std::endl;
        }
        else {
            double t = tdub(lowestUBIdx);
            exitVar = lowestUBIdx;
            std::cout << "Exiting var: " << x_b[exitVar]+1 << std::endl;
            exitidx = x_b[exitVar];
            ans(entidx) += direction*t;
            ans(x_b) += -direction*t*d;
            if (ans.hasNaN()) {
                std::cout << count << " nan\n";
                return;
            }

            std::cout << "ans\n"  << ans.transpose() << std::endl;

            // std::swap(c_b(exitVar), c_n(enteringVar));
            // B.col(exitVar).swap(A_n.col(enteringVar));
            // B_inv = B.inverse();
            std::swap(x_b[exitVar], x_n[enteringVar]);
        }

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
            }
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
