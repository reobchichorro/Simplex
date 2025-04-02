#include <iostream>
#include <string>

#include "Utils.cpp"
#include "mpsReader.h"
#include "Simplex.h"

typedef std::string str;

int main(int argc, char **argv) {
    str filename = argv[1];
    str filetype;
    // getFileType(filename, filetype);

    filetype = filename.substr(filename.find('.'));

    mpsReader mps;

    int pp = 0; // TODO: parametrizar

    if (filetype == ".mps") {
        mps.read(filename, pp);

        std::cout << mps.Name << std::endl;
        std::cout << mps.n_rows << std::endl;
        std::cout << mps.n_rows_eq << std::endl;
        std::cout << mps.n_rows_inq << std::endl;
        std::cout << mps.n_cols << std::endl;
        std::cout << mps.m << std::endl;
        std::cout << mps.n << std::endl;
    
        std::cout << "A" << std::endl;
        std::cout << mps.A << std::endl;
        std::cout << "b" << std::endl;
        std::cout << mps.b << std::endl;
        std::cout << "lb" << std::endl;
        std::cout << mps.lb.transpose() << std::endl;
        std::cout << "ub" << std::endl;
        std::cout << mps.ub.transpose() << std::endl;
        std::cout << "c" << std::endl;
        std::cout << mps.c.transpose() << std::endl;

        for (std::vector<int>::size_type i=0; i < mps.restricoes.size(); i++)
            std::cout << mps.restricoes[i] << " ";
        std::cout << std::endl;
        {
            // std::cout << mps.row_labels << std::endl;
            // std::cout << mps.col_labels << std::endl;
            // std::cout << mps.row_list << std::endl;
            // std::cout << mps.col_list << std::endl;
            // std::cout << mps.restricoes << std::endl;
        
            // std::cout << preprocess << std::endl;

            // l = mps.lb;
            // u = mps.ub;
            // A_dense = mps.A;
            // b = mps.b;
            // c = mps.c;
            // m = mps.n_rows_eq + mps.n_rows_inq;
            // n = mps.n_cols + mps.n_rows_inq + mps.n_rows_eq;
        }
    }
    
    std::cout << "Simplex\n";
    Simplex solver(mps);
    solver.Revised();
}