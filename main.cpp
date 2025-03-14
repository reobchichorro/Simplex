#include <iostream>
#include <string>

#include "utils.cpp"
#include "mpsReader.h"

typedef std::string str;

int main(int argc, char **argv) {
    str filename = argv[1];
    str filetype;
    getFileType(filename, filetype);

    mpsReader mps;

    if (filetype == ".mps") {
        mps.read(filename, pp);

        // l = mps.lb;
        // u = mps.ub;
        // A_dense = mps.A;
        // b = mps.b;
        // c = mps.c;
        // m = mps.n_rows_eq + mps.n_rows_inq;
        // n = mps.n_cols + mps.n_rows_inq + mps.n_rows_eq;
    }
}