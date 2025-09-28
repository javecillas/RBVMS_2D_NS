#include <iostream>
#include <Eigen/Dense>
#include "readsolmech.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// int main() {
//     // Create material
//     Elas2D steel(1, 1.1, 0.3);
//     // Create nodes
//     // Create elements
//     // Solve the problem
//     // Define a 2x2 matrix A
//     MatrixXd A(2, 2);
//     A(0, 0) = 2; A(0, 1) = -1;
//     A(1, 0) = -3; A(1, 1) = 4;

//     // Define a vector b
//     VectorXd b(2);
//     b(0) = 1;
//     b(1) = 7;

//     // Solve Ax = b
//     VectorXd x = A.colPivHouseholderQr().solve(b);

//     // Output the solution
//     std::cout << "The solution x is:\n" << x << std::endl;

//     return 0;
// }

int main() {
    // Read input file
    char in_file_name[51];
    std::cout << "Enter the name of the input file: ";
    std::cin >> in_file_name;
    std::cout << "The input file name is: " << in_file_name << std::endl;
    ReadSolMech rsm(in_file_name);

    return 0;
}