#include "Matrix.h"
Matrix<double> matrix;
void slauTest(){
    std::cout<<"Testing slau solving enter size and vector\n";
    std::vector<double>b;
    int n;
    std::cin>>n;
    for (int i = 0; i < n; ++i) {
        double elem;
        std::cin>>elem;
        b.push_back(elem);
    }try {
        Matrix<double> gauss_res = matrix.gauss(b);
        Matrix<double> cramer_res = matrix.crammer(b);
        std::cout<<gauss_res<<'\n';
        std::cout << matrix * gauss_res.transpose()<<'\n';
        std::cout<<cramer_res<<'\n';
        std::cout << matrix * cramer_res.transpose()<<'\n';
    }
    catch (std::runtime_error &e){
        std::cout<<e.what()<<'\n';
        slauTest();
    }
}
void checkReverse(){
    std::cout<<"Checking reverse matrix"<<'\n';
    try{std::cout<<matrix.reverseMatrix()<<'\n';
        std::cout<<matrix*matrix.reverseMatrix();}
    catch(std::runtime_error&e){
        std::cout<<e.what();
        std::cout.flush();
        matrix = {};
        std::cout<<"Enter deminition and matrix\n";
        std::cin>>matrix;
        checkReverse();
    }
}
void operationsCheck(){
    std::cout<<"Checking matrix operations enter matrix";
    Matrix<double> matrix1;
    std::cin>>matrix1;
    try {
        std::cout<<matrix+matrix1<<'\n';
        std::cout<<matrix-matrix1<<'\n';
    }
    catch (std::runtime_error&e){
        std::cout<<e.what()<<'\n';
        operationsCheck();
    }
}

int main() {
    std::cout<<"Enter deminition and matrix\n";
    std::cin>>matrix;
    checkReverse();
    slauTest();
    operationsCheck();
    return 0;
}
