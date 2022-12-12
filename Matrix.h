#include <bits/stdc++.h>
#define EPS 0.0000001

template<typename T> class Matrix;

template <typename T>
std::istream &operator>>(std::istream &is,Matrix<T> &matrix) {
    if (matrix.data.empty()) {
        int n, m;
        is >> n >> m;
        matrix.mkVector(n, m);
    }
    for (auto &i: matrix.data) {
        for (int j = 0; j < matrix.data[0].size(); ++j) {
            is >> i[j];
        }
    }
    return is;
}
template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix) {
    for (auto &i: matrix.data) {
        for (auto j: i) {
            os << j<<' ';
        }
        os<<'\n';
    }
    return os;
}

template <typename T>

class Matrix {
public:
    Matrix();

    Matrix(int n, int m);

    explicit Matrix(int n);

    explicit Matrix(std::vector<std::vector<T>>&input_data);
    explicit Matrix(std::vector<T>&input_data);
    ~Matrix();
    T det();
    Matrix gauss(const std::vector<T>&b);
    Matrix crammer(const std::vector<T>&b);
    Matrix transpose();
    Matrix operator+(const Matrix&matrix);
    Matrix operator-(const Matrix&matrix);
    Matrix operator*(const Matrix&matrix);
    Matrix operator*(T c);
    Matrix &operator=(const Matrix &matrix);
    T algebraicCompliment(int i,int j);
    Matrix reverseMatrix();

    friend std::istream & operator >> <T>(std::istream &is,Matrix<T>&matrix);
    friend std::ostream & operator << <T>(std::ostream &os,const Matrix<T>& matrix);
private:
    void mkVector(int n, int m);
    void triangulateMatrix();

    std::vector<std::vector<T>> data;
    std::vector<std::vector<T>> triangulated_data;
    std::vector<int> where;
};



template<typename T> Matrix<T>::Matrix() = default;

template<typename T> Matrix<T>::Matrix(int n, int m) {
    mkVector(n, m);
}

template<typename T> Matrix<T>::Matrix(int n) {
    mkVector(n, n);
}


template<typename T> Matrix<T>::Matrix(std::vector<std::vector<T>> &input_data) {
    data.assign(input_data.begin(),input_data.end());
    triangulated_data.assign(data.begin(),data.end());
}

template<typename T> Matrix<T>::Matrix(std::vector<T>&input_data) {
    data.resize(1);
    data[0].assign(input_data.begin(),input_data.end());
}

template<typename T> void Matrix<T>::mkVector(int n, int m) {
    data.resize(n);
    for (int i = 0; i < n; ++i) {
        data[i].resize(m);
    }
    triangulated_data.assign(data.begin(),data.end());
}

template<typename T> T Matrix<T>::det() {
    if (data.size() != data[0].size()) {
        throw std::runtime_error("Invalid dimension of the matrix");
    } else {
        triangulateMatrix();
        T ans = 1;
        for (int i = 0; i < data.size(); ++i) {
            ans*=triangulated_data[i][i];
        }
        return ans;
    }
}


template<typename T> void Matrix<T>::triangulateMatrix() {
    int n = (int) data.size();
    int m = (int) data[0].size();
    triangulated_data.assign(data.begin(),data.end());
    where.assign(m, -1);
    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i)
            if (abs(triangulated_data[i][col]) > abs(triangulated_data[sel][col]))
                sel = i;
        if (abs(triangulated_data[sel][col]) < EPS)
            continue;
        for (int i = col; i <= m; ++i) {
            T t = triangulated_data[sel][i];
            triangulated_data[sel][i] = triangulated_data[row][i];
            triangulated_data[row][i] = t;
        }
        where[col] = row;

        for (int i = 0; i < n; ++i)
            if (i != row) {
                T c = triangulated_data[i][col] / triangulated_data[row][col];
                for (int j = col; j <= m; ++j)
                    triangulated_data[i][j] -= triangulated_data[row][j] * c;
            }
        ++row;
    }
}

template<typename T> Matrix<T> Matrix<T>::gauss(const std::vector<T> &b) {

    if(b.size()!=data.size()){
        throw std::runtime_error("INVALID");
    }
    int n = (int) data.size();
    int m = (int) data[0].size();
    std::vector<T>ans;
    ans.assign (m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = b[where[i]] / triangulated_data[where[i]][i];
    for (int i=0; i<n; ++i) {
        T sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * triangulated_data[i][j];
        if (abs (sum - b[i]) > EPS)
            return {};
    }
    return Matrix(ans);
}

template<typename T>Matrix<T> Matrix<T>::operator+(const Matrix<T> &matrix) {
    if(data.size()==matrix.data.size()&& matrix.data[0].size()==data[0].size()){
        std::vector<std::vector<T>>res(data.size());
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[0].size(); ++j) {
                res[i].push_back(data[i][j]+matrix.data[i][j]);
            }
        }
        return Matrix(res);
    }
    throw std::runtime_error("Invalid dimension of the matrix");
}

template<typename T> Matrix<T> Matrix<T>::operator-(const Matrix<T> &matrix) {
    if(data.size()==matrix.data.size()&& matrix.data[0].size()==data[0].size()){
        std::vector<std::vector<T>>res(data.size());
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[0].size(); ++j) {
                res[i].push_back(data[i][j]-matrix.data[i][j]);
            }
        }
        return Matrix(res);
    }
    throw std::runtime_error("Invalid dimension of the matrix");
}

template<typename T> Matrix<T> Matrix<T>::operator*(const Matrix<T> &matrix) {
    if (data[0].size() != matrix.data.size())throw std::runtime_error("Invalid dimension of the matrix");
    else {
        std::vector<std::vector<T>>res(data.size());
        for (int i = 0; i < data.size(); ++i) {
            res[i].resize(matrix.data[0].size());
        }
        for(int i=0;i<data.size();++i){
            for(int j=0;j<matrix.data[0].size();++j){
                T elem=0;
                for(int k=0;k<data[0].size();++k){
                    elem+=data[i][k]*matrix.data[k][j];
                }
                res[i][j]=elem;
            }
        }
        return Matrix(res);
    }
}

template<typename T> Matrix<T> Matrix<T>::operator*(T c) {
    std::vector<std::vector<T>>res(data.size());
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[0].size(); ++j) {
            res[i].push_back(data[i][j]*c);
        }
    }
    return Matrix(res);
}

template<typename T> Matrix<T> &Matrix<T>::operator=(const Matrix<T> &matrix) {
    data.assign(matrix.data.begin(),matrix.data.end());
    triangulated_data.assign(matrix.triangulated_data.begin(),matrix.triangulated_data.end());
    where.assign(matrix.where.begin(),matrix.where.end());
    return *this;
}

template<typename T> Matrix<T> Matrix<T>::crammer(const std::vector<T>&b) {
    if (det()==0)throw std::runtime_error("determinant is 0");
    std::vector<T> res;
    if (data.size()!=data[0].size())throw std::runtime_error("Invalid dimension of the matrix");
    for (int i = 0; i < data[0].size(); ++i) {
        std::vector<std::vector<T>>new_data;
        new_data.assign(data.begin(),data.end());
        for (int j = 0; j < data.size(); ++j) {
            new_data[j][i]=b[j];
        }
        T determinant = Matrix(new_data).det();
        res.push_back(determinant/det());
    }


    return Matrix(res);
}

template<typename T> T Matrix<T>::algebraicCompliment(int i, int j) {
    std::vector<std::vector<T>>new_data(data.size()-1);
    for (int k = 0; k < data.size(); ++k) {

        if (k!=i){
            for (int l = 0; l < data[k].size(); ++l) {
                if(l!=j){
                    new_data[k-(k>i)].push_back(data[k][l]);
                }
            }
        }

    }
    return pow(-1,i+j)*Matrix(new_data).det();
}

template<typename T> Matrix<T> Matrix<T>::reverseMatrix() {
    if(det()==0)throw std::runtime_error("Determinant is 0");
    std::vector<std::vector<T>>added_data(data.size());
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[0].size(); ++j) {
            auto compliment = algebraicCompliment(j,i);
            added_data[i].push_back(compliment);
        }
    }
    return Matrix(added_data)*(1/det());
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    auto ans = Matrix(data[0].size(),data.size());
    for (int i = 0; i < data[0].size(); ++i) {
        for (int j = 0; j < data.size(); ++j) {
            ans.data[i][j] = data[j][i];
        }
    }
    return ans;
}

template<typename T>
Matrix<T>::~Matrix() = default;

