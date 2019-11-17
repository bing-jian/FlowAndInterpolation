// Author: Ce Liu (c) Dec, 2009; celiu@mit.edu
// Modified By: Deepak Pathak (c) 2016; pathak@berkeley.edu

#ifndef MATRIX_H_
#define MATRIX_H_

#include "Vector.h"
#include "project.h"
#include "stdio.h"
#ifdef _QT
#include <QFile>
#endif
#include <cstring>
#include <iostream>

using namespace std;

template <class T = double>
class Matrix {
 private:
  int nRow, nCol;
  double* pData;
  static bool IsDispInfo;

 public:
  Matrix(void);
  Matrix(int _nrow, int _ncol, double* data = NULL);
  Matrix(const Matrix<T>& matrix);
  ~Matrix(void);
  void releaseData();
  void copyData(const Matrix<T>& matrix);
  void allocate(const Matrix<T>& matrix);
  void allocate(int _nrow, int _ncol);
  void reset();
  bool dimMatch(const Matrix<T>& matrix) const;
  bool dimcheck(const Matrix<T>& matrix) const;
  void loadData(int _nrow, int _ncol, T* data);
  static void enableDispInfo(bool dispInfo = false) { IsDispInfo = dispInfo; };
  // display the matrix
  void printMatrix();
  void identity(int ndim);

  // function to access the member variables
  inline int nrow() const { return nRow; };
  inline int ncol() const { return nCol; };
  inline double* data() { return pData; };
  inline const double* data() const { return (const double*)pData; };
  inline double operator[](int index) const { return pData[index]; };
  inline double& operator[](int index) { return pData[index]; };
  inline double data(int row, int col) const {
    return pData[row * nCol + col];
  };
  inline double& data(int row, int col) { return pData[row * nCol + col]; };
  bool matchDimension(int _nrow, int _ncol) const {
    if (nRow == _nrow && nCol == _ncol)
      return true;
    else
      return false;
  };
  bool matchDimension(const Matrix<T>& matrix) const {
    return matchDimension(matrix.nrow(), matrix.ncol());
  };

  // functions to check dimensions
  bool checkDimRight(const Vector<T>& vector) const;
  bool checkDimRight(const Matrix<T>& matrix) const;
  bool checkDimLeft(const Vector<T>& vector) const;
  bool checkDimLeft(const Matrix<T>& matrix) const;

  // functions for matrix computation
  void Multiply(Vector<T>& result, const Vector<T>& vect) const;
  void Multiply(Matrix<T>& result, const Matrix<T>& matrix) const;

  void transpose(Matrix& result) const;
  void fromVector(const Vector<T>& vect);
  double norm2() const;
  double sum() const {
    double total = 0;
    for (int i = 0; i < nCol * nRow; i++) total += pData[i];
    return total;
  }
  // operators
  Matrix& operator=(const Matrix<T>& matrix);

  Matrix& operator+=(double val);
  Matrix& operator-=(double val);
  Matrix& operator*=(double val);
  Matrix& operator/=(double val);

  Matrix& operator+=(const Matrix<T>& matrix);
  Matrix& operator-=(const Matrix<T>& matrix);
  Matrix& operator*=(const Matrix<T>& matrix);
  Matrix& operator/=(const Matrix<T>& matrix);

  friend Vector<T> operator*<> (const Matrix<T>& matrix, const Vector<T>& vect);
  friend Matrix<T> operator*<> (const Matrix<T>& matrix1,
                             const Matrix<T>& matrix2);

  // solve linear systems
  void SolveLinearSystem(Vector<T>& result, const Vector<T>& b) const;
  void ConjugateGradient(Vector<T>& result, const Vector<T>& b) const;

#ifdef _QT
  bool writeMatrix(QFile& file) const;
  bool readMatrix(QFile& file);
#endif
#ifdef _MATLAB
  void readMatrix(const mxArray* prhs);
  void writeMatrix(mxArray*& prhs) const;
#endif
};

template <class T>
bool Matrix<T>::IsDispInfo = false;

#include "Matrix.cc"

#endif  // MATRIX_H_
