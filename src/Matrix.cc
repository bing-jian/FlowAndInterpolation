#include "Matrix.h"

template <class T>
Matrix<T>::Matrix(void) {
  nRow = nCol = 0;
  pData = NULL;
}

template <class T>
Matrix<T>::Matrix(int nrow, int ncol, double* data) {
  nRow = nrow;
  nCol = ncol;
  pData = new T[nRow * nCol];
  if (data == NULL)
    memset(pData, 0, sizeof(T) * nRow * nCol);
  else
    memcpy(pData, data, sizeof(T) * nRow * nCol);
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& matrix) {
  nRow = nCol = 0;
  pData = NULL;
  copyData(matrix);
}

template <class T>
Matrix<T>::~Matrix(void) {
  releaseData();
}

template <class T>
void Matrix<T>::releaseData() {
  if (pData != NULL) delete pData;
  pData = NULL;
  nRow = nCol = 0;
}

template <class T>
void Matrix<T>::copyData(const Matrix<T>& matrix) {
  if (!dimMatch(matrix)) allocate(matrix);
  memcpy(pData, matrix.pData, sizeof(T) * nRow * nCol);
}

template <class T>
bool Matrix<T>::dimMatch(const Matrix<T>& matrix) const {
  if (nCol == matrix.nCol && nRow == matrix.nRow)
    return true;
  else
    return false;
}

template <class T>
bool Matrix<T>::dimcheck(const Matrix<T>& matrix) const {
  if (!dimMatch(matrix)) {
    cout << "The dimensions of the matrices don't match!" << endl;
    return false;
  }
  return true;
}

template <class T>
void Matrix<T>::reset() {
  if (pData != NULL) memset(pData, 0, sizeof(T) * nRow * nCol);
}

template <class T>
void Matrix<T>::allocate(int nrow, int ncol) {
  releaseData();
  nRow = nrow;
  nCol = ncol;
  if (nRow * nCol > 0) {
    pData = new T[nRow * nCol];
    memset(pData, 0, sizeof(T) * nRow * nCol);
  }
}

template <class T>
void Matrix<T>::allocate(const Matrix<T>& matrix) {
  allocate(matrix.nRow, matrix.nCol);
}

template <class T>
void Matrix<T>::loadData(int _nrow, int _ncol, T* data) {
  if (!matchDimension(_nrow, _ncol)) allocate(_nrow, _ncol);
  memcpy(pData, data, sizeof(T) * nRow * nCol);
}

template <class T>
void Matrix<T>::printMatrix() {
  for (int i = 0; i < nRow; i++) {
    for (int j = 0; j < nCol; j++) cout << pData[i * nCol + j] << " ";
    cout << endl;
  }
}

template <class T>
void Matrix<T>::identity(int ndim) {
  allocate(ndim, ndim);
  reset();
  for (int i = 0; i < ndim; i++) pData[i * ndim + i] = 1;
}

//--------------------------------------------------------------------------------------------------
// functions to check dimensionalities
//--------------------------------------------------------------------------------------------------
template <class T>
bool Matrix<T>::checkDimRight(const Vector<T>& vect) const {
  if (nCol == vect.dim())
    return true;
  else {
    cout << "The matrix and vector don't match in multiplication!" << endl;
    return false;
  }
}

template <class T>
bool Matrix<T>::checkDimRight(const Matrix<T>& matrix) const {
  if (nCol == matrix.nrow())
    return true;
  else {
    cout << "The matrix and matrix don't match in multiplication!" << endl;
    return false;
  }
}

template <class T>
bool Matrix<T>::checkDimLeft(const Vector<T>& vect) const {
  if (nRow == vect.dim())
    return true;
  else {
    cout << "The vector and matrix don't match in multiplication!" << endl;
    return false;
  }
}

template <class T>
bool Matrix<T>::checkDimLeft(const Matrix<T>& matrix) const {
  if (nRow == matrix.ncol())
    return true;
  else {
    cout << "The matrix and matrix don't match in multiplication!" << endl;
    return false;
  }
}

//--------------------------------------------------------------------------------------------------
// functions for numerical computation
//--------------------------------------------------------------------------------------------------
template <class T>
void Matrix<T>::Multiply(Vector<T>& result, const Vector<T>& vect) const {
  checkDimRight(vect);
  if (result.dim() != nRow) result.allocate(nRow);
  for (int i = 0; i < nRow; i++) {
    double temp = 0;
    for (int j = 0; j < nCol; j++) temp += pData[i * nCol + j] * vect.data()[j];
    result.data()[i] = temp;
  }
}

template <class T>
void Matrix<T>::Multiply(Matrix<T>& result, const Matrix<T>& matrix) const {
  checkDimRight(matrix);
  if (!result.matchDimension(nRow, matrix.nCol))
    result.allocate(nRow, matrix.nCol);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < matrix.nCol; j++) {
      double temp = 0;
      for (int k = 0; k < nCol; k++)
        temp += pData[i * nCol + k] * matrix.pData[k * matrix.nCol + j];
      result.pData[i * matrix.nCol + j] = temp;
    }
}

template <class T>
void Matrix<T>::transpose(Matrix<T>& result) const {
  if (!result.matchDimension(nCol, nRow)) result.allocate(nCol, nRow);
  for (int i = 0; i < nCol; i++)
    for (int j = 0; j < nRow; j++)
      result.pData[i * nRow + j] = pData[j * nCol + i];
}

template <class T>
void Matrix<T>::fromVector(const Vector<T>& vect) {
  if (!matchDimension(vect.dim(), 1)) allocate(vect.dim(), 1);
  memcpy(pData, vect.data(), sizeof(double) * vect.dim());
}

template <class T>
double Matrix<T>::norm2() const {
  if (pData == NULL) return 0;
  double temp = 0;
  for (int i = 0; i < nCol * nRow; i++) temp += pData[i] * pData[i];
  return temp;
}

//--------------------------------------------------------------------------------------------------
// operators
//--------------------------------------------------------------------------------------------------
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& matrix) {
  copyData(matrix);
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(double val) {
  for (int i = 0; i < nCol * nRow; i++) pData[i] += val;
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(double val) {
  for (int i = 0; i < nCol * nRow; i++) pData[i] -= val;
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(double val) {
  for (int i = 0; i < nCol * nRow; i++) pData[i] *= val;
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(double val) {
  for (int i = 0; i < nCol * nRow; i++) pData[i] /= val;
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& matrix) {
  dimcheck(matrix);
  for (int i = 0; i < nCol * nRow; i++) pData[i] += matrix.pData[i];
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& matrix) {
  dimcheck(matrix);
  for (int i = 0; i < nCol * nRow; i++) pData[i] -= matrix.pData[i];
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& matrix) {
  dimcheck(matrix);
  for (int i = 0; i < nCol * nRow; i++) pData[i] *= matrix.pData[i];
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& matrix) {
  dimcheck(matrix);
  for (int i = 0; i < nCol * nRow; i++) pData[i] /= matrix.pData[i];
  return *this;
}

template <class T>
Vector<T> operator*(const Matrix<T>& matrix, const Vector<T>& vect) {
  Vector<T> result;
  matrix.Multiply(result, vect);
  return result;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& matrix1, const Matrix<T>& matrix2) {
  Matrix<T> result;
  matrix1.Multiply(result, matrix2);
  return result;
}

//--------------------------------------------------------------------------------------------------
// function for conjugate gradient method
//--------------------------------------------------------------------------------------------------
template <class T>
void Matrix<T>::ConjugateGradient(Vector<T>& result, const Vector<T>& b) const {
  if (nCol != nRow) {
    cout << "Error: when solving Ax=b, A is not square!" << endl;
    return;
  }
  checkDimRight(b);
  if (!result.matchDimension(b)) result.allocate(b);

  Vector<T> r(b), p, q;
  result.reset();

  int nIterations = nRow * 5;
  Vector<T> rou(nIterations);
  for (int k = 0; k < nIterations; k++) {
    rou[k] = r.norm2();
    if (IsDispInfo) cout << rou[k] << endl;

    if (rou[k] < 1E-20) break;
    if (k == 0)
      p = r;
    else {
      double ratio = rou[k] / rou[k - 1];
      p = r + p * ratio;
    }
    Multiply(q, p);
    double alpha = rou[k] / innerproduct(p, q);
    result += p * alpha;
    r -= q * alpha;
  }
}

template <class T>
void Matrix<T>::SolveLinearSystem(Vector<T>& result, const Vector<T>& b) const {
  if (nCol == nRow) {
    ConjugateGradient(result, b);
    return;
  }
  if (nRow < nCol) {
    cout << "Not enough observations for parameter estimation!" << endl;
    return;
  }
  Matrix<T> AT, ATA;
  transpose(AT);
  AT.Multiply(ATA, *this);
  Vector<T> ATb;
  AT.Multiply(ATb, b);
  ATA.ConjugateGradient(result, ATb);
}

#ifdef _QT

template <class T>
bool Matrix<T>::writeMatrix(QFile& file) const {
  file.write((char*)&nRow, sizeof(int));
  file.write((char*)&nCol, sizeof(int));
  if (file.write((char*)pData, sizeof(double) * nRow * nCol) !=
      sizeof(double) * nRow * nCol)
    return false;
  return true;
}

template <class T>
bool Matrix<T>::readMatrix(QFile& file) {
  releaseData();
  file.read((char*)&nRow, sizeof(int));
  file.read((char*)&nCol, sizeof(int));
  if (nRow * nCol > 0) {
    allocate(nRow, nCol);
    if (file.read((char*)pData, sizeof(double) * nRow * nCol) !=
        sizeof(double) * nRow * nCol)
      return false;
  }
  return true;
}
#endif

#ifdef _MATLAB

template <class T>
void Matrix<T>::readMatrix(const mxArray* prhs) {
  if (pData != NULL) delete pData;
  int nElements = mxGetNumberOfDimensions(prhs);
  if (nElements > 2) mexErrMsgTxt("A matrix is expected to be loaded!");
  const int* dims = mxGetDimensions(prhs);
  allocate(dims[0], dims[1]);
  double* data = (double*)mxGetData(prhs);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nCol; j++) pData[i * nCol + j] = data[j * nRow + i];
}

template <class T>
void Matrix<T>::writeMatrix(mxArray*& plhs) const {
  int dims[2];
  dims[0] = nRow;
  dims[1] = nCol;
  plhs = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double* data = (double*)mxGetData(plhs);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nCol; j++) data[j * nRow + i] = pData[i * nCol + j];
}
#endif
