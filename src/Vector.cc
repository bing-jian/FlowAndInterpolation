#include "Vector.h"

#include <iostream>
#include <vector>
// template <class T>
// double innerproduct(const Vector<T>& vect1,const Vector<T>& vect2)
//{
//	double result = 0;
//	for(int i = 0;i<vect1.dim();i++)
//		result += vect1[i]*vect2[i];
//	return result;
//}

template <class T>
Vector<T>::Vector(void) {
  nDim = 0;
  pData = NULL;
}

template <class T>
Vector<T>::Vector(int ndim, const T* data) {
  nDim = ndim;
  pData = new T[nDim];
  if (data != NULL)
    memcpy(pData, data, sizeof(T) * nDim);
  else
    memset(pData, 0, sizeof(T) * nDim);
}

template <class T>
Vector<T>::Vector(const Vector& vect) {
  nDim = 0;
  pData = NULL;
  copyData(vect);
}

template <class T>
Vector<T>::~Vector(void) {
  releaseData();
}

template <class T>
void Vector<T>::releaseData() {
  if (pData != NULL) delete[] pData;
  pData = NULL;
  nDim = 0;
}

template <class T>
void Vector<T>::allocate(int ndim) {
  releaseData();
  nDim = ndim;
  if (nDim > 0) {
    pData = new T[nDim];
    reset();
  }
}

template <class T>
void Vector<T>::copyData(const Vector& vect) {
  if (nDim != vect.nDim) {
    releaseData();
    nDim = vect.nDim;
    pData = new T[nDim];
  }
  memcpy(pData, vect.pData, sizeof(T) * nDim);
}

template <class T>
void Vector<T>::dimcheck(const Vector& vect) const {
  if (nDim != vect.nDim)
    std::cout << "The dimensions of the vectors don't match!" << std::endl;
}

template <class T>
void Vector<T>::reset() {
  if (pData != NULL) memset(pData, 0, sizeof(T) * nDim);
}

template <class T>
T Vector<T>::sum() const {
  T total = 0;
  for (int i = 0; i < nDim; i++) total += pData[i];
  return total;
}

template <class T>
double Vector<T>::norm2() const {
  double temp = 0;
  for (int i = 0; i < nDim; i++) temp += pData[i] * pData[i];
  return temp;
}

template <class T>
void Vector<T>::printVector() {
  for (int i = 0; i < nDim; i++) {
    std::cout << pData[i] << ' ';
  }
  std::cout << std::endl;
}

//----------------------------------------------------------------------------------
// operators
//----------------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& vect) {
  copyData(vect);
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& vect) {
  dimcheck(vect);
  for (int i = 0; i < nDim; i++) pData[i] += vect.data()[i];
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator*=(const Vector<T>& vect) {
  dimcheck(vect);
  for (int i = 0; i < nDim; i++) pData[i] *= vect.data()[i];
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& vect) {
  dimcheck(vect);
  for (int i = 0; i < nDim; i++) pData[i] -= vect.data()[i];
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator/=(const Vector<T>& vect) {
  dimcheck(vect);
  for (int i = 0; i < nDim; i++) pData[i] /= vect.data()[i];
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator+=(double val) {
  for (int i = 0; i < nDim; i++) pData[i] += val;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator*=(double val) {
  for (int i = 0; i < nDim; i++) pData[i] *= val;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator-=(double val) {
  for (int i = 0; i < nDim; i++) pData[i] -= val;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator/=(double val) {
  for (int i = 0; i < nDim; i++) pData[i] /= val;
  return *this;
}

template <class T>
const Vector<T> operator+(const Vector<T>& vect1, const Vector<T>& vect2) {
  vect1.dimcheck(vect2);
  Vector<T> result(vect1);
  result += vect2;
  return result;
}

template <class T>
const Vector<T> operator-(const Vector<T>& vect1, const Vector<T>& vect2) {
  vect1.dimcheck(vect2);
  Vector<T> result(vect1);
  result -= vect2;
  return result;
}

template <class T>
const Vector<T> operator*(const Vector<T>& vect1, const Vector<T>& vect2) {
  vect1.dimcheck(vect2);
  Vector<T> result(vect1);
  result *= vect2;
  return result;
}

template <class T>
const Vector<T> operator/(const Vector<T>& vect1, const Vector<T>& vect2) {
  vect1.dimcheck(vect2);
  Vector<T> result(vect1);
  result /= vect2;
  return result;
}

template <class T>
Vector<T> operator+(const Vector<T>& vect, double val) {
  Vector<T> result(vect);
  result += val;
  return result;
}

template <class T>
Vector<T> operator-(const Vector<T>& vect, double val) {
  Vector<T> result(vect);
  result -= val;
  return result;
}

template <class T>
Vector<T> operator*(const Vector<T>& vect, double val) {
  Vector<T> result(vect);
  result *= val;
  return result;
}

template <class T>
Vector<T> operator/(const Vector<T>& vect, double val) {
  Vector<T> result(vect);
  result /= val;
  return result;
}

template <class T>
double innerproduct(const Vector<T>& vect1, const Vector<T>& vect2) {
  vect1.dimcheck(vect2);
  double result = 0;
  for (int i = 0; i < vect1.nDim; i++) {
    result += vect1.pData[i] * vect2.pData[i];
  }
  return result;
}

template <class T>
void Vector<T>::concatenate(const std::vector<Vector<T> >& vect) {
  releaseData();
  nDim = 0;
  for (int i = 0; i < vect.size(); i++) nDim += vect[i].dim();
  allocate(nDim);
  int dim = 0;
  for (int i = 0; i < vect.size(); i++) {
    for (int j = 0; j < vect[i].dim(); j++) pData[dim + j] = vect[i][j];
    dim += vect[i].dim();
  }
}

#ifdef _QT

bool Vector::writeVector(QFile& file) const {
  file.write((char*)&nDim, sizeof(int));
  if (file.write((char*)pData, sizeof(double) * nDim) != sizeof(double) * nDim)
    return false;
  return true;
}

bool Vector::readVector(QFile& file) {
  releaseData();
  file.read((char*)&nDim, sizeof(int));
  if (nDim < 0) return false;
  if (nDim > 0) {
    allocate(nDim);
    if (file.read((char*)pData, sizeof(double) * nDim) != sizeof(double) * nDim)
      return false;
  }
  return true;
}

#endif

#ifdef _MATLAB

template <class T>
void Vector<T>::readVector(const mxArray* prhs) {
  if (pData != NULL) delete[] pData;
  int nElements = mxGetNumberOfDimensions(prhs);
  if (nElements > 2) mexErrMsgTxt("A vector is expected to be loaded!");
  const int* dims = mxGetDimensions(prhs);
  nDim = dims[0] * dims[1];
  pData = new T[nDim];
  double* ptr = (double*)mxGetData(prhs);
  for (int i = 0; i < nDim; i++) pData[i] = ptr[i];
}

template <class T>
void Vector<T>::writeVector(mxArray*& plhs) const {
  int dims[2];
  dims[0] = nDim;
  dims[1] = 1;
  plhs = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double* ptr = (double*)mxGetData(plhs);
  for (int i = 0; i < nDim; i++) ptr[i] = pData[i];
}

#endif
