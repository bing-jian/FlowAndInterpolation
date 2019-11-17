// Author: Ce Liu (c) Dec, 2009; celiu@mit.edu
// Modified By: Deepak Pathak (c) 2016; pathak@berkeley.edu

#ifndef VECTOR_H_
#define VECTOR_H_

#include <fstream>
#include <vector>

// TODO(bingjian): Replace it with Eigen.
template <class T>
class Vector {
 protected:
  int nDim;
  T* pData;

 public:
  Vector(void);
  Vector(int ndim, const T* data = NULL);
  Vector(const Vector<T>& vect);
  ~Vector(void);
  void releaseData();
  void allocate(int ndim);
  void allocate(const Vector<T>& vect) { allocate(vect.nDim); };
  void copyData(const Vector<T>& vect);
  void dimcheck(const Vector<T>& vect) const;
  void reset();
  double norm2() const;

  T sum() const;

  void printVector();

  // access the members
  const T* data() const { return (const T*)pData; };
  T* data() { return pData; };
  int dim() const { return nDim; };
  inline bool matchDimension(int _ndim) const { return (nDim == _ndim); }
  inline bool matchDimension(const Vector<T>& vect) const {
    return matchDimension(vect.nDim);
  }

  // operators
  inline T operator[](int index) const { return pData[index]; };
  inline T& operator[](int index) { return *(pData + index); };
  Vector<T>& operator=(const Vector<T>& vect);

  Vector<T>& operator+=(const Vector<T>& vect);
  Vector<T>& operator*=(const Vector<T>& vect);
  Vector<T>& operator-=(const Vector<T>& vect);
  Vector<T>& operator/=(const Vector<T>& vect);

  Vector<T>& operator+=(double val);
  Vector<T>& operator*=(double val);
  Vector<T>& operator-=(double val);
  Vector<T>& operator/=(double val);

  friend double innerproduct(const Vector<T>& vect1, const Vector<T>& vect2) {
    double result = 0;
    for (int i = 0; i < vect1.dim(); i++) {
      result += vect1[i] * vect2[i];
    }
    return result;
  }

  void concatenate(const std::vector<Vector<T> >& vect);

  // friend const Vector<T> concatenate(const vector<Vector<T>>& vect){Vector<T>
  // result; result.concatenate(vect); return result;};
  bool write(std::ofstream& myfile) {
    myfile.write((char*)&nDim, sizeof(int));
    myfile.write((char*)pData, sizeof(T) * nDim);
    return true;
  }
  bool read(std::ifstream& myfile) {
    myfile.read((char*)&nDim, sizeof(int));
    allocate(nDim);
    myfile.read((char*)pData, sizeof(T) * nDim);
    return true;
  }
  T mean(int N = -1) const {
    if (N == -1) {
      N = nDim;
    }
    T result = 0;
    for (int i = 0; i < N; ++i) {
      result += pData[i];
    }
    return result / N;
  }
#ifdef _MATLAB
  void readVector(const mxArray* prhs);
  void writeVector(mxArray*& prhs) const;
#endif
};

#include "Vector.cc"

#endif  // VECTOR_H_
