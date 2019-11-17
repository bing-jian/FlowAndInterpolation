// Author: Ce Liu (c) Dec, 2009; celiu@mit.edu
// Modified By: Deepak Pathak (c) 2016; pathak@berkeley.edu

#ifndef _ImageProcessing_H_
#define _ImageProcessing_H_

#include <typeinfo>
#include "math.h"
//#include "stdio.h"
#include "stdlib.h"
#include "string.h"

//----------------------------------------------------------------------------------
// class to handle basic image processing functions
// this is a collection of template functions. These template functions are
// used in other image classes such as BiImage, IntImage and FImage
//----------------------------------------------------------------------------------

class ImageProcessing {
 public:
  ImageProcessing(void);
  ~ImageProcessing(void);

 public:
  // basic functions
  template <class T>
  static inline T EnforceRange(const T& x, const int& MaxValue) {
    return std::min(std::max(x, 0), MaxValue - 1);
  };

  //---------------------------------------------------------------------------------
  // function to interpolate the image plane
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static inline void BilinearInterpolate(const T1* pImage, int width,
                                         int height, int nChannels, double x,
                                         double y, T2* result);

  template <class T1>
  static inline T1 BilinearInterpolate(const T1* pImage, int width, int height,
                                       double x, double y);

  // the transpose of bilinear interpolation
  template <class T1, class T2>
  static inline void BilinearInterpolate_transpose(const T1* pImage, int width,
                                                   int height, int nChannels,
                                                   double x, double y,
                                                   T2* result);

  template <class T1>
  static inline T1 BilinearInterpolate_transpose(const T1* pImage, int width,
                                                 int height, double x,
                                                 double y);

  template <class T1, class T2>
  static void ResizeImage(const T1* pSrcImage, T2* pDstImage, int SrcWidth,
                          int SrcHeight, int nChannels, double Ratio);

  template <class T1, class T2>
  static void ResizeImage(const T1* pSrcImage, T2* pDstImage, int SrcWidth,
                          int SrcHeight, int nChannels, int DstWidth,
                          int DstHeight);

  //---------------------------------------------------------------------------------
  // functions for 1D filtering
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static void hfiltering(const T1* pSrcImage, T2* pDstImage, int width,
                         int height, int nChannels, const double* pfilter1D,
                         int fsize);

  template <class T1, class T2>
  static void vfiltering(const T1* pSrcImage, T2* pDstImage, int width,
                         int height, int nChannels, const double* pfilter1D,
                         int fsize);

  template <class T1, class T2>
  static void hfiltering_transpose(const T1* pSrcImage, T2* pDstImage,
                                   int width, int height, int nChannels,
                                   const double* pfilter1D, int fsize);

  template <class T1, class T2>
  static void vfiltering_transpose(const T1* pSrcImage, T2* pDstImage,
                                   int width, int height, int nChannels,
                                   const double* pfilter1D, int fsize);

  //---------------------------------------------------------------------------------
  // functions for 2D filtering
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static void filtering(const T1* pSrcImage, T2* pDstImage, int width,
                        int height, int nChannels, const double* pfilter2D,
                        int fsize);

  template <class T1, class T2>
  static void filtering_transpose(const T1* pSrcImage, T2* pDstImage, int width,
                                  int height, int nChannels,
                                  const double* pfilter2D, int fsize);

  template <class T1, class T2>
  static void Laplacian(const T1* pSrcImage, T2* pDstImage, int width,
                        int height, int nChannels);

  //---------------------------------------------------------------------------------
  // functions for sample a patch from the image
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static void getPatch(const T1* pSrcImgae, T2* pPatch, int width, int height,
                       int nChannels, double x, double y, int wsize);

  //---------------------------------------------------------------------------------
  // function to warp image
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static void warpImage(T1* pWarpIm2, const T1* pIm1, const T1* pIm2,
                        const T2* pVx, const T2* pVy, int width, int height,
                        int nChannels);

  template <class T1, class T2>
  static void warpImageFlow(T1* pWarpIm2, const T1* pIm1, const T1* pIm2,
                            const T2* pFlow, int width, int height,
                            int nChannels);

  template <class T1, class T2>
  static void warpImage(T1* pWarpIm2, const T1* pIm2, const T2* pVx,
                        const T2* pVy, int width, int height, int nChannels);

  template <class T1, class T2>
  static void warpImage_transpose(T1* pWarpIm2, const T1* pIm2, const T2* pVx,
                                  const T2* pVy, int width, int height,
                                  int nChannels);

  template <class T1, class T2>
  static void warpImage(T1* pWarpIm2, const T1* pIm2, const T2* flow, int width,
                        int height, int nChannels);

  template <class T1, class T2>
  static void warpImage_transpose(T1* pWarpIm2, const T1* pIm2, const T2* flow,
                                  int width, int height, int nChannels);

  template <class T1, class T2, class T3>
  static void warpImage(T1* pWarpIm2, T3* pMask, const T1* pIm1, const T1* pIm2,
                        const T2* pVx, const T2* pVy, int width, int height,
                        int nChannels);

  //---------------------------------------------------------------------------------
  // function to crop an image
  //---------------------------------------------------------------------------------
  template <class T1, class T2>
  static void cropImage(const T1* pSrcImage, int SrcWidth, int SrcHeight,
                        int nChannels, T2* pDstImage, int Left, int Top,
                        int DstWidth, int DstHeight);
  //---------------------------------------------------------------------------------

  //---------------------------------------------------------------------------------
  // function to generate a 2D Gaussian
  //---------------------------------------------------------------------------------
  template <class T>
  static void generate2DGaussian(T*& pImage, int wsize, double sigma = -1);

  template <class T>
  static void generate1DGaussian(T*& pImage, int wsize, double sigma = -1);
};

#include "ImageProcessing.cc"

#endif  // _ImageProcessing_H_
