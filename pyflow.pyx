# distutils: language = c++
# distutils: sources = src/Coarse2FineFlowWrapper.cpp
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
import numpy as np
cimport numpy as np
# Author: Deepak Pathak (c) 2016

cdef extern from "src/Coarse2FineFlowWrapper.h":
    void Coarse2FineFlowWrapper(double * vx, double * vy, double * warpI2,
                                const double * Im1, const double * Im2,
                                double alpha, double ratio, int minWidth,
                                int nOuterFPIterations, int nInnerFPIterations,
                                int nSORIterations, int colType,
                                int h, int w, int c);

cdef extern from "src/Coarse2FineFlowWrapper.h":
    void SplatMotionsWrapper(double *vx, double *vy,
                             const double *vxForward, const double *vyForward,
                             const double *vxBackward, const double *vyBackward,
                             const double * Im1, const double * Im2,
                             int colType, int h, int w, int c, double t);

cdef extern from "src/Coarse2FineFlowWrapper.h":
    void ColorTransferWrapper(double *dest, const double *flow,
                             const double *Im1, const double *Im2,
                             const double *forward, const double *backward,
                             int colType, int h, int w, int c, double t);


def coarse2fine_flow(np.ndarray[double, ndim=3, mode="c"] Im1 not None,
                     np.ndarray[double, ndim=3, mode="c"] Im2 not None,
                     double alpha=1, double ratio=0.5, int minWidth=40,
                     int nOuterFPIterations=3, int nInnerFPIterations=1,
                     int nSORIterations=20, int colType=0):
    """
    Input Format:
      double * vx, double * vy, double * warpI2,
      const double * Im1 (range [0,1]), const double * Im2 (range [0,1]),
      double alpha (1), double ratio (0.5), int minWidth (40),
      int nOuterFPIterations (3), int nInnerFPIterations (1),
      int nSORIterations (20),
      int colType (0 or default:RGB, 1:GRAY)
    Images Format: (h,w,c): float64: [0,1]
    """
    cdef int h = Im1.shape[0]
    cdef int w = Im1.shape[1]
    cdef int c = Im1.shape[2]
    cdef np.ndarray[double, ndim=2, mode="c"] vx = \
        np.ascontiguousarray(np.zeros((h, w), dtype=np.float64))
    cdef np.ndarray[double, ndim=2, mode="c"] vy = \
        np.ascontiguousarray(np.zeros((h, w), dtype=np.float64))
    cdef np.ndarray[double, ndim=3, mode="c"] warpI2 = \
        np.ascontiguousarray(np.zeros((h, w, c), dtype=np.float64))
    Im1 = np.ascontiguousarray(Im1)
    Im2 = np.ascontiguousarray(Im2)

    Coarse2FineFlowWrapper(&vx[0, 0], &vy[0, 0], &warpI2[0, 0, 0],
                            &Im1[0, 0, 0], &Im2[0, 0, 0],
                            alpha, ratio, minWidth, nOuterFPIterations,
                            nInnerFPIterations, nSORIterations, colType,
                            h, w, c)
    return vx, vy, warpI2

def splat_motions(np.ndarray[double, ndim=2, mode="c"] vxForward not None,
                  np.ndarray[double, ndim=2, mode="c"] vyForward not None,
                  np.ndarray[double, ndim=2, mode="c"] vxBackward not None,
                  np.ndarray[double, ndim=2, mode="c"] vyBackward not None,
                  np.ndarray[double, ndim=3, mode="c"] Im1 not None,
                  np.ndarray[double, ndim=3, mode="c"] Im2 not None,
                  double t = 0.5, int colType=0):
    """
    Output Format:
      double * vx, double * vy,
      const double * vxForward, const double * vyForward,
      const double * vxBackward, const double * vyBackward,
      const double * Im1 (range [0,1]), const double * Im2 (range [0,1]),
      int colType (0 or default:RGB, 1:GRAY),
      flost t,
    Images Format: (h,w,c): float64: [0,1]
    """
    cdef int h = Im1.shape[0]
    cdef int w = Im1.shape[1]
    cdef int c = Im1.shape[2]
    cdef np.ndarray[double, ndim=2, mode="c"] vx = \
        np.ascontiguousarray(np.zeros((h, w), dtype=np.float64))
    cdef np.ndarray[double, ndim=2, mode="c"] vy = \
        np.ascontiguousarray(np.zeros((h, w), dtype=np.float64))
    vxForward = np.ascontiguousarray(vxForward)
    vyForward = np.ascontiguousarray(vyForward)
    vxBackward = np.ascontiguousarray(vxBackward)
    vyBackward = np.ascontiguousarray(vyBackward)
    Im1 = np.ascontiguousarray(Im1)
    Im2 = np.ascontiguousarray(Im2)

    SplatMotionsWrapper(&vx[0, 0], &vy[0, 0],
                        &vxForward[0, 0], &vyForward[0, 0],
                        &vxBackward[0, 0], &vyBackward[0, 0],
                        &Im1[0, 0, 0], &Im2[0, 0, 0],
                        colType, h, w, c, t)
    return vx, vy


def color_transfer(np.ndarray[double, ndim=3, mode="c"] Im1 not None,
                   np.ndarray[double, ndim=3, mode="c"] Im2 not None,
                   np.ndarray[double, ndim=3, mode="c"] forward not None,
                   np.ndarray[double, ndim=3, mode="c"] backward not None,
                   np.ndarray[double, ndim=3, mode="c"] flow not None,
                   double t = 0.5, int colType=0):



    """
    Output Format:
      double * vx, double * vy,
      const double * vxForward, const double * vyForward,
      const double * vxBackward, const double * vyBackward,
      const double * Im1 (range [0,1]), const double * Im2 (range [0,1]),
      int colType (0 or default:RGB, 1:GRAY),
      flost t,
    Images Format: (h,w,c): float64: [0,1]
    """
    cdef int h = Im1.shape[0]
    cdef int w = Im1.shape[1]
    cdef int c = Im1.shape[2]
    cdef np.ndarray[double, ndim=3, mode="c"] dest = \
        np.ascontiguousarray(np.zeros((h, w, 3), dtype=np.float64))
    flow = np.ascontiguousarray(flow)
    Im1 = np.ascontiguousarray(Im1)
    Im2 = np.ascontiguousarray(Im2)
    forward = np.ascontiguousarray(forward)
    backward = np.ascontiguousarray(backward)

    ColorTransferWrapper(&dest[0, 0, 0], &flow[0, 0, 0],
                         &Im1[0, 0, 0], &Im2[0, 0, 0],
                         &forward[0, 0, 0], &backward[0, 0, 0],
                         colType, h, w, c, t)
    return dest
