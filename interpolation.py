from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
import os, time

import numpy as np
import cv2
import pyflow


# Flow Options:
alpha = 0.012
ratio = 0.75
minWidth = 20
nOuterFPIterations = 7
nInnerFPIterations = 1
nSORIterations = 30
colType = 0  # 0 or default:RGB, 1:GRAY (but pass gray image with shape (h,w,1))


def apply_flow(img, u, v):
  flow = np.concatenate((u[..., None], v[..., None]), axis=2)
  h, w = img.shape[0:2]
  x, y = np.meshgrid(np.arange(w), np.arange(h))
  x_new = (x + u).astype(np.float32)
  y_new = (y + v).astype(np.float32)

  x_new[np.where(x_new >= w - 1)] = w - 1
  y_new[np.where(y_new >= h - 1)] = h - 1
  x_new[np.where(x_new < 0)] = 0
  y_new[np.where(y_new < 0)] = 0
  img_new = cv2.remap(img, x_new, y_new, cv2.INTER_LINEAR)
  return img_new


def interpolation(im1, im2, ts):
  ## be careful about the ordering
  uBackward, vBackward, im2WBackward = pyflow.coarse2fine_flow(im2, im1,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  interpolated = []
  for t in ts:
    ut = uBackward * t
    vt = vBackward * t
    interp = apply_flow(im1, ut, vt)
    interpolated.append(interp)
  return interpolated


def test_interp(f1, f2, num_interp=1):
  im1_orig = cv2.imread(f1)
  im2_orig = cv2.imread(f2)
  im1 = im1_orig.astype(float) / 255.
  im2 = im2_orig.astype(float) / 255.
  s = time.time()
  ts = np.linspace(0.0, 1.0, num_interp + 2)
  interpolated = interpolation(im1, im2, ts)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  return [(out*255).astype(int) for out in interpolated], im1_orig, im2_orig


import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input',
            action='store', dest='input_list',
            type=str, nargs='*', help='Input list: -i input1 input2')
    parser.add_argument('-n', '--num_interp', required=True,
            type=int, help='Number of interpolated frames')
    parser.add_argument('-o', '--output',
            type=str, help='Output dir', default='tmp_interp_test')
    args = parser.parse_args()
    assert(len(args.input_list)==2)
    f1 = args.input_list[0]
    f2 = args.input_list[1]
    num_interp = args.num_interp
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    interpolated, im1, im2 = test_interp(f1, f2, num_interp)
    combined = np.concatenate([im1] + interpolated + [im2], axis=1)
    cv2.imwrite('%s/combined.jpg'%output_dir, combined)
    ts = np.linspace(0.0, 1.0, num_interp + 2)
    for k, interp in zip(ts, interpolated):
        cv2.imwrite('%s/interpolated_%03d.jpg'%(output_dir, int(100*k)), interp)
