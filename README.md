This project is 2d mpm engine for mac.

<div>
<img src="https://github.com/buaaaaang/MPM_2D/assets/69184903/5089a794-5043-4765-b911-71e9261d676a" width="45%">
<img src="https://github.com/buaaaaang/MPM_2D/assets/69184903/165483c0-9aa8-4587-982e-464b7f2cd470" width="45%">
</div>

The simulation implements the method as detailed in the paper.

> Hu, Yuanming, et al. "A moving least squares material point method with displacement discontinuity and two-way rigid body coupling." ACM Transactions on Graphics (TOG) 37.4 (2018): 1-14.

I used their [simplified implementation](https://github.com/yuanming-hu/taichi_mpm/blob/master/mls-mpm88-explained.cpp) to calculate particles in each frame.

Above that, I made header file for vector and matrix calculation, and applied metal-cpp frameork for graphic engine.

Each frame will be saved as .png. Opening sequence of images by Quicktime Player will give us resulting video.
