This package contains the code that implements the following two papers about bounded distortion harmonic mapping

Renjie Chen and Ofir Weber
Bounded Distortion Harmonic Mappings in the Plane 
ACM Transactions on Graphics, 34(4) (SIGGRAPH 2015)

Edward Chien and Renjie Chen* and Ofir Weber
Bounded Distortion Harmonic Shape Interpolation 
ACM Transactions on Graphics, 35(4) (SIGGRAPH 2016)

The app is built with a combination of MATLAB (for core computation) and C++ code (for the UI).
Source C++ with MS VS C++ project files is in the glvu_bdh folder. The C++ code needs the freeglut and anttweakbar libraries to build. A precompiled binary is provided with the package.

Requirements:

- MS Windows (Windows 7/8/10)
- MATLAB (>2014a)
- CVX and Mosek (for deformation)
- A GLSL 3.3 compatible GPU.
- The OpenGL UI (glvu_bdh.exe)

To run the software:

1. Start MATLAB
2. cd to the code folder
3. Call bdh_main.m within MATLAB. This will open the main GUI, and load the data for the giraffe shape

The User Interface (main options):

For deformation, user edits the p2p constraint, i.e. 
    moving the p2p target by dragging any p2p constraint (the red/purple points), 
    adding constraints by left clicking on the shape,
    removing constraints by right clicking the red/purple points,
and then right click on the blank area to run the deformation algorithm.

1. glvu widget

- Continuous deform: by default (when the option is off), the user needs to right click on the blank area to apply the deformation algorithm. When this option is on, the algorithm automatically runs whenever the p2p constraints get edited by the user.

- Iteration: the deformation algorithm is iterative in nature, when this option is on, the app will try run as many iterations as possible to improve the results.

- Clear P2P: remove all the p2p constraints.

- Data Path: set to any subfolder under the "data" folder, and then the relevant data and settings for BDH deformation/interpolation, including the cage of the shape, texture, p2p constraints, and some pre-computed keyframes will be loaded. These data sets are created and bundled with the code for convenience purpose. In order to use the app for other shape, a cage of the shape needs to be extracted and a texture file needs to be provided. (This functionality may be provided later).

2. BDH Deformer widget
* deformation
- numFixedSample, numDenseEvalutionSample, numActiveSetPoolSample: parameters for the active set method used in the deformation algorithm
- Deformer: Sigma1, Sigma2, k: bounds for the conformal and isometric distortions for deformation
- solver: direct Mosek call provides better performance for solving conic optimization than CVX

* interpolation
- Add keyframe: add the current deformation result as a key frame for interpolation
- Set as keyframe: replace the current key frame with the deformation result
- Save all keyframes: save all key frames to disk
- View Keyframe: switch the current key frame
- t: run interpolation with the given time parameter
- InterAlg: choose interpolation algorithm
- numFrame: number of frames for producing continuous interpolation sequence
- Generate sequence: generate a sequence with numFrame frames, that interpolates all the key frames

Note: for interpolation, 2 anchor points need to be picked to determine the global pose of the interpolation result. These anchors are shown as cyan points, which can be set by CTRL + left click on the shape. 




Changelog
08/08/16	 Added system/matlab check for gpu (CUDA) computing suport, and use gpu when possible, otherwise falls back on CPU. The performance may be quite different, especially for interpolation.

10/08/16	 Fix problems with Matlab version < 2015a, where the triangulation class does not have nearestNeighbor method
Improve the performance for interpolation with BDHI/eta and BDHI/nu
Fix interpAlignPose to allow duplicated (single) anchors
