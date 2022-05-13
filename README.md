# CMinMax
A Fast algorithm to find the corners of an N-dimensional Convex Polytope

 - [Paper](https://arxiv.org/abs/2011.14035)
 - [Video](https://www.youtube.com/watch?v=Ug313Nf-S-A) 

This the experimental code for the cMinMax algorithm. 
We did not try higher than 3D dimensions because we did not have such point clouds.
For higher than 3D dimensions we believe that the random version will be more appropriate and it can be extended fairly easy.
A possible avenue for rotating in more than 3 dimentional-space could be implemented as described in Ref #2 “Aguilera,  A.  and  Perez-Aguila,  R.  (2004).General  n-dimensional rotations. 

There are 3 Matlab functions for the cMinMax algorith:
- function [X,Y]=cMinMax(Im_in,Nc)                               % This is for 2D images
- function [ptCloud_Corners]=cMinMax3D(ptCloud, Dphi, Dtheta)    % This is for 3D point clouds
- function [ptCloud_Corners]=cMinMax3Drandom(ptCloud,N)          % This is the random sampling method applie to 3D point clouds

There 3 main files for running the above functions as well some pointclouds and images for testing the algorithm.
The code was developed and tested with Matlab 2017b while the point clouds were generated from MeshLab.

If you use this code in your work please cite the corresponding paper
``` 
@inproceedings{ChamzasCMinMax,
  doi = {10.5220/0010259002290236},
  url = {https://doi.org/10.5220%2F0010259002290236},
  year = 2021,
  publisher = {{SCITEPRESS} - Science and Technology Publications},
  author = {Dimitrios Chamzas and Constantinos Chamzas and Konstantinos Moustakas},
  title = {{cMinMax}: A Fast Algorithm to Find the Corners of an N-dimensional Convex Polytope},
  booktitle = {Proceedings of the 16th International Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications}
}
```

