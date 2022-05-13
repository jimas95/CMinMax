This is an experimental code that  we developed for cMinMax. If we find time we hope to optimize it for speed and functionality but currently is a code just for proof of concept.
We did not try higher than 3D dimensions because we did not have such point clouds.
For higher than 3D dimensions we believe that the random version will be more appropriate and it can be extended fairly easy.
Maybe Rotating in higher than 3D spaces could be implemented as described in Ref #2 “Aguilera,  A.  and  Perez-Aguila,  R.  (2004).General  n-dimensional rotations, but we did not implemented.
The new version (v2) uploaded in arxiv https://arxiv.org/pdf/2011.14035v2.pdf  is more detailed and also you can watch our conference oral presentation in https://www.youtube.com/watch?v=Ug313Nf-S-A .

There are 3 matlab functions for the cMinMax algorithm.
function [X,Y]=cMinMax(Im_in,Nc)                               % This is for 2D images
function [ptCloud_Corners]=cMinMax3D(ptCloud, Dphi, Dtheta)    % This is for 3D point clouds
function [ptCloud_Corners]=cMinMax3Drandom(ptCloud,N)          % This is the random sampling method applie to 3D point clouds

There are also 3 Main files for running the above functions as well some pointclouds and images for testing the algorithm

We did use Matlab 2017b while the point clouds were generated from MeshLab

Please reference our publication 
Chamzas, D., Chamzas, C. and Moustakas, K. "cMinMax: A Fast Algorithm to Find the Corners of an N-dimensional Convex Polytope", 
Proceedings of the 16th International Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications (VISIGRAPP 2021) - Volume 1: GRAPP, pages 229-236 , DOI: 10.5220/0010259002290236

Dimitris Chamzas   chamzas95@gmail.com
Constantinos Chamzas chamzask@gmail.com














