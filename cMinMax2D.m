function [X,Y]=cMinMax(Im_in,Nc)
% find corners in an convex polygon 
%   using cMinMax
% Im_in : The input image with the mask.
%         It is bilevel. Pixels in the
%         mask have high values
% Nc    : The number of expected corners
% X,Y   : contain the detected corners
%         for each rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Authors Dimitrios & Constantinos Chamzas 2020 Matlab 2017b
% This work is described in "cMinMax: A Fast Algorithm to Find the Corners
% of an N-dimensional Convex Polytope"
% by D Chamzas, C Constantinos, K Moustakas
% GRAPP 2021: 16th International Joint Conference on Computer Vision, Imaging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% if image is RGB convert to gray 
  if size(Im_in,3)==3
      Im_in=rgb2gray(Im_in);
  end
  % finds the pixels that are in the mask 
  % and places their coordinates in Xm,Ym
  Im_in_max=max(max(Im_in));
  Im_in_min=min(min(Im_in));
  % T is a threshold
  T=(Im_in_max+Im_in_min)/2;
  [Xm,Ym]=ind2sub(size(Im_in), ...
    find(Im_in>T));

  % finds first four corners in Im_in
  [Xc,Yc] = FndCrnr(Xm,Ym) ;
  X(1,:)=Xc'; Y(1,:)=Yc';
  for k=1:round(Nc/2)-1
    Dtheta=k*pi/Nc; 
    Trot=[cos(Dtheta) -sin(Dtheta);
        sin(Dtheta) cos(Dtheta)];
    % rotate the pixels of the mask
    Indx_rot=Trot*[Xm Ym]';
    [Xcr,Ycr]= FndCrnr(Indx_rot(1,:), ...
      Indx_rot(2,:)); 
    % Xcr,Ycr coordinates of 4 rotated corners
    % Indx_c, coordinates in original img
    % Inverse rotation is the transpose of Trot
    Indx_c=Trot'*[Xcr' Ycr']';
    X(k+1,:)=Indx_c(1,:);
    Y(k+1,:)=Indx_c(2,:);
  end
end

function [Xc, Yc]= FndCrnr(Xm,Ym)
% [Xm,Ym] the coordinates of the mask
% [Xc,Yc] coordinates of the four corners

  % find xmin
  I=find(Xm==min(Xm));
  % if I has more than one points then
  %    for min we select the first and 
  %    for max the last
  Xc(1)=Xm(I(1)); Yc(1)=Ym(I(1));
  % find xmax
  I=find(Xm==max(Xm)); 
  Xc(2)=Xm(I(end)); Yc(2)=Ym(I(end));
  % find ymin
  I=find(Ym==min(Ym)); 
  Xc(3)=Xm(I(1)); Yc(3)=Ym(I(1));
  % find ymax
  I=find(Ym==max(Ym)); 
  Xc(4)=Xm(I(end)); Yc(4)=Ym(I(end));
end 