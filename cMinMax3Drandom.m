function [ptCloud_Corners]=cMinMax3Drandom(ptCloud,N)
% find corners in an convex polyhedron
%   using cMinMax3D and Random Sampling
% ptCloud : The input pointcloud
% N is the number of random rotations
% Phi and Theta are the random values for Theta and Phi
%       
% Nc         : The number of expected corners
% pTCloudCrn : contain the detected corners
%             for each rotation

% Generates N random pairs for rotation angle phi and theta  
% Phi is in [-pi, pi] and Theta in [-pi/2,pi/2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Authors Dimitrios & Constantinos Chamzas 2020 Matlab 2017b
% This work is described in "cMinMax: A Fast Algorithm to Find the Corners
% of an N-dimensional Convex Polytope"
% by D Chamzas, C Constantinos, K Moustakas
% GRAPP 2021: 16th International Joint Conference on Computer Vision, Imaging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % % % % % % For Debugging  These are the points from Meshlab
% % DodecaHedron_corners=[
% %  0.0000   -1.61803   -4.23607;  -2.61803    2.61803    2.61803;
% %  2.61803   -2.61803   -2.61803;   0.0000    1.61803    4.23607;
% % -2.61803   -2.61803   -2.61803;   2.61803    2.61803    2.61803;
% %  1.61803   -4.23607    0.00000;  1.61803    4.23607   0.00000;
% %  -2.61803   -2.61803    2.61803;  2.61803    2.61803   -2.61803;
% %  0.00000    1.61803   -4.23607; -1.5471    4.23607    0.00000;
% %  -1.5708   -4.23607   0.00000;  -0.00000   -1.6246    4.23607;
% %  2.61803   -2.61803    2.61803; -2.61803    2.61803   -2.61803;
% %  -4.23607    0.0334    1.5692; 4.23607    0.00000   -1.61803;
% %  -4.23607   0.00000   -1.61803; 4.23607   0.00000    1.61803 ]; 

x=randn([1 N]);y=randn([1 N]);z=randn([1 N]);
[Phi,Theta,r]=cart2sph(x,y,z);
% generates N random integers 1 to 6 for random of rotation around the axis
R_order=randi(6,1,N);

% Obtains the point locations
 pcts_init=ptCloud.Location';
 Corners=[];  N_Corners=zeros(1,3);
 for k=1:N
   % generates a random rotation matrix
   
   Rt= RotationMatrix(Phi(k),Theta(k),R_order(k));
    
    pcts_rot=Rt*pcts_init;
% pcts_rot=round(pcts_rot,1);
    Xm=pcts_rot(1,:); 
    Ym=pcts_rot(2,:);
    Zm=pcts_rot(3,:);
    [Xc,Yc,Zc] = FndCrnr(Xm,Ym,Zm) ;
    A = Rt'*[Xc; Yc; Zc] ;   %Rt' transpose is the inverse of Rt
   % plot3(A(1,:),A(2,:),A(3,:),'yo','MarkerSize',5,'MarkerFaceColor', 'r' )
    Corners=[Corners ; A' ];
    % estimate detected number of corners  (Remove it for speed estmation)
    D_c=pdist2(A',N_Corners);
    for i=1:6
      if min(D_c(i,:))> 0.5
        N_Corners=[N_Corners ; A(:,i)'];
        D_c=pdist2(A',N_Corners);
        N_needed=k;
       end
    end
% %  % for debugging Finds how many of the Dodocaedron corners have been already detected
% %      if (length(find(min(pdist2(Corners, DodecaHedron_corners))<0.4))==20)
% %        break
% %      end
 end
 fprintf('Number of rotations %d \n Corners were detected at first %d rotations \n Number of points detected as corners %d \n', k, N_needed, (size(N_Corners,1)-1));
 ptCloud_Corners=pointCloud(Corners);   
end

function [Xc, Yc, Zc]= FndCrnr(Xm,Ym,Zm)
% [Xm,Ym,Zm] the coordinates of the mask
% [Xc,Yc,Zc] coordinates of the six corners
  Xc(1:6)=0;Yc(1:6)=0;Zc(1:6)=0;
  % for speed it could move to main for original data
  xm=Xm/max(Xm);ym=Ym/max(Ym);zm=Zm/max(Zm); % normalization
  acc=1; % denoising
  xm=round(Xm,acc); ym=round(Ym,acc); zm=round(Zm,acc);  
  % find xmin
  I=find(xm==min(xm));
  % if I has more than one points then
  %    for min we select the first and 
  %    for max the last
   if length(I)==1
     Xc(1)=Xm(I(1)); Yc(1)=Ym(I(1)); Zc(1)=Zm(I(1));
   end
  % find xmax
  I=find(xm==max(xm)); 
  
   if length(I)==1
     Xc(2)=Xm(I(end)); Yc(2)=Ym(I(end)); Zc(2)=Zm(I(end));
   end
   
  % find ymin
  I=find(ym==min(ym)); 
  if length(I)==1
    Xc(3)=Xm(I(1)); Yc(3)=Ym(I(1)); Zc(3)=Zm(I(1));
  end
  % find ymax
  I=find(ym==max(ym)); 
   if length(I)==1
     Xc(4)=Xm(I(end)); Yc(4)=Ym(I(end)); Zc(4)=Zm(I(end));
   end
  
  % find zmin
  I=find(zm==min(zm)); 
   if length(I)==1
     Xc(5)=Xm(I(1)); Yc(5)=Ym(I(1)); Zc(5)=Zm(I(1));
    end
  % find zmax
  I=find(zm==max(zm)); 
   if length(I)==1
     Xc(6)=Xm(I(end)); Yc(6)=Ym(I(end)); Zc(6)=Zm(I(end));
   end
%    % debug
%    figure(2); plot3(Xm,Ym,Zm,'bo','MarkerSize',6); hold on
%    plot3(Xc,Yc,Zc,'ro','MarkerSize',4,'MarkerFaceColor', 'r' ); hold off
%    aa=1;
end 

function Rt=RotationMatrix(ph,th,Nr)

  switch Nr
            % the first rotation is with phi (-p, p)
            case 1    % xy
              theta=ph; omega=th;  
              Rx=[ 1        0          0;
                  0    cos(omega) -sin(omega);
                  0    sin(omega) cos(omega)];
              Ry=[cos(theta) 0     sin(theta);
                 0          1        0;
                -sin(theta) 0     cos(theta)];
              Rt=Rx*Ry;
                
            case 2  % yx
              omega=ph; theta=th;  
              Rx=[ 1        0          0;
                  0    cos(omega) -sin(omega);
                  0    sin(omega) cos(omega)];
              Ry=[cos(theta) 0     sin(theta);
                 0          1        0;
                -sin(theta) 0     cos(theta)];
              Rt=Ry*Rx;
            case 3  % xz
             phi=ph; omega=th;   
              Rx=[ 1        0          0;
                  0    cos(omega) -sin(omega);
                  0    sin(omega) cos(omega)];
              Rz=[cos(phi) -sin(phi) 0;
                sin(phi)   cos(phi) 0;
                0          0        1];
              Rt=Rx*Rz;
            case 4  % zx
              omega=ph; phi=th;  
              Rx=[ 1        0          0;
                  0    cos(omega) -sin(omega);
                  0    sin(omega) cos(omega)];
              Rz=[cos(phi) -sin(phi) 0;
                sin(phi)   cos(phi) 0;
                0          0        1];
              Rt=Rz*Rx;
            case 5  % yz
              phi=ph; theta=th;   
              Ry=[cos(theta) 0     sin(theta);
                 0          1        0;
                -sin(theta) 0     cos(theta)];
              Rz=[cos(phi) -sin(phi) 0;
                sin(phi)   cos(phi) 0;
                0          0        1];
              Rt=Ry*Rz;
              
            case 6  % zy
              theta=ph; phi=th;  
              Ry=[cos(theta) 0     sin(theta);
                 0          1        0;
                -sin(theta) 0     cos(theta)];
              Rz=[cos(phi) -sin(phi) 0;
                sin(phi)   cos(phi) 0;
                0          0        1];
              Rt=Rz*Ry;
  end
end

