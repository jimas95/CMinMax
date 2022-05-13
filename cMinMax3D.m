function [ptCloud_Corners]=cMinMax3D(ptCloud,Dphi, Dtheta)
% find corners in an convex polyhedron
%   using cMinMax3D
% ptCloud : The input pointcloud
%       
% Dpi, Dtheta  : The rotation angle step in polar coordinates in degrees
% pTCloudCrn : contain the detected corners for each rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Authors Dimitrios & Constantinos Chamzas 2020 Matlab 2017b
% This work is described in "cMinMax: A Fast Algorithm to Find the Corners
% of an N-dimensional Convex Polytope"
% by D Chamzas, C Constantinos, K Moustakas
% GRAPP 2021: 16th International Joint Conference on Computer Vision, Imaging 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Converts degrees to radians for Dphi and Dthera    
  Dtheta=Dtheta*pi/180; Dphi=Dphi*pi/180;  % theta in [-pi/2,pi/2], phi in [-pi,pi]
  N1=round((pi/(2*Dtheta))+1);  N2=round((pi/(Dphi))+1);
  fprintf('Number of Rotations %d (%d,%d) \n',N1*N2, N1, N2);
 pcts_init=ptCloud.Location';
 Corners=[];  N_Corners=zeros(1,3);
 for k=1:N1
   % k
   for m=1:N2
     phi=(k-1)*Dphi; theta=(m-1)*Dtheta;
     Rt_z=[cos(phi) -sin(phi) 0; 
        sin(phi) cos(phi) 0;
        0        0        1];
      Rt_y=[cos(theta) 0 sin(theta) ; 
            0        1        0;
      -sin(theta) 0 cos(theta)  ];
    Rt_zy=Rt_z*Rt_y;
    pcts_rot=Rt_zy*pcts_init;
% pcts_rot=round(pcts_rot,1);
    Xm=pcts_rot(1,:); 
    Ym=pcts_rot(2,:);
    Zm=pcts_rot(3,:);
    [Xc,Yc,Zc] = FndCrnr(Xm,Ym,Zm) ;
    % it applies inverse rotation For Speed you can avoid it and just
    % return the position of the corners
     A = Rt_zy'*[Xc; Yc; Zc] ;
   % plot3(A(1,:),A(2,:),A(3,:),'yo','MarkerSize',5,'MarkerFaceColor', 'r' )
    Corners=[Corners ; A' ];
%     %%% DEBUG STARTS%%%
%     % estimate detected number of corners  (Remove it for speed estmation)
%     D_c=pdist2(A',N_Corners);
%     for i=1:6
%       if min(D_c(i,:))> 0.5
%         N_Corners=[N_Corners ; A(:,i)'];
%         D_c=pdist2(A',N_Corners);
%         N_needed=k*m;   % this is when it finds the last distinct corner
%       end
%     end
% %     %%% DEBUG  END%%%
   end
 end
%  fprintf('Number of rotations %d \n Corners were detected at first %d rotations \n Number of points detected as corners %d \n', k*m, N_needed, (size(N_Corners,1)-1));
 ptCloud_Corners=pointCloud(Corners);   
end

function [Xc, Yc, Zc]= FndCrnr(Xm,Ym,Zm)
% [Xm,Ym,Zm] the coordinates of the mask
% 
  Xc(1:6)=0;Yc(1:6)=0;Zc(1:6)=0;   %Dummy Corners
  % for speed it could move to main for original data
  xm=Xm/max(Xm);ym=Ym/max(Ym);zm=Zm/max(Zm); % normalization
  acc=2; % denoising
  xm=round(Xm,acc); ym=round(Ym,acc); zm=round(Zm,acc);  
  % find xmin
  I=find(xm==min(xm));
  % if I has more than one points then
  % 
  %  We select a corner only if we have ONE max or min
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
%    plot3(Xc,Yc,Zc,'ro','MarkerSize',10,'MarkerFaceColor', 'r' ); hold off
%    [Xc;Yc;Zc]
%    aa=1;

end 