close all ; clear all;
% create solids from meshlab Filters-> Create New Mesh
ptCloud = pcread('dodecahedron.ply');
% ptCloud = pcread('icosahedron.ply');

% % % % % For Debugging  These are the points from Meshlab
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
% ptCloud12=pointCloud(DodecaHedron_corners);
% ptCloud=ptCloud12;
% ptCloud = pcread('icosahedron.ply');
% ptCloud = pcread('cubeh.ply');
% ptCloud = pcread('teapot.ply');
pcshow(ptCloud, 'MarkerSize',6);
% bb=ptCloud.Location;
hold on
Dtheta=8; Dphi=7;
Dtheta=41; Dphi=42;
Dtheta=90/5; Dphi=90/5;
[ptCloud_Corners]=cMinMax3D(ptCloud,Dphi, Dtheta);

aa=ptCloud_Corners.Location;
aa( ~any(aa,2),:)=[] ; % removes zero rows
plot3(aa(:,1),aa(:,2),aa(:,3),'bo','MarkerSize',5,'MarkerFaceColor', 'r' );
%finds number of corners and their centroid
A_Dist=squareform(pdist(aa)); D_max=max(max(A_Dist));
C=(A_Dist<max(max(A_Dist)/5));
D=unique(C,'rows');
[M,N]=size(D);
Corner_Final=[];
% finds the mean of corners that are close
for i=1:length(D(:,1))
  I=find(D(i,:)==1);
  if length(I) == 1
    Corner_Final(i,:)=aa(I,:);
  else    
   Corner_Final(i,:)=mean(aa(I,:));
  end
end
% plot3(DodecaHedron_corners(:,1),DodecaHedron_corners(:,2),DodecaHedron_corners(:,3),'ko','MarkerSize',15) % for debuging

plot3(Corner_Final(:,1),Corner_Final(:,2),Corner_Final(:,3),'bo','MarkerSize',10 );
plot3(aa(:,1),aa(:,2),aa(:,3),'bo','MarkerSize',5,'MarkerFaceColor', 'r' );
fprintf('Size of the point cloude %d\n',length(ptCloud.Location));
fprintf('Total number of Detected Corners %d \nEstimated number of Centroid Corners %d\n',length(aa),size(Corner_Final,1));

