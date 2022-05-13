clear all; close all;
% generate a bilevel test image 
% with a N-regular polygon mask 
Nc=6; M=1000;
t = (0:1/Nc:1)'*2*pi;t=t+.01;
x = cos(t); y = sin(t);
BW=poly2mask(M*x/4+M/2,M*y/3+M/2,M,M);
Im_in = uint8(255 * BW);
% Im_in=imread('hexagon.png'); Nc=6;
% Im_in=imread('6_imgContour6.png'); Nc=4;
Im_in=imread('convexPolygone01.png'); Im_in=rgb2gray(Im_in); Im_in=255-Im_in; Nc=7;
T1=min(min(Im_in)) ; T2=max(max(Im_in))/2; T=(T1+T2)/2
Im_in(find(Im_in>T))=255; Im_in(find(Im_in<=T))=0; %apply simple threshold

Nc=6
tic
[X,Y]=cMinMax2D(Im_in,Nc);
toc

figure(1);imshow(Im_in); axis on; hold on; 
for i=1:Nc/2
  plot(Y(i,:),X(i,:),'yo','MarkerSize',10,'MarkerFaceColor', 'r' );
end
  
% plot(Y(1,:),X(1,:),'yo','MarkerSize',15,'MarkerFaceColor', 'r' );
% plot(Y(2,:),X(2,:),'yo','MarkerSize',10,'MarkerFaceColor', 'g' );
% plot(Y(3,:),X(3,:),'yo','MarkerSize',5,'MarkerFaceColor', 'b' );