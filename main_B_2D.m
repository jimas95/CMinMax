clear all; close all;
% generate a bilevel test image 
% with a N-regular polygon mask 
Nc=5; M=1000;
t = (0:1/Nc:1)'*2*pi;t=t+.01;
x = cos(t); y = sin(t);
BW=poly2mask(M*x/2+M/1,M*y/4+M/2,M,2*M);
Im_in = uint8(255 * BW);
% Im_in=imread('hexagon.png'); Nc=6;
T1=min(min(Im_in)) ; T2=max(max(Im_in))/2;
T=(T1+T2)/2; %apply simple threshold
Im_in(find(Im_in>T))=255;
Im_in(find(Im_in<=T))=0; 
[X,Y]=cMinMax2D(Im_in,Nc);
figure(1);imshow(Im_in); axis on; hold on; 
% plots the detecetd corners
for i=1:Nc/2
  plot(Y(i,:),X(i,:),'yo','MarkerSize',5,...
    'MarkerFaceColor', 'r' );
end
  
plot(Y(1,:),X(1,:),'yo','MarkerSize',15,'MarkerFaceColor', 'r' );
plot(Y(2,:),X(2,:),'yo','MarkerSize',10,'MarkerFaceColor', 'g' );
plot(Y(3,:),X(3,:),'yo','MarkerSize',5,'MarkerFaceColor', 'b' );