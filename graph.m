%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab computer program plots different model second order elliptic
% problems so one can verify a numerical approximation plot is
% visually correct. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=[0:.05:1]; 
y=[0:.05:1];


[xx,yy] = meshgrid(x,y);

%example 1
%zz = xx.*(1-xx).*yy.*(1-yy);       

%example 2
%z1 = sin(pi*xx).*cos(pi*yy);       

%example 3
%z1 = .5*cos(pi*(xx+yy));            

%example 4
%z1 = xx.*(1-xx).^2 + yy.*(1-yy).^2;

%example 5
%z1 = xx.*yy.*(xx - 1/sqrt(2)).*(yy-1/sqrt(2));  %wg-fem ex1


figure
surf(xx,yy,z1)
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
colorbar
