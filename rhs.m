%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function contains different f.
% The reader may choose different model second order ellitpitc problems but
% one has to match rhs.m, exactu.m, and gradientu.m function to the same
% second order elliptic problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z=rhs(node)

x=node(:,1);
y=node(:,2);  

%example 1
%z=2.*y.*(1-y)+2*x.*(1-x);  

%example 2
%z=2*pi^2*sin(pi*x).*cos(pi*y);  %ex2 from paper 1

%example 3
%z = (pi^2)*cos(pi*(x+y));

%example 4
%z = -6*x -6*y + 8;

%example 5
z = -2*x.*(x - 1/sqrt(2)) - 2*y.*(y - 1/sqrt(2));

