%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the exact solution of u.
% The reader may choose different model second order ellitpitc problems but
% one has to match rhs.m, exactu.m, and gradientu.m function to the same
% second order elliptic problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=exactu(node)

x=node(:,1);
y=node(:,2);

%example 1
%u=x.*(1-x).*y.*(1-y);  

%example 2
%u = sin(pi*x).*cos(pi*y);  %ex 2 from paper 1

%example 3
%u = .5*cos(pi*(x+y));

%example 4
%u = x.*(1-x).^2 + y.*(1-y).^2;

%exaample 5
u = x.*(x - 1/sqrt(2)).*y.*(y - 1/sqrt(2));




