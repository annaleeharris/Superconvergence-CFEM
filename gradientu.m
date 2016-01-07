%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function contains gradients of u.
% The reader may choose different model second order ellitpitc problems but
% one has to match rhs.m, exactu.m, and gradientu.m function to the same
% second order elliptic problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gu=gradientu(node)

x=node(:,1);
y=node(:,2);

% example 1
%gux=(1-x).*y.*(1-y)-x.*y.*(1-y);  
%guy=x.*(1-x).*(1-y)-x.*(1-x).*y;

%example 2
%gux = pi*cos(pi*y).*cos(pi*x);  
%guy = -pi*sin(pi*x).*sin(pi*y);

%example 3
%gux = -.5*pi*sin(pi*(x+y));  
%guy = -.5*pi*sin(pi*(x+y));

%example 4
%gux = (1-x).*(1-3*x);
%guy = (1-y).*(1-3*y);

%example 5
gux = y.*(y - 1/sqrt(2)).*(2*x - 1/sqrt(2));
guy = x.*(x - 1/sqrt(2)).*(2*y - 1/sqrt(2));

gu=[gux,guy];
