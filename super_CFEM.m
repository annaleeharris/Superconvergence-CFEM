%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This super_CFEM.m Matlab program computes convergence rates for H1 norm and L2
%norms.  It also produces surface plots for u_h, numerical FEM
%approximation.  Readers may change different model funtions to test this
%program by choosing differen c (see exactu.m function
%for model second order ellptic problems)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
format short

n=2;
elem=zeros(n^2*2,3);
node=zeros((n+1)^2,2);
k=1;

for nn=1:3
    n=2^nn;
    
   
%------------The Mesh Information------------------------------------------
 for i=1:n+1
    for j=1:n+1
        
    node((i-1)*(n+1)+j,2)=1-(i-1)*1/n;
    node((i-1)*(n+1)+j,1)=0+(j-1)*1/n;
    
    end
    
end

for i=1:n
    for j=1:n
      
     elem((i-1)*2*n+2*(j-1)+1,1)=(i-1)*(n+1)+j;
     elem((i-1)*2*n+2*(j-1)+1,2)=(i)*(n+1)+j;
     elem((i-1)*2*n+2*(j-1)+1,3)=(i-1)*(n+1)+j+1;
     elem((i-1)*2*n+2*(j-1)+2,1)=(i)*(n+1)+j+1;
     elem((i-1)*2*n+2*(j-1)+2,2)=(i-1)*(n+1)+j+1;
     elem((i-1)*2*n+2*(j-1)+2,3)=(i)*(n+1)+j;
     
    end
end

BDNode=zeros(4*n,2);

BDNode(1:n,:)=[2:n+1;1:n]'; % top
BDNode(n+1:2*n,:)=[2*(n+1):n+1:(n+1)*(n+1);n+1:n+1:n*(n+1)]'; % right
BDNode(2*n+1:3*n,:) = [(n+1)*(n+1)-1:-1:n*(n+1)+1;(n+1)*(n+1):-1:n*(n+1)+2]'; % bottom
BDNode(3*n+1:4*n,:) = [ (n-1)*(n+1)+1:-n-1:1;n*(n+1)+1:-n-1:n+2]'; % left


N = size(node,1); NT = size(elem,1); u = zeros(N,1);

%-------- Assemble stiffness matrix -----------------------------------
aa = zeros(NT,3);
bb = zeros(NT,3);
cc = zeros(NT,3);

n1 = elem(:,1);	 n2 = elem(:,2);   n3 = elem(:,3);
x1 = node(n1,1); y1 = node(n1,2);   
x2 = node(n2,1); y2 = node(n2,2);    
x3 = node(n3,1); y3 = node(n3,2);

aa(:,1) = y2-y3;	 bb(:,1) = x3-x2;	
aa(:,2) = y3-y1;	 bb(:,2) = x1-x3;	
aa(:,3) = y1-y2;	 bb(:,3) = x2-x1;	
cc(:,1)=x2.*y3-x3.*y2;
cc(:,2)=x3.*y1-x1.*y3;
cc(:,3)=x1.*y2-x2.*y1;

area = (-x2.*y1+x3.*y1+x1.*y2-x3.*y2-x1.*y3+x2.*y3)/2.0;

A = sparse(N,N);
for i = 1:3
    for j = 1:3
       
        Aij = (aa(:,i).*aa(:,j)+bb(:,i).*bb(:,j))./(4*area);
        A = A + sparse(elem(:,i),elem(:,j), Aij, N,N);
    end
end
%-------- Assembing right hand side by 3-point rule -----------------
mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2;
bt1 = area.*(rhs(mid2)+rhs(mid3))/6;
bt2 = area.*(rhs(mid3)+rhs(mid1))/6;
bt3 = area.*(rhs(mid1)+rhs(mid2))/6;

b = accumarray(elem(:),[bt1;bt2;bt3],[N 1]);


%--------- Handle the boundary condition ----------------------------
for E=1:size(BDNode,1)
    p=BDNode(E,1);
    A(p,:)=0;
    A(p,p)=1;
    b(p)=exactu(node(p,:));
end

uh=A\b;

%----Compute the H1error--------------------------------------------
LocalH1f1error=zeros(NT,1);
l2f1error=0;
val2f1=0;
quad_num=7;
[ quad_w, quad_xy ] = quad_rule ( quad_num );
quad_w=quad_w';
quad_xy=quad_xy';
  
for i=1:3
  n1=elem(i,1);
  n2=elem(i,2);
  n3=elem(i,3);
  
  xx=node(n1,1)*quad_xy(:,1)+node(n2,1)*quad_xy(:,2)+...
     node(n3,1)*(1-quad_xy(:,1)-quad_xy(:,2));
  yy=node(n1,2)*quad_xy(:,1)+node(n2,2)*quad_xy(:,2)+...
     node(n3,2)*(1-quad_xy(:,1)-quad_xy(:,2));
  w(1:quad_num,1) = area(i)* quad_w(1:quad_num);
end

Mid1=(node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
guf1=gradientu(Mid1);

for g1=1:NT
  n1=elem(g1,1);
  n2=elem(g1,2);
  n3=elem(g1,3);
  xx=node(n1,1)*quad_xy(:,1)+node(n2,1)*quad_xy(:,2)+...
     node(n3,1)*(1-quad_xy(:,1)-quad_xy(:,2));
  yy=node(n1,2)*quad_xy(:,1)+node(n2,2)*quad_xy(:,2)+...
     node(n3,2)*(1-quad_xy(:,1)-quad_xy(:,2));
  w(1:quad_num,1) = area(g1)* quad_w(1:quad_num);
  
  val2f1=((aa(g1,1)*xx+bb(g1,1)*yy+cc(g1,1)).*uh(n1)/(2*area(g1))+...
         (aa(g1,2)*xx+bb(g1,2)*yy+cc(g1,2)).*uh(n2)/(2*area(g1))+...
         (aa(g1,3)*xx+bb(g1,3)*yy+cc(g1,3)).*uh(n3)/(2*area(g1)));
    
  trueuf1=exactu([xx,yy]);
  l2f1error=l2f1error+sum(w.*(trueuf1-val2f1).^2);        
end

valf1= [guf1(:,1)-(uh(elem(:,1)).*aa(:,1)+uh(elem(:,2)).*aa(:,2)+uh(elem(:,3)).*aa(:,3))./(2*area),...
        guf1(:,2)-(uh(elem(:,1)).*bb(:,1)+uh(elem(:,2)).*bb(:,2)+uh(elem(:,3)).*bb(:,3))./(2*area)];

LocalH1f1error(:)=(valf1(:,1).*valf1(:,1)+valf1(:,2).*valf1(:,2)).*area;

hh=1/n;
H11error(k)= sum(LocalH1f1error)^(1/2)
L2f1error(k)= sqrt(sum(l2f1error))
h(k)=1/n
k=k+1;

trisurf(elem, node(:,1),node(:,2), uh');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
%axis tight
%shading interp
colorbar
pause(1)
end

%*******Plot the convergence rate****************************************** 
figure
set(gca,'YDir','normal');   
set(gca,'fontsize',20);
loglog(1./h,H11error,'r--s',1./h,h.^(1.3),'g-h','LineWidth',2,'MarkerSize',10)
legend('|u-u_h|_1','O(h)');


figure
set(gca,'YDir','normal');   
set(gca,'fontsize',20);
loglog(1./h,L2f1error,'r--s',1./h,h.^2,'g-h','LineWidth',2,'MarkerSize',10)
legend('||u-u_h||','O(h^2)');

disp('The convergence rate for grad_uh in H1-norm is:')
polyfit(log(h),log(H11error),1)    

disp('The convergence rate for uh in L2-norm is:')
polyfit(log(h),log(L2f1error),1)
