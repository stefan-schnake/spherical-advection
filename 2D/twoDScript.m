%2D Advection problem in spherical coordinates (r,th);
% u_t+div(au) = 0
% with a = a_r\hat{r}+a_th\hat{th} and a_r,a_th>0. 
% Dirichlet data at \{r=r_0\}\cup\{th=th_0\}.
% Upwinding flux


%%%------- Data -------------------
%Number of cells in each direction
N = 8;

%Polynomial degree (k=0 -> constant, etc...)
%   can be as high as you want.
k = 2;

%dt
dt = 0.05;

%Flow vector (each entry must be non-negative
a = [2,1];

%Solution and source vector
soln = @(r,th,t) 1./sqrt(r).*(exp(th).^(-3)).*cos(-.5*r+t)./sin(th);
source = @(r,th,t) 0*r;

%Upper and lower bounds for r and theta
rr = [.5,2*pi];
thth = [pi/6,5*pi/6];
%%%--------------------------------

%Create uniform mesh in each direction
r = rr(1):(rr(2)-rr(1))/N:rr(2);
th = thth(1):(thth(2)-thth(1))/N:thth(2);

%Build Advection and Mass matrices
L = buildAdvection(r,th,k,a);
M = buildMass(r,th,k);

%Project initial data
u0 = M\(buildNonSeparableSource(r,th,k,@(r,th) soln(r,th,0)));
u = u0;

t = 0;
for i=1:50
    plotVec(r,th,k,u,@(r,th) soln(r,th,t));
    t = t + dt;
    bc = buildDirichletBC(r,th,k,a,@(r,th) soln(r,th,t));
    F = buildNonSeparableSource(r,th,k,@(r,th) source(r,th,t));
    u = (M+dt*L)\(M*u-dt*bc+dt*F);
    pause(0.05)

end