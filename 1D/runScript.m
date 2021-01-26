%DG Method for u_t + div(au) = 0
%where \W=B_b(0)\setminus B_a(0)\subset\R^3
%The problem is posed in spherical coordinates 
%and the solution u(r,phi,theta) is constant in (phi,theta).

%Endpoints of domain
a = 0;
b = 2*pi;

%Number of elements
N = 16;

%dt and maxtime (initial time is 0)
dt = .001;
T = 50*dt;

%Bilinear form for <\hat{au},[v]>
L = boundTerms(a,b,N);

%Mesh
R = a:(b-a)/N:b;
%Meshes used for plotting the soln
x = (R(2:end)+R(1:end-1))/2;
xplot = sort([R,R(2:end-1),x]);
xfine = a:(b-a)/(5*N):b;

%Mass matrix with respect the L^2 inner product in spherical coordinates
mass = diag((R(2:end).^3-R(1:end-1).^3)/3);

%True soln and its antiderivivative
soln = @(r,t) exp(-3*t) * cos(exp(-t)*r);
soln_anti = @(r,t) ( (exp(-2*t)*r.^2-2).*sin(exp(-t)*r) + 2*exp(-t)*r.*cos(exp(-t)*r) )';

%Initial condition and it's projection.  U is used for plotting
u0 = soln_anti(R(2:end),0) - soln_anti(R(1:end-1),0);
u = mass\u0;
U = repmat(u,1,3)'; U = U(:);

%Start timestepping
figure(31);
t=0;
plot(xplot,U,'o',xfine,soln(xfine,t));
L2 = [];
while t<T
    t = t + dt;
    bc = zeros(N,1);
    %Setting the Dirichlet BC
    bc(1) = a^3*soln(a,t);
    %Solving using BE
    u = (mass-dt*L)\(mass*u + dt*bc);
    %%%Plot solution
    U = repmat(u,1,3)'; U = U(:); %Replicating for plotting
    l2soln = mass\( soln_anti(R(2:end),t) - soln_anti(R(1:end-1),t));
    %plot(x,u',x,soln(x,t));    
    plot(xplot,U,'o',xfine,soln(xfine,t));
    l2 = sqrt((u-l2soln)'*mass*(u-l2soln));
    L2 = [L2 l2];
    title("L2 norm of error is " + num2str(l2,'%e'));
    ylim([-1,1]);
    xlim([a,b]);
    pause(0.05)
    %%%
end