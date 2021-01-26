function [F] = buildDirichletBC(r,th,k,a,u)
%Build Dirichlet BC for (div_x(au),v) using upwind flux in spherical
%coordinates

%a = a_r \hat{r} + a_th \hat{th} = [a_r,a_th]. We assume a_r,a_th > 0.
a_r = a(1);
a_th = a(2);

num_r = numel(r)-1;
num_th = numel(th)-1;

jac_r = (r(2)-r(1))/2;
jac_th = (th(2)-th(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';



F = zeros((k+1)^2*num_r*num_th,1);

[leg_vals,~,leg_edge_vals,~] = buildLegendre(10,k);

for i=1:num_r %Integrals on theta=theta_0 edge first, i.e, j=1;
    blockstart = (k+1)^2*((i-1)*num_th);
    indices = blockstart+1:blockstart+(k+1)^2;
    
    quad_r = quad_ref*(r(i+1)-r(i))/2 + (r(i+1)+r(i))/2;
    
    u_val = u(quad_r,th(1));
    
    temp_integral = (u_val.*w_ref'.*quad_r)*(leg_vals/sqrt(jac_r))'*jac_r;
    
    F(indices) = F(indices) + ...
        kron(temp_integral',-leg_edge_vals(:,1)/sqrt(jac_th)*sin(th(1))*a_th);
end

for j=1:num_r %Integrals on r=r_0 edge first, i.e, i=1;
    blockstart = (k+1)^2*(j-1);
    indices = blockstart+1:blockstart+(k+1)^2;
    
    quad_th = quad_ref*(th(j+1)-th(j))/2 + (th(j+1)+th(j))/2;

    u_val = u(r(1),quad_th);
    
    temp_integral = (u_val.*w_ref'.*sin(quad_th))*(leg_vals/sqrt(jac_th))'*jac_th;
    
    F(indices) = F(indices) + ...
        kron(-leg_edge_vals(:,1)/sqrt(jac_r)*(r(1)^2)*a_r,temp_integral');
end


end