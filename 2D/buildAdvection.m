function [L] = buildAdvection(r,th,k,a)
%Stiffness matrix for (div_x(au),v) using upwind flux in spherical
%coordinates

%Volume integrals are (a_r*u*dv/dr+a_th*u*1/r*dv/dth *r^2sin(th) )_\W 

%a = a_r \hat{r} + a_th \hat{th} = [a_r,a_th]. We assume a_r,a_th > 0.
a_r = a(1);
a_th = a(2);

num_r = numel(r)-1;
num_th = numel(th)-1;

%Get jacobians/scalings for FE functions on each domain
jac_r = (r(2)-r(1))/2;
jac_th = (th(2)-th(1))/2;

%Create Gaussian quadrature points and weights
[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

%Build Legendre Polynomials of each polynomial degree
% and scale by jac in r direction
% to make orthonormanl in L^2 on a cartesian grid. 
[leg_vals,leg_der_vals,leg_edge_vals,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_r);
leg_der_vals = leg_der_vals/(sqrt(jac_r)*jac_r);
leg_edge_vals = leg_edge_vals/sqrt(jac_r);

%Multiply weights in for speed
test_ref = repmat(w_ref',k+1,1).*leg_vals;
test_der_ref = repmat(w_ref',k+1,1).*leg_der_vals;

%Build necessary volume and edge integrals in the r direction
block_r = cell(11,num_r);
for i=1:num_r
    %Create quadrature points on r cell
    quad_r = quad_ref*(r(i+1)-r(i))/2 + (r(i+1)+r(i))/2;

    %%%Volume integrals
    %(phi_i,phi_j*r)_{\W_r}
    block_r{1,i} = test_ref*(leg_vals.*quad_r)'*jac_r;

    %(phi_i,\grad phi_j*r^2)_{\W_r}
    block_r{2,i} = a_r*test_der_ref*(leg_vals.*quad_r.^2)'*jac_r;
    
    %%%Edge integrals
    
    %%%Lower edge integrals (r = r(i))
    %Sharing volume
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,1)/2;
    jump = -leg_edge_vals(:,1);
    block_r{3,i} = a_r*jump*avg'*r(i)^2;
    
    %<[phi_i],[phi_j]>_e
    jump_i = -leg_edge_vals(:,1);
    jump_j = -leg_edge_vals(:,1);  
    block_r{4,i} = abs(a_r)/2*jump_j*jump_i'*r(i)^2;
       
    %Share edge only
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,1)/2;
    jump = leg_edge_vals(:,2);
    block_r{5,i} =  a_r*jump*avg'*r(i)^2;
    
    %<[phi_i],[phi_j]>_e
    jump_i = -leg_edge_vals(:,1);
    jump_j =  leg_edge_vals(:,2);  
    block_r{6,i} = abs(a_r)/2*jump_j*jump_i'*r(i)^2;   
    
    %%%Upper edge integrals  r = r(i+1)
    %Sharing volume
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,2)/2;
    jump = leg_edge_vals(:,2);
    block_r{7,i} = a_r*jump*avg'*r(i+1)^2;   
    
    %<[phi_i],[phi_j]>_e
    jump_i = leg_edge_vals(:,2);
    jump_j = leg_edge_vals(:,2);  
    block_r{8,i} = abs(a_r)/2*jump_j*jump_i'*r(i+1)^2;

    %<phi_i^-,phi_j^->_e %Used for upper BC
    phi_i = leg_edge_vals(:,2);
    phi_j = leg_edge_vals(:,2);
    block_r{9,i} = a_r*phi_j*phi_i'*r(i+1)^2;

    %Share edge only
    %<{phi_i},[phi_j]>_e
    avg = leg_edge_vals(:,2)/2;
    jump = -leg_edge_vals(:,1);
    block_r{10,i} = a_r*jump*avg'*r(i+1)^2;
    
    %<[phi_i],[phi_j]>_e
    jump_i =  leg_edge_vals(:,2);
    jump_j = -leg_edge_vals(:,1);  
    block_r{11,i} = abs(a_r)/2*jump_j*jump_i'*r(i+1)^2;     

end

%Build Legendre Polynomials of each polynomial degree
% and scale by jac in th direction
% to make orthonormanl in L^2 on a cartesian grid. 
[leg_vals,leg_der_vals,leg_edge_vals,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_th);
leg_der_vals = leg_der_vals/(sqrt(jac_th)*jac_th);
leg_edge_vals = leg_edge_vals/sqrt(jac_th);

%Multiply weights in for speed
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions
test_der_ref = repmat(w_ref',k+1,1).*leg_der_vals;

%Build necessary volume and edge integrals in the theta direction
block_th = cell(2,num_th);
for i=1:num_th

    %Create quadrature points on th cell
    quad_th = quad_ref*(th(i+1)-th(i))/2 + (th(i+1)+th(i))/2;
    
    %Volume integrals
    %(phi_i,phi_j sin(theta))_{\W_th}
    block_th{1,i} = leg_vals*(sin(quad_th).*test_ref)'*jac_th;
    
    %(phi_i,\grad phi_j sin(theta))_{\W_v}
    block_th{2,i} = a_th*leg_der_vals*(sin(quad_th).*test_ref)'*jac_th;

    %%%Edge integrals
    
    %%%Lower edge integrals, th = th(i)
    %Sharing volume
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,1)/2;
    jump = -leg_edge_vals(:,1);
    block_th{3,i} = a_th*jump*avg'*sin(th(i));
    
    %<[phi_i],[phi_j]>_e
    jump_i = -leg_edge_vals(:,1);
    jump_j = -leg_edge_vals(:,1);  
    block_th{4,i} = abs(a_th)/2*jump_j*jump_i'*sin(th(i));
       
    %Share edge only
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,1)/2;
    jump = leg_edge_vals(:,2);
    block_th{5,i} = a_th*jump*avg'*sin(th(i));
    
    %<[phi_i],[phi_j]>_e
    jump_i = -leg_edge_vals(:,1);
    jump_j =  leg_edge_vals(:,2);  
    block_th{6,i} = abs(a_th)/2*jump_j*jump_i'*sin(th(i));   
    
    %%%Upper edge integrals, th = th(i+1)
    %Sharing volume
    %<{phi_i},[phi_j]>_e
    avg  = leg_edge_vals(:,2)/2;
    jump = leg_edge_vals(:,2);
    block_th{7,i} = a_th*jump*avg'*sin(th(i+1));   
    
    %<[phi_i],[phi_j]>_e
    jump_i = leg_edge_vals(:,2);
    jump_j = leg_edge_vals(:,2);  
    block_th{8,i} = abs(a_th)/2*jump_j*jump_i'*sin(th(i+1));

    %<phi_i^-,phi_j^->_e %Used for upper BC
    phi_i = leg_edge_vals(:,2);
    phi_j = leg_edge_vals(:,2);
    block_th{9,i} = a_th*phi_j*phi_i'*sin(th(i+1));

    %Share edge only
    %<{phi_i},[phi_j]>_e
    avg = leg_edge_vals(:,2)/2;
    jump = -leg_edge_vals(:,1);
    block_th{10,i} = a_th*jump*avg'*sin(th(i+1));
    
    %<[phi_i],[phi_j]>_e
    jump_i =  leg_edge_vals(:,2);
    jump_j = -leg_edge_vals(:,1);  
    block_th{11,i} = abs(a_th)/2*jump_j*jump_i'*sin(th(i+1));  
end


%Assemble stiffness matrix
I = zeros((k+1)^2*5*max([num_r num_th]),1);
J = zeros(size(I));
S = zeros(size(I));

count = 1;
for i=1:num_r
    for j=1:num_th
        blockstart = (k+1)^2*((i-1)*num_th+(j-1));
        indices = blockstart+1:blockstart+(k+1)^2;
        %Add volume integrals
        temp = zeros((k+1)^2);
        temp = temp + kron(-block_r{2,i},block_th{1,j});
        temp = temp + kron(-block_r{1,i},block_th{2,j});
        %Add edge integrals.  First shared volume
        if i==1 %We are on bottom.  Do not add bottom edge integrals
            temp_r = block_r{7,i} + block_r{8,i};
        elseif i==num_r %We are on the top.  Add bottom edge integrals 
                        %and modified top edge integral
            temp_r = block_r{3,i} + block_r{4,i} + block_r{9,i};
        else
            temp_r = block_r{3,i} + block_r{4,i} + block_r{7,i} + block_r{8,i};
        end
        temp = temp + kron(temp_r,block_th{1,j});
        if j==1 %We are on bottom.  Do not add bottom edge integrals
            temp_th = block_th{7,j} + block_th{8,j};
        elseif j==num_th %We are on the top.  Add bottom edge integrals 
                        %and modified top edge integral
            temp_th = block_th{3,j} + block_th{4,j} + block_th{9,j};
        else
            temp_th = block_th{3,j} + block_th{4,j} + block_th{7,j} + block_th{8,j};
        end
        temp = temp + kron(block_r{1,i},temp_th);
        J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
        I(count:count+(k+1)^4-1) = reshape(repmat(indices,1,(k+1)^2),[],1);
        S(count:count+(k+1)^4-1) = temp(:);
        count = count + (k+1)^4;

        %Next sharing edge in r direction.
        %First - lower edge
        if i > 1
            indices_mod = indices - (k+1)^2*num_th;
            temp = zeros((k+1)^2);
            temp = temp + kron(block_r{5,i}+block_r{6,i},block_th{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_mod,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        %Upper edge
        if i < num_r
            indices_mod = indices + (k+1)^2*num_th;
            temp = zeros((k+1)^2);
            temp = temp + kron(block_r{10,i}+block_r{11,i},block_th{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_mod,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end

        %Next sharing edge in th direction.
        %First - lower edge
        if j > 1
            indices_mod = indices - (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(block_r{1,i},block_th{5,j}+block_th{6,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_mod,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        %Upper edge
        if j < num_th
            indices_mod = indices + (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(block_r{1,i},block_th{10,j}+block_th{11,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_mod,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end

    end
end
I(count+1:end) = [];
J(count+1:end) = [];
S(count+1:end) = [];

L = sparse(I,J,S);




end