function [L] = buildMass(r,th,k)
%Stiffness matrix in spherical coordinates, jacobian = r^2sin(th)

num_r = numel(r)-1;
num_th = numel(th)-1;

jac_r = (r(2)-r(1))/2;
jac_th = (th(2)-th(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals_r = leg_vals/sqrt(jac_r);
test_ref = repmat(w_ref',k+1,1).*leg_vals_r; %Weights included with the test functions

block_r = cell(1,num_r);
for i=1:num_r
    quad_r = quad_ref*(r(i+1)-r(i))/2 + (r(i+1)+r(i))/2;

    %%%Volume integrals
    %(phi_i,phi_j*r^2)_{\W_r}
    block_r{1,i} = test_ref*(leg_vals_r.*(quad_r.^2))'*jac_r;

end

leg_vals_th = leg_vals/sqrt(jac_th);
test_ref = repmat(w_ref',k+1,1).*leg_vals_th; %Weights included with the test functions

block_th = cell(2,num_th);
for i=1:num_th

    %Create quadrature points on x cell
    quad_th = quad_ref*(th(i+1)-th(i))/2 + (th(i+1)+th(i))/2;
    
    %Volume integrals
    %(phi_i,phi_j sin(theta))_{\W_th}
    block_th{1,i} = leg_vals_th*(sin(quad_th).*test_ref)'*jac_th;

    
end

%%Assemble stiffness matrix
% L = zeros((k+1)^2*num_x*num_v);
% for i=1:num_x
%     for j=1:num_v
%         blockstart = (k+1)^2*((i-1)*num_v+(j-1));
%         indices = blockstart+1:blockstart+(k+1)^2;
%         L(indices,indices) = L(indices,indices) + kron(-block_x{1,i},block_v{1,j});
%         L(indices,indices) = L(indices,indices) + kron( block_x{2,i},block_v{1,j});
%         L(indices,indices) = L(indices,indices) + kron( block_x{3,i},block_v{2,j});
%         %Interations with x-block "below" element
%         indices_dn = indices - (k+1)^2*num_v;
%         if indices_dn(1) <= 0
%             indices_dn = indices_dn + (k+1)^2*num_x*num_v;
%         end        
%         L(indices_dn,indices) = L(indices_dn,indices) + kron(block_x{4,i},block_v{1,j});
%         L(indices_dn,indices) = L(indices_dn,indices) + kron(block_x{5,i},block_v{2,j});
%         %Interations with x-block "above" element
%         indices_up = indices + (k+1)^2*num_v;
%         if indices_up(1) >= (k+1)^2*num_x*num_v
%             indices_up = indices_up - (k+1)^2*num_x*num_v;
%         end
%         L(indices_up,indices) = L(indices_up,indices) + kron(block_x{6,i},block_v{1,j});
%         L(indices_up,indices) = L(indices_up,indices) + kron(block_x{7,i},block_v{2,j});
%     end
% end

%sparsify
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
        temp = temp + kron(block_r{1,i},block_th{1,j});
        J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
        I(count:count+(k+1)^4-1) = reshape(repmat(indices,1,(k+1)^2),[],1);
        S(count:count+(k+1)^4-1) = temp(:);
        count = count + (k+1)^4;
    end
end
I(count+1:end) = [];
J(count+1:end) = [];
S(count+1:end) = [];

L = sparse(I,J,S);




end