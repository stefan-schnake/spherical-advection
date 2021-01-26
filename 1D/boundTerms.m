function L = boundTerms(a,b,N)
%Constant DG scheme for 
% u_t + div(au) = 0
% where a = r and an upwind flux is used
% and the problem is posed in spherical coordinates


%No interior integration because constants - only flux values
R = a:(b-a)/N:b;

L = zeros(numel(R)-1);

%Because a=r>0 and we are using an upwind flux, only the right edge on each
%cell is taken into consideration.  
L(1,1) = -R(2)^3;
L(2,1) =  R(2)^3;
for i=2:N-1
    L(i,i) = -R(i+1)^3;
    L(i+1,i) = R(i+1)^3;
end
L(N,N) = -b^3;

end

