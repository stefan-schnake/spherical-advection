function [] = plotVec(x,v,k,u,utrue)
%Plots u

num_x = numel(x)-1;
num_v = numel(v)-1;
poi = 4;

[quad_ref, ~]  = lgwt(poi,-1,1);
quad_ref = quad_ref';
quad_ref = fliplr(quad_ref);

[leg_vals,~,~,~] = buildLegendre(poi,k);
leg_vals = fliplr(leg_vals);

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

X = zeros(poi*num_x,1);
V = zeros(poi*num_v,1);
for i=1:num_x
    X((i-1)*poi+1:i*poi) = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
end
for i=1:num_v
    V((i-1)*poi+1:i*poi) = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
end

Z = zeros(numel(X),numel(V));
count = 1;
for i=1:num_x
    for j=1:num_v
        coeff = u(count:count+(k+1)^2-1);
        temp = zeros(poi);
        for ii=1:(k+1)
            for jj=1:(k+1)
                temp = temp + coeff((ii-1)*(k+1)+jj)*(leg_vals(ii,:)/sqrt(jac_x))'*(leg_vals(jj,:)/sqrt(jac_v));
            end
        end
        Z((i-1)*poi+1:i*poi,(j-1)*poi+1:j*poi) = temp;
        count = count + (k+1)^2;
    end
end
Z = Z';

if nargin <= 4
    figure(5)
    surf(X,V,Z);
    xlabel('r');
    ylabel('th');
    title('Plot of numerical soln')
else
    subplot(2,2,1)
    surf(X,V,Z);
    xlabel('r');
    ylabel('th');
    title('Plot of numerical soln')
    [XX,VV] = meshgrid(X,V);
    subplot(2,2,2)
    surf(XX,VV,utrue(XX,VV));
    xlabel('r');
    ylabel('th');
    title('Plot of true soln');
    subplot(2,2,3)
    surf(XX,VV,Z-utrue(XX,VV));
    xlabel('r');
    ylabel('th');
    error = Z-utrue(XX,VV);
    title("Plot of error");
    subplot(2,2,4)
    h = findobj(gca,'Type','line');
    T = get(h,'Xdata');
    Linf = get(h,'Ydata');
    if isempty(T)
        T = 1;
        Linf = norm(reshape(error,[],1),inf);
    else
        T = [T,T(end)+1];
        Linf = [Linf norm(reshape(error,[],1),inf)];
    end
    plot(T,Linf);
    title("L^{inf} Norm");
end
%surf(X,V,Z-sin(2*pi*X)*(V.^2-1)');
%figure()
%Y = sort(X);
%surf(Y,Y,sin(2*pi*Y)*(Y.^2-1)');

end

