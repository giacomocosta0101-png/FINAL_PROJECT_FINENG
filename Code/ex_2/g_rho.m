function g = g_rho(theta,p,U,z_1,z_2,Z)

L = [1 0;
    sin(theta) cos(theta)];
R = L*L';
rho = R(1,2);

no_jumps = (U(:,1)<=(1-p(1)) & U(:,2)<=(1-p(2)));

g = mvnpdf([z_1 z_2], zeros(1,2),R)/mvncdf([z_1 z_2],...
            zeros(1, 2), R)*nnz(no_jumps);

invMillsRatio = @(x)normpdf(x)./normcdf(-x);

A_jumps_B_doesnt = (U(:,1)>(1-p(1)) & U(:,2)<=(1-p(2)));

g = g + sum(rho*((z_2-Z(A_jumps_B_doesnt,1))./(1-rho^2)^(3/2)).*...
    invMillsRatio(-(z_2-rho.*Z(A_jumps_B_doesnt,1))./sqrt(1-rho^2)));

B_jumps_A_doesnt = (U(:,2)>(1-p(2)) & U(:,1)<=(1-p(1)));


g = g + sum(rho*((z_1-Z(B_jumps_A_doesnt,2))./(1-rho^2)^(3/2)).*...
    invMillsRatio(-(z_1-rho.*Z(B_jumps_A_doesnt,2))./sqrt(1-rho^2)));

both_jumps = (U(:,1)>(1-p(1)) & U(:,2)>(1-p(2)));

g = g +sum((rho/(1-rho^2)-(rho/(1-rho^2)^2)*...
    (Z(both_jumps,1).^2+Z(both_jumps,2).^2)+...
    ((1+rho^2)/(1-rho^2)^2).*(Z(both_jumps,1).*Z(both_jumps,2))));

end
