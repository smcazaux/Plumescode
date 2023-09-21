% Reynolds, Nusselt, Prandtl, Rayleig numbers
function [h Re Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d)
mu = 0.925e-5*(T_d./300).^1.1; % dynamic viscosity (ingersoll)
cp = 1850; % (ingersoll)
k = 3; % thermal conductivity (ingersoll)
Re = rho_d.*V_d.*2.*r_d./mu; % Reynolds
Pr = cp.*mu./k; % Prandlt
h = 0.023.*Re.^0.8.*Pr.^0.4.*k./(2.*r_d); % convective heat transfer coefficient
Nu = h.*2.*r_d./k; % Nusselt
end