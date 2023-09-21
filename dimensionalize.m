function [rho_d T_d p_d V_d A_d x_d dx_d r_d] = dimensionalize(rho, T, p, M, y, Rg, A, x, dx, rho_ref, T_ref, p_ref, r_res_d, L, throat, r_d, As, circShape, crackLength)
rho_d = rho * rho_ref;
T_d = T * T_ref;
p_d = p * p_ref;
V_d = M .* sqrt(y * Rg .* T_d);
%r_d = r.*sqrt(As);
%A_d = pi*r_d.^2;
%fprintf('rd=%f',r(150))
%r_dnew=r_d
%r_d = r.*(As./L);
%fprintf('rdd=%f',r_d(150)) 
%A_d = pi*r_d.^2;
if circShape == false
    A_d = 2.*r_d.*crackLength;
else
    A_d = pi.*r_d.^2;
end
x_d = x .* L;
dx_d = dx .* L;
%r_d=r_dnew
end