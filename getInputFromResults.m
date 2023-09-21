%% getInputFromResults
function [x, r, A, L, dx, Length, n, rho_ref, T_ref, V_ref,  M_ref, p_ref, rho, T, V, M, p] = getInputFromResults(oldResults_xl, circShape)
[T_d rho_d p_d V_d M A_d x_d L L_h cv T_ref Rg rho_ref y As p_ref V_ref M_ref n crackLength] = getOldResults(oldResults_xl);

% generate non-dimensional input %
T = T_d./T_ref;
V = M./ sqrt(T);
rho = rho_d./rho_ref;
p = p_d./p_ref;
if circShape == true
    r_d = sqrt(A_d./pi);
else
    r_d = A_d./(crackLength.*2);
end
r = r_d./(As./L);
A = A_d./As;
x = x_d./L;
dx_d = L/(n-1);
dx = 1/(n-1);
Length = x(end);
end