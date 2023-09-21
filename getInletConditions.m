% Inlet Conditions
function [p_ref rho_ref T_ref M_ref V_ref] = getInletConditions(M_inlet, p_res, T_res, cv, Rg, y)
% Energy Equation
%T_ref = T_res - 0.5*M_inlet^2/(cv+Rg);
T_ref = T_res/(1+0.5*M_inlet^2*y*Rg/(cv+Rg));
% Adiabatic Equation of State
p_ref = ((T_res^y*p_res^(1-y))/T_ref^y)^(1/(1-y));
% Ideal gas law
rho_ref = p_ref/(Rg*T_ref);
M_ref = M_inlet;
V_ref = M_ref*sqrt(y*Rg*T_ref);
end