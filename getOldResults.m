function [T_d rho_d p_d V_d M A_d x_d L L_h cv T_ref Rg rho_ref y As p_ref V_ref M_ref n crackLength] = getOldResults(oldResults_xl)
results =load([oldResults_xl, '\runData.txt']);
fullResults = load([oldResults_xl, '\runFullData.txt']);
constants = load([oldResults_xl, '\constantsData.txt']);
% results = xlsread([oldResults_xl, '\runData.xlsx']);
% fullResults = xlsread([oldResults_xl, '\runFullData.xlsx']);
% constants = xlsread([oldResults_xl, '\constantsData.xlsx']);
% Dimensional Results %
T_d = fullResults(:,1);
rho_d = fullResults(:,2);
p_d = fullResults(:,3);
V_d = fullResults(:,4);
M = fullResults(:,5);
f = fullResults(:,6);
R = fullResults(:,7);
S = fullResults(:,8);
gn = fullResults(:,9);
fp = fullResults(:,10);
dRdx = fullResults(:,11);
E_dotL = fullResults(:,12);
E = fullResults(:,13);
A_d = fullResults(:,14);
x_d = fullResults(:,15);
% Transpose %
T_d = transpose(T_d); rho_d = transpose(rho_d); p_d = transpose(p_d); V_d = transpose(V_d); M = transpose(M);
f = transpose(f); R = transpose(R); S = transpose(S); gn = transpose(gn); fp = transpose(fp);
dRdx = transpose(dRdx); E_dotL = transpose(E_dotL); A_d = transpose(A_d); x_d = transpose(x_d);
% constants %
n = length(x_d);
L = constants(1);
L_h = constants(2);
cv = constants(3);
T_ref = constants(4);
Rg = constants(5);
rho_ref = constants(6);
y = constants(7);
As = constants(8);
crackLength = constants(9); 
m_stick = constants(10); 
%r_sub = constants(11); 
%x_d_throat = constants(12);
p_ref = p_d(1);
V_ref = V_d(1);
M_ref = M(1);
end