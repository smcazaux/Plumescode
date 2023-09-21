% SealingTime
%function days = sealingTime()
%[wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, run_name, k, nt, predictor, m_stick, T_w, c_d, r_sub) % [kg/s]

[wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, 'xxx', 1, 2000, 1, 3, T_w, c_d, r_sub); % [kg/s]
%%
%E_max = max(E); % Maximum accretion rate per unit area [kg/m2s]
%x = find(E_max == E);
%D = r_d(x)*2
% E_max = 0.06;
% D = 2;
E_max = 0.0117;
D = 0.195;
m0 = 2.992*10^(-26); % weight of 1 water molecule [kg]
n_molecules_m2s = E_max/m0; % nr. molecules accreted per m2 per s
A = 10^-10; % Angstrom
n_molecules_m2 = 1/(3*A)^2; % nr. of molecules per m2 (layer)
layers_s = n_molecules_m2s/n_molecules_m2; % layers accreted per second
layer = 3*A; % layer thickness [m]
m_s = layers_s*layer;
time = D/m_s;
days = time/(24*3600)