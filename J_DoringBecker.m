%%%%%% VERIFICATION: Nucleation rate %%%%%%
clear all, %close all, clc
% Constants %
m = 2.988*10^(-26); % Molecular mass of water [kg] (Karagiannis)
water_molarmass = 18.01528; % g/mol 
N_A = 6.022140857*10^23; % Avogadro constant, molecules/mol
m = water_molarmass/N_A*10^-3; % kg per molecule
k = 1.38064852*10^(-23); % Boltzmann constant [m^2kgs^-2K^-1]
T_C = 647.096; % Critical Temperature [K] "Revised Release on Surface Tension of Ordinary Water Substance" 
T_C = 647.14; % 647.15? Cite: n.4 (DataPaper_WölkStrey)
% Inputs %%
%T_range = [217.1 222.6 228.2 233.5 238.8 244.1 248.5 253.7 259]; % Viisanen range 
T_range = [259.86 249.75 239.58 229.52 218.94]; % Temperature range of Wölk & STrey
%T_range = [242.391]%%
T = T_range;
x = (T-225)./46.2; % DataPaper
t_r = (T_C - T)./T_C; % DataPaper
%rho = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32;
rho = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32; %[g.cm3]
rho = rho.*10^3;%[kg/m3]
v_m = m./rho; % [1/m3]
% v_m = m_m./rho; % Molecular volume of water [cm^3]
% v_m = v_m*10^(-6);
% v_m = m./rho; % Molecular volume of water [cm^3]
% v_m = v_m*10^(-3);
% Viisanen Strey Reiss 1993
J_data_259 = [9.2e8 3.9e8 1.4e8 8.8e7 3.6e7 1.4e7 5.3e6 1.6e6 9.9e5 1.9e5 7.9e8 6.8e8 2.8e8 1.6e8 6.8e7 2.6e7 1.2e7 3e6 2.8e6 7.5e5];
J_data_2537 = [3.2e8 4.6e8 1.4e8 5.8e7 3.2e7 1.2e7 4.5e6 1.2e6 8.2e8 4.3e8 1.5e8 6.2e7 3.5e7 1.9e7 5e6 2.1e6 7.6e5];  
J_data_2485 = [6.2e8 2.5e8 1e8 5.6e7 2.7e7 1.3e7 6.3e6 2.9e6 6.8e5 6.6e8 4.5e8 1.8e8 1e8 3.5e7 1e7 4.9e6 8.6e8 7.7e8 2.3e8 1.1e8 6.1e7 2.6e7 1.4e7 7.4e6 2.1e6];
J_data_2441 = [7.1e8 3.5e8 1.9e8 8.7e7 3.9e7 2e7 8.3e6 4.5e6 2.1e6 1.1e6 4.2e5 2.8e8 1.4e8 7.5e7 7.7e7 3.7e7 3.4e7 1.5e7 7.4e6 3e6];
J_data_2388 = [1.8e8 1.7e8 9.6e7 4.4e7 2.9e7 1.6e7 7.5e6 2.8e6 1.4e6 4e8 2e8 9.7e7 4.5e7 2.1e7 1.3e7 5.6e6 2.9e6 1.6e6 6.7e5];
J_data_2335 = [1.6e8 1.6e8 8.4e7 3.7e7 2.6e7 6.2e6 1.8e6 2.9e8 9.5e7 2.6e7 1.1e7 7.2e6 3.6e6 1.8e6 9.7e5 4.6e5]; 
J_data_2282 = [1.2e8 8.5e7 4.3e7 1.9e7 6.3e5 1.3e8 1.4e8 6.7e7 4.1e7 2.9e7 1.6e7 8.1e6 5.4e6 1.4e6 1.2e8 8.4e7 3.9e7];
J_data_2226 = [1.3e8 5.7e7 3.5e7 2.5e7 1.2e7 1.3e8 4.4e7 2.3e7 1.1e7 5.6e6 3.9e6 2.2e6 1.4e6 1e6 4.4e5];
J_data_2171 = [5.4e7 3e7 1.9e7 1.1e7 5.2e7 3.3e7 1.5e7 5.5e7 2.1e7 9.5e6 6.3e6 3.1e6];
S_data_259 = [8.73 8.53 8.32 8.2 7.99 7.84 7.68 7.51 7.36 7.22 8.75 8.71 8.48 8.36 8.19 7.97 7.85 7.67 7.55 7.32];
S_data_2537 = [9.59 9.76 9.4 9.17 8.97 8.78 8.58 8.27 9.92 9.61 9.39 9.15 9.09 8.88 8.68 8.47 8.25];
S_data_2485 = [10.94 10.73 10.52 10.25 10.07 9.88 9.66 9.45 9.21 11.1 10.92 10.69 10.49 10.17 9.81 9.63 10.99 11.04 10.68 10.51 10.28 10.05 9.93 9.68 9.42];
S_data_2441 = [12.21 11.94 11.77 11.5 11.26 11.07 10.84 10.65 10.4 10.23 9.98 11.79 11.6 11.41 11.38 11.2 11.16 10.95 10.72 10.49];
S_data_2388 = [13.49 13.45 13.09 12.84 12.69 12.42 12.13 11.83 11.63 13.75 13.38 13.17 12.87 12.57 12.42 12.05 11.87 11.61 11.39];
S_data_2335 = [15.47 15.44 15.16 14.84 14.65 14 13.58 15.77 15.09 14.52 14.16 13.94 13.65 13.41 13.05 12.85];
S_data_2282 = [18.16 17.89 17.53 17.01 15.24 18.24 18.33 17.85 17.5 17.26 16.92 16.53 16.36 15.81 18.23 17.96 17.31];
S_data_2226 = [21.59 20.79 20.49 20.3 19.68 21.51 20.68 20.26 19.65 19.34 19.1 18.61 18.38 18.04 17.57];
S_data_2171 = [24.71 23.96 23.64 23.31 24.56 23.93 23.55 24.47 23.66 22.98 22.76 22.16];
% Nucleation data points, of Wölk and Strey paper: Homogeneous nucleation
% for homogeneous water nucleation.
J_data_25986 = [2.17 2.67 2.9 2.2 2.09 56.7 62.0 63 13.2 13.4 12.4 28.5 28.9 27.5 26.2 5.53 5.3 5.74 ...
    58.5 91.5 91.8 89.9 71.2 13.5 15.2 19.9 20.6 21.5 3.57 31.9 26.1 29.1 4.66 5.29 5 2.12 2.83 2.68 2.62 3.2 ...
    29.0 31.5 35 151 181] .* 10^7; % T = 259.87 K (left-side), Where is the data from? What experiment?
J_data_24975 = [24.2 23.4 26.9 9.08 7.03 7.12 1.56 1.36 1.14 0.983 13.7 13.1 13.8 13.7 3.84 3.99 3.38 3.95 0.692 0.822 0.829 0.823 ...
    34.9 30.6 2.76 8.95 7.55 7.64 6.79 1.23 1.33 1.33 0.43] .* 10^8; % T = 
J_data_23958 = [0.238 .215 .212 .302 .262 .297 .236 .25 .964 2.94 2.42 2.81 2.39 5.86 12.3 12.7 13.8 3.53 3.18 ...
    3.43 4.28 3.71 1.33 1.17 1.21 0.259 0.615 .703 .724 .138 .112 .148] .* 10^8;
J_data_22952 = [2.11 1.91 1.7 .572 .613 .592 .142 .118 1.23 1.36 1.55 .356 .411 .373 .108 .11 1.25 1.11 1.56 4.92 5.51 4.72 5.95 .483 .383 ...
    4.85 3.94 4.41 3.01 2.49 3.05 0.792 .854 .832] .* 10^8; % (left-side)
J_data_21894 = [17.3 17 11.6 22.6 4.31 5.63 5.45 19.9 20.4 15.2 19.2 8.82 10.4 9.65 33.6 35.4 44 34.6 22.6 16.3 14 19 16.9 17.8 5.72 4.63 4.34 1.78 1.27 1.36 2.83 1.43 0.399 ...
    8.61 4.87 6.78] .* 10^7;
S_data_25986 = [7.19 7.18 7.19 7.57 7.58 7.94 7.94 7.93 7.58 7.58 7.61 7.92 7.87 7.9 7.87 7.56 7.52 7.52 7.96 8.04 8.05 8.05 8.02 7.61 ...
    7.78 7.8 7.79 7.8 7.4 8.01 7.95 7.94 7.58 7.56 7.54 7.5 7.49 7.54 7.48 7.53 7.99 7.99 7.98 8.39 8.44]; % T = 259.87 K
S_data_24975 = [10.45 10.56 10.49 10.04 10.08 10 9.52 9.49 9.53 9.53 10.45 10.35 10.33 10.38 9.8 9.85 9.83 9.9 9.31 9.44 9.37 9.41 ...
    10.59 10.55 10.6 10.03 10.06 9.96 10.01 9.47 9.49 9.45 9.24]; % (left-side)
S_data_23958 = [11.34 11.24 11.38 11.43 11.38 11.4 11.37 11.4 11.92 12.28 12.24 12.36 12.26 12.69 13.01 13.09 13.01 12.34 12.26 ...
    12.49 12.52 12.42 11.94 11.93 11.92 11.28 11.82 11.80 11.91 11.12 11.13 11.08];
S_data_22952 = [16.06 16.06 16.28 15.37 15.39 15.31 14.58 14.51 15.97 15.94 16.02 15.19 15.08 14.97 14.37 14.29 ...
    16 15.9 15.92 16.75 16.87 16.63 16.79 15.12 15.07 16.64 16.65 16.67 16.28 16.21 16.29 15.31 15.4 15.44];
S_data_21894 = [21.38 21.64 21.42 21.68 20.27 20.22 20.47 22.38 22.3 22.03 22.19 20.92 21.15 20.94 23.29 23.15 23.35 22.98 21.91 21.83 21.74 ...
    21.83 21.84 22 20.9 20.54 20.54 19.4 19.59 19.58 19.82 19.74 18.53 21.2 21 21.04] ; % FILL IN

% Data T = 259.86 K case %
w_data_25986 = [ 0.02641 .02646 .02623 0.022653 .02587 .0264 .02624 .02651];
p0_data_25986 = {[88.24 88.24 88.24 92.88 92.88],[97.77 97.77 97.77 92.88 92.88 92.88],[97.77 97.77 97.77 97.77 92.88 92.88 92.88], ...
    [97.77 97.77 97.77 97.77 97.77 92.88],[97.77 97.77 97.77 97.77 92.88] [97.77 97.77 97.77 92.88 92.88 92.88],[92.88 92.88 92.88 92.88 92.88],[97.77 97.77 97.77 102.92 102.92]};
p_v_data_25986 = [];
for j = 1:length(w_data_25986)
    p_v_data_25986 = [p_v_data_25986, w_data_25986(j) .* p0_data_25986{j}];
end
p_v_data_25986 = p_v_data_25986*10^3;% kPa -> Pa5    
% Data T = 249.75 K case %
w_data_24975 = [.01677 .01652 .01677];
p0_data_24975 = {[ 91.11 91.11 91.11 86.55 86.55 86.55 82.23 82.23 82.23 82.23],[91.11 91.11 91.11 91.11  86.55 86.55 86.55 86.55 82.23 82.23 82.23 82.23],[91.11 91.11 91.11 86.55 86.55 86.55 86.55 82.23 82.23 82.23 80.17]};
p_v_data_24975 = [];
for j = 1:length(w_data_24975)
    p_v_data_24975 = [p_v_data_24975, w_data_24975(j) .* p0_data_24975{j}];
end
p_v_data_24975 = p_v_data_24975*10^3;% kPa -> Pa5    
% Data T = 239.58 K case %
w_data_23958 = [.006797 .006808 .00665 .006862 .006917 .006827];
p0_data_23958 = {[99.45 99.45 99.45 99.45 99.45], [99.45 99.45 99.45 104.68],[ 110.19 110.19 110.19 110.19 113.02],[ 113.02 113.02 113.02 107.37 107.37],[ 107.37 107.37 107.37 103.08 103.08 103.08 97.92],[ 103.08 103.08 103.08 96.9 96.9 96.9]}; 
p_v_data_23958 = [];
for j = 1:length(w_data_23958)
    p_v_data_23958 = [p_v_data_23958, w_data_23958(j) .* p0_data_23958{j}];
end
p_v_data_23958 = p_v_data_23958*10^3;% kPa -> Pa5    
S{1} = min(S_data_25986):0.01:max(S_data_25986);
S{2} = min(S_data_24975):0.01:max(S_data_24975);
S{3} = min(S_data_23958):0.01:max(S_data_23958);
S{4} = min(S_data_22952):0.01:max(S_data_22952);
S{5} = min(S_data_21894):0.01:max(S_data_21894);
S{1} = min(S_data_25986)-2:0.01:max(S_data_25986)+2;
S{2} = min(S_data_24975)-2:0.01:max(S_data_24975)+2;
S{3} = min(S_data_23958)-2:0.01:max(S_data_23958)+2;
S{4} = min(S_data_22952)-2:0.01:max(S_data_22952)+2;
S{5} = min(S_data_21894)-2:0.01:max(S_data_21894)+2;
%S{1} = [7.95735];%%
% Equations of Wölk and Strey paper
for i=1:length(T_range)
    T = T_range(i);
    x = (T-225)/46.2; % DataPaper
    t_r = (T_C - T)/T_C; % DataPaper
    rho(i) = 0.08*tanh(x) + 0.7415*t_r^0.33 + 0.32;
    rho(i) = rho(i).*10^3;
    %sigma = (235.8*(1-T/T_C)^1.256*(1-0.625*(1-T/T_C)))*10^(-3); % Surface tension [N/m] "Revised Release on Surface Tension of Ordinary Water Substance"
    sigma(i) = (93.6635 + 0.009133*T - 0.000275*T^2)*10^(-3); % cite: n.4,7 (DataPaper_WölkStray)
    
    %p_eq = 610.8.*exp(-5.1421.*log(T./273.15)-6828.77*(1./T-1/273.15)); % Eq. Vapor pressure [pa]
    p_eq(i) = exp(77.34491-7235.42465/T - 8.2*log(T) + 0.0057113*T); % cite: n.66 (DataPaper_WölkStrey)
    %p_eq = 133.3224*10^(19.301142-2892.3693/T - 2.892736*log10(T)-0.0049369728*T+5.606905*10^(-6)*T^2-4.645869*10^(-9)*T^3+3.7874*10^(-12)*T^4); % Viisanen eq.
    %S = {[4:0.1:10],[5:0.1:13],[7:0.1:16],[9:0.1:19],[8:0.1:27]};
    %S = {[linspace(7.18,8.44,length(p_v_data_25986))],[linspace(9.24,10.59,length(p_v_data_24975))],[linspace(11.08,13.09,length(p_v_data_23958))],[14.29:0.1:16.87]};
    %S = p_v ./ p_eq; % Super-saturation
%     if i == 1
%         p_v = p_v_data_25986;
%     elseif i == 2
%         p_v = p_v_data_24975;
%     elseif i == 3
%         p_v = p_v_data_23958;
%     else
%         p_v = p_eq(i) .* S{i};
%     end
    p_v = p_eq(i) .* S{i};
    %p_v = 0.345*[106.51:2.1:127.73]*10^3;
    J_DB = sqrt((2.*sigma(i))./(pi.*m)).*v_m(i).*(p_v ./ (k.*T)).^2.*exp((-16.*pi.*v_m(i).^2.*sigma(i).^3)./(3.*(k.*T).^3 .* log(S{i}).^2)); 
    J_H2O = J_DB .* exp(-27.56 + 6.5*10^3/T)
    figure(21);
    %semilogy(S{i},J_DB/(10^6),'r--')
    %semilogy(S{i},J_DB/(10^6),'r*')%%
    semilogy(S{i},J_H2O/(10^6),'b-','linewidth', 2)
    hold on
    
    %semilogy(S{i},J_H2O/(10^6),'b*')%%
%     legend('J_{DB}','J_{H2O}')
    ylabel('Nucleation [cm^{-3}s^{-1}]')
    xlabel('Super-Saturation')
    
end
grid on
% ylim([10e5 10e8])
semilogy(S_data_25986, J_data_25986, 'k*')
semilogy(S_data_24975, J_data_24975, 'k*')
semilogy(S_data_23958, J_data_23958, 'k*')
semilogy(S_data_22952, J_data_22952, 'k*')
semilogy(S_data_21894, J_data_21894, 'k*')
%title('\gamma_{nuc}')
ylabel('\gamma_{nuc} [1/(cm^3s)]')
xlabel('Super-Saturation [-]')
%legend('T = 260 K','T = 250 K','T = 240 K','T = 230 K','T = 220 K')
ylim([10^6 10^10])
xlim([5 25])
set(gca,'FontSize',14)
%figure(21)
hold on
%%
close all
T = [259.86 249.75 239.58 229.52 218.94]; % Temperature range of Wölk & STrey
S{1} = min(S_data_25986):0.01:max(S_data_25986);
S{2} = min(S_data_24975):0.01:max(S_data_24975);
S{3} = min(S_data_23958):0.01:max(S_data_23958);
S{4} = min(S_data_22952):0.01:max(S_data_22952);
S{5} = min(S_data_21894):0.01:max(S_data_21894);
%T = [259 253.7 248.5 244.1 238.8 233.5 228.2 222.6 217.1];
%T = [230];
%S = {[10:0.01:20]};
x = (T-225)./46.2; % DataPaper
t_r = (T_C - T)./T_C; % DataPaper
rho = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32;
v_m = m./rho; % Molecular volume of water [cm^3]

v_m = v_m*10^(-6);
% S{6} = min(S_data_2335)-2:0.01:max(S_data_2335)+2;
% S{7} = min(S_data_2282)-2:0.01:max(S_data_2282)+2;
% S{8} = min(S_data_2226)-2:0.01:max(S_data_2226)+2;
% S{9} = min(S_data_2171)-2:0.01:max(S_data_2171)+2;

for i=1:length(T)
    x = (T(i)-225)/46.2; % DataPaper
    t_r = (T_C - T(i))/T_C; % DataPaper
    rho(i) = 0.08*tanh(x) + 0.7415*t_r^0.33 + 0.32;
    %sigma = (235.8*(1-T/T_C)^1.256*(1-0.625*(1-T/T_C)))*10^(-3); % Surface tension [N/m] "Revised Release on Surface Tension of Ordinary Water Substance"
    sigma(i) = (93.6635 + 0.009133*T(i) - 0.000275*T(i)^2)*10^(-3); % cite: n.4,7 (DataPaper_WölkStray)
    %p_eq = 610.8.*exp(-5.1421.*log(T./273.15)-6828.77*(1./T-1/273.15)); % Eq. Vapor pressure [pa]
    p_eq(i) = exp(77.34491-7235.42465/T(i) - 8.2*log(T(i)) + 0.0057113*T(i)); % cite: n.66 (DataPaper_WölkStrey)
    %p_eq = 133.3224*10^(19.301142-2892.3693/T - 2.892736*log10(T)-0.0049369728*T+5.606905*10^(-6)*T^2-4.645869*10^(-9)*T^3+3.7874*10^(-12)*T^4); % Viisanen eq.
    %S = {[linspace(7.18,8.44,length(p_v_data_25986))],[linspace(9.24,10.59,length(p_v_data_24975))],[linspace(11.08,13.09,length(p_v_data_23958))],[14.29:0.1:16.87]};
    %S = p_v ./ p_eq; % Super-saturation
%     if i == 1
%         p_v = p_v_data_25986;
%     elseif i == 2
%         p_v = p_v_data_24975;
%     elseif i == 3
%         p_v = p_v_data_23958;
%     else
%         p_v = p_eq(i) .* S{i};
%     end
    p_v = p_eq(i) .* S{i};
    %p_v = 0.345*[106.51:2.1:127.73]*10^3;
    J_DB = sqrt((2.*sigma(i))./(pi.*m)).*v_m(i).*(p_v ./ (k.*T(i))).^2.*exp((-16.*pi.*v_m(i).^2.*sigma(i).^3)./(3.*(k.*T(i)).^3 .* log(S{i}).^2)); 
    J_H2O = J_DB .* exp(-27.56 + 6.5*10^3/T(i));
    figure(1);
    loglog(S{i},J_H2O/(10^6),'r-')
    hold on
    %loglog(S{i},J_H2O/(10^6),'b-')
    legend('J_{DB}')
    ylabel('Nucleation [cm^{-3}s^{-1}]')
    xlabel('Super-Saturation')
    
end


% T = [259 253.7 248.5 244.1 238.8 233.5 228.2 222.6 217.1];
% S{1} = min(S_data_259)-2:0.01:max(S_data_259)+2;
% S{2} = min(S_data_2537)-2:0.01:max(S_data_2537)+2;
% S{3} = min(S_data_2485)-2:0.01:max(S_data_2485)+2;
% S{4} = min(S_data_2441)-2:0.01:max(S_data_2441)+2;
% S{5} = min(S_data_2388)-2:0.01:max(S_data_2388)+2;
% S{6} = min(S_data_2335)-2:0.01:max(S_data_2335)+2;
% S{7} = min(S_data_2282)-2:0.01:max(S_data_2282)+2;
% S{8} = min(S_data_2226)-2:0.01:max(S_data_2226)+2;
% S{9} = min(S_data_2171)-2:0.01:max(S_data_2171)+2;
% S_data_25986 = p_v_data_25986 ./ p_eq(1);
% S_data_24975 = p_v_data_24975 ./ p_eq(2);
% S_data_23958 = p_v_data_23958 ./ p_eq(3);
%S_data_22952 = p_v_data_24975 ./ p_eq(2);
% Figures
%% Viisanen 1993 validation %%
%close all
T = [259 253.7 248.5 244.1 238.8 233.5 228.2 222.6 217.1];
S{1} = min(S_data_259):0.01:max(S_data_259);
S{2} = min(S_data_2537):0.01:max(S_data_2537);
S{3} = min(S_data_2485):0.01:max(S_data_2485);
S{4} = min(S_data_2441):0.01:max(S_data_2441);
S{5} = min(S_data_2388):0.01:max(S_data_2388);
S{6} = min(S_data_2335):0.01:max(S_data_2335);
S{7} = min(S_data_2282):0.01:max(S_data_2282);
S{8} = min(S_data_2226):0.01:max(S_data_2226);
S{9} = min(S_data_2171):0.01:max(S_data_2171);
for i = 1:length(T)
    B = 10.^(-32.025880 - 0.070074889 * T(i) + 0.00040240774 * T(i).^2);
    n = 1.1743362 + 0.011626596 .* T(i) + 0.00047597152 .* T(i).^2;
    g_nuc = (B.*(S{i}-1).^n);
    figure(2)
    loglog(S{i}, g_nuc, '-b')
    hold on
    
end
grid on
% ylim([10e5 10e8])
semilogy(S_data_259, J_data_259, 'k*')
semilogy(S_data_2537, J_data_2537, 'k*')
semilogy(S_data_2485, J_data_2485, 'k*')
semilogy(S_data_2441, J_data_2441, 'k*')
semilogy(S_data_2388, J_data_2388, 'k*')
semilogy(S_data_2335, J_data_2335, 'k*')
semilogy(S_data_2282, J_data_2282, 'k*')
semilogy(S_data_2226, J_data_2226, 'k*')
semilogy(S_data_2171, J_data_2171, 'k*')
% ylim([10^5 10^9])
% xlim([5 30])
ylabel('y_{nuc} [1/(cm^3s)]')
xlabel('Super-Saturation [-]')
%figure(21)
hold on
% semilogy(S_data_25986, J_data_25986, 'r*')
% %semilogy(S_data_259, J_data_259, 'b*')
% semilogy(S_data_24975, J_data_24975, 'r*')
% semilogy(S_data_23958, J_data_23958, 'r*')
% semilogy(S_data_22952, J_data_22952, 'r*')
% semilogy(S_data_21894, J_data_21894, 'r*')
grid on
legend('259.0 K','253.7 K','248.5 K','244.1 K','238.8 K','233.5 K', '228.2 K','222.6 K','217.1 K')
%ylim([10^5 10^10]) 
%xlim([5 30])

% figure(2)
% ax1 = axes; 
% yyaxis left                 % see [1]
% plot(T_range,p_eq)
% ylabel('p_{eq}')
% hold on
% pause(0.1)                  % see [3]
% ax1.XTickMode = 'manual'; 
% ax1.YTickMode = 'manual'; 
% ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
% ax1.XLimMode = 'manual'; 
% grid(ax1,'on')
% ytick = ax1.YTick;  
% yyaxis right                % see [1]
% plot(T_range, sigma)
% ylabel('\sigma')
% ax2 = axes('position', ax1.Position);
% plot(ax2,T_range, rho, 'k')
% ylabel('\rho')
% pause(0.1)                 % see [3]
% ax2.Color = 'none'; 
% grid(ax2, 'on')
% ax2.XLim = ax1.XLim; 
% ax2.XTick = ax1.XTick; 
% ax2.YLimMode = 'manual'; 
% yl = ax2.YLim; 
% ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% % horzontally offset y tick labels
% ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}); 
% %legend('p_{eq}','\sigma','\rho')
% 
% 
% T = [220:2.5:320];
% T_C = 647.15; % Cite: n.4 (DataPaper_WölkStrey)
% x = (T-225)./46.2; % DataPaper
% t_r = (T_C - T)./T_C; % DataPaper
% rho = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32;
% sigma = (93.6635 + 0.009133*T - 0.000275*T.^2).*10^(-3); % cite: n.4,7 (DataPaper_WölkStray)
% p_eq = exp(77.34491-7235.42465./T - 8.2*log(T) + 0.0057113*T); % cite: n.66 (DataPaper_WölkStrey)
% figure(2)
% ax1 = axes; 
% yyaxis left                 % see [1]
% semilogy(T,p_eq)
% ylabel('p_{eq}')
% ylim([10^0 10^5])
% hold on
% pause(0.1)                  % see [3]
% ax1.XTickMode = 'manual'; 
% ax1.YTickMode = 'manual'; 
% ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
% ax1.XLimMode = 'manual'; 
% grid(ax1,'on')
% ytick = ax1.YTick;  
% yyaxis right                % see [1]
% plot(T, sigma)
% ylabel('\sigma')
% ylim([0.065 0.085])
% ax2 = axes('position', ax1.Position);
% plot(ax2,T, rho, 'k')
% ylabel('\rho')
% ylim([0.9 1.15])
% pause(0.1)                 % see [3]
% ax2.Color = 'none'; 
% grid(ax2, 'on')
% ax2.XLim = ax1.XLim; 
% ax2.XTick = ax1.XTick; 
% ax2.YLimMode = 'manual'; 
% yl = ax2.YLim; 
% ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% % horzontally offset y tick labels
% ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}); 
% %legend('p_{eq}','\sigma','\rho')        
