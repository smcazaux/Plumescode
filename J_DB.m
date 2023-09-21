%%%%%% VERIFICATION: Nucleation rate %%%%%%
clear all, close all, clc
display('run')
% Constants %
m = 2.988*10^(-26); % Molecular mass of water [kg] (Karagiannis)
v_m = 2.989*10^(-29); % Molecular volume of water [ml][m^3]
k = 1.38064852*10^(-23); % Boltzmann constant [m^2kgs^-2K^-1]
T_C = 647.15; % Cite: n.4 (DataPaper_WölkStrey)
% Inputs %
% Temperature data retrieved from "DataPaper"
T_range = [259.86 249.75 239.58 229.52 218.94];
%T_range = [260 250 240 230 220];
% Nucleation data retrieved from "DataPaper"
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
% Super-Saturation data retrieved from "DataPaper"
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
w_data_25986 = [ 0.02641 .02646 .02623 0.022653 .02587 .0264 .02624 .02651];
% p0 data retrieved from "DataPaper"
p0_data_25986 = {[88.24 88.24 88.24 92.88 92.88],[97.77 97.77 97.77 92.88 92.88 92.88],[97.77 97.77 97.77 97.77 92.88 92.88 92.88], ...
    [97.77 97.77 97.77 97.77 97.77 92.88],[97.77 97.77 97.77 97.77 92.88] [97.77 97.77 97.77 92.88 92.88 92.88],[92.88 92.88 92.88 92.88 92.88],[97.77 97.77 97.77 102.92 102.92]};
p_v_data_25986 = [];
% Vapor pressure computation according to Eq. 15 from "Data paper" 
% Data 260 K case %
for j = 1:length(w_data_25986)
    p_v_data_25986 = [p_v_data_25986, w_data_25986(j) .* p0_data_25986{j}];
end
p_v_data_25986 = p_v_data_25986*10^3;% kPa -> Pa    
% Data 250 K case %
w_data_24975 = [.01677 .01652 .01677];
p0_data_24975 = {[ 91.11 91.11 91.11 86.55 86.55 86.55 82.23 82.23 82.23 82.23],[91.11 91.11 91.11 91.11  86.55 86.55 86.55 86.55 82.23 82.23 82.23 82.23],[91.11 91.11 91.11 86.55 86.55 86.55 86.55 82.23 82.23 82.23 80.17]};
p_v_data_24975 = [];
for j = 1:length(w_data_24975)
    p_v_data_24975 = [p_v_data_24975, w_data_24975(j) .* p0_data_24975{j}];
end
p_v_data_24975 = p_v_data_24975*10^3;% kPa -> Pa    
% Data 240 K case %
w_data_23958 = [.006797 .006808 .00665 .006862 .006917 .006827];
p0_data_23958 = {[99.45 99.45 99.45 99.45 99.45], [99.45 99.45 99.45 104.68],[ 110.19 110.19 110.19 110.19 113.02],[ 113.02 113.02 113.02 107.37 107.37],[ 107.37 107.37 107.37 103.08 103.08 103.08 97.92],[ 103.08 103.08 103.08 96.9 96.9 96.9]}; 
p_v_data_23958 = [];
for j = 1:length(w_data_23958)
    p_v_data_23958 = [p_v_data_23958, w_data_23958(j) .* p0_data_23958{j}];
end
p_v_data_23958 = p_v_data_23958*10^3;% kPa -> Pa    

for i=1:length(T_range)
    T = T_range(i); % Temperature
    x = (T-225)/46.2; % Formula retrieved from "DataPaper"
    t_r = (T_C - T)/T_C; % Formula retrieved from "DataPaper"
    rho(i) = 0.08*tanh(x) + 0.7415*t_r^0.33 + 0.32; % Formula retrieved from "DataPaper"
    v_m = m_m./rho; % Molecular volume of water [cm^3]
    v_m = v_m*10^(-6);
    %sigma = (235.8*(1-T/T_C)^1.256*(1-0.625*(1-T/T_C)))*10^(-3); % Surface tension [N/m] "Revised Release on Surface Tension of Ordinary Water Substance"
    sigma(i) = (93.6635 + 0.009133*T - 0.000275*T^2)*10^(-3); % Formula retrieved from "DataPaper"
    %p_eq = 610.8.*exp(-5.1421.*log(T./273.15)-6828.77*(1./T-1/273.15)); % Eq. Vapor pressure [pa]
    p_eq(i) = exp(77.34491-7235.42465/T - 8.2*log(T) + 0.0057113*T); % Formula retrieved from "DataPaper" cite: n.66 (DataPaper_WölkStrey)
    %p_eq = 133.3224*10^(19.301142-2892.3693/T - 2.892736*log10(T)-0.0049369728*T+5.606905*10^(-6)*T^2-4.645869*10^(-9)*T^3+3.7874*10^(-12)*T^4); % Viisanen eq.
    S = {[4:0.1:10],[5:0.1:13],[7:0.1:16],[9:0.1:19],[8:0.1:27]};
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
    p_v = p_eq(i) .* S{i}; % Vapor pressure [pa] 
    %p_v = 0.345*[106.51:2.1:127.73]*10^3;
    J_BD = sqrt((2.*sigma(i))./(pi.*m)).*v_m.*(p_v ./ (k.*T)).^2.*exp((-16.*pi.*v_m.^2.*sigma(i).^3)./(3.*(k.*T).^3 .* log(S{i}).^2)); % Classical Nucleation Rate (Döring Becker)
    J_H2O = J_BD .* exp(-27.56 + 6.5*10^3/T); % Corrected Nucleation rate 
    figure(1);
    loglog(S{i},J_BD/(10^6),'r--')
    hold on
    loglog(S{i},J_H2O/(10^6),'b-')
    legend('J_{BD}','J_{H2O}')
    ylabel('Nucleation [cm^{-3}s^{-1}]')
    xlabel('Super-Saturation')
    
end
% S_data_25986 = p_v_data_25986 ./ p_eq(1);
% S_data_24975 = p_v_data_24975 ./ p_eq(2);
% S_data_23958 = p_v_data_23958 ./ p_eq(3);
%S_data_22952 = p_v_data_24975 ./ p_eq(2);
figure(1)
hold on
semilogy(S_data_25986, J_data_25986, 'k*')
semilogy(S_data_24975, J_data_24975, 'k*')
semilogy(S_data_23958, J_data_23958, 'k*')
semilogy(S_data_22952, J_data_22952, 'k*')
semilogy(S_data_21894, J_data_21894, 'k*')
grid on
ylim([10^5 10^10]) 
xlim([5 30])
figure(2)
ax1 = axes; 
yyaxis left                 % see [1]
plot(T_range,p_eq)
ylabel('p_{eq}')
hold on
pause(0.1)                  % see [3]
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
grid(ax1,'on')
ytick = ax1.YTick;  
yyaxis right                % see [1]
plot(T_range, sigma)
ylabel('\sigma')
ax2 = axes('position', ax1.Position);
plot(ax2,T_range, rho, 'k')
ylabel('\rho')
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
grid(ax2, 'on')
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}); 
%legend('p_{eq}','\sigma','\rho')

%% Verification plot of parameters %%
T = [220:2.5:320];
T_C = 647.15; % Formula retrieved from "DataPaper" Cite: n.4 (DataPaper_WölkStrey)
x = (T-225)./46.2; % Formula retrieved from "DataPaper"
t_r = (T_C - T)./T_C; % Formula retrieved from "DataPaper"
rho = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32; % Formula retrieved from "DataPaper"
v_m = m_m./rho; % Molecular volume of water [cm^3]
v_m = v_m*10^(-6);
sigma = (93.6635 + 0.009133*T - 0.000275*T.^2).*10^(-3); % Formula retrieved from "DataPaper", cite: n.4,7 (DataPaper_WölkStray)
p_eq = exp(77.34491-7235.42465./T - 8.2*log(T) + 0.0057113*T); % Formula retrieved from "DataPaper", cite: n.66 (DataPaper_WölkStrey)
figure(3)
ax1 = axes; 
yyaxis left                 % see [1]
semilogy(T,p_eq)
ylabel('p_{eq}')
ylim([10^0 10^5])
hold on
pause(0.1)                  % see [3]
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
grid(ax1,'on')
ytick = ax1.YTick;  
yyaxis right                % see [1]
plot(T, sigma)
ylabel('\sigma')
ylim([0.065 0.085])
ax2 = axes('position', ax1.Position);
plot(ax2,T, rho, 'k')
ylabel('\rho')
ylim([0.9 1.15])
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
grid(ax2, 'on')
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'       '}); 
%legend('p_{eq}','\sigma','\rho')        
figure(3)
semilogy(T,p_eq)
display('run finished')