clear all, close all
% Get verification data
channelNr = 100; % Verification Example_0 = 100, Example_1 = 101
[x, A, rho, T, V, c, R, f, gn, S, data, r, L, n, dx, Length, M, p, x_d, dx_d] = getData_Verification(channelNr)
% R = [m]
r_e = R(end) - R; % Eq. (S15) Schmidt; Gives particle size at exit [m], nucleated from location z
figure(7)
semilogy(x_d,r_e);
ylabel('Particle Radius exit, generated from nucleation at z [m]')
xlabel('z')
% gn = [1/(m^3s)], dx_d = [m], V = [m/s]
n_z = gn.*dx_d./V; % Eq. (S16) Schmidt; Concentration of such particles the point of their origin
figure(9)
semilogy(x_d, n_z)
ylabel('Particle concentration at the point of their origin [n/m^3]')
xlabel('z')
drdz = gradient(r_e(:)) ./ gradient(x_d(:)); % change in exit particle size versus change in location along z
for i = 1:length(r_e)-1
    dr(i) = r_e(i+1)-r_e(i); % difference between exit radius 'dr' [m]
end
dr = [dr 0]; % [m]

% Three different Equations, first one has unit [1/m^4], the other two [1/m^3]

% (S18) divided by dr
P = gn ./ V(end) .* abs(transpose(drdz.^(-1))) .* A./A(end); % [1/m^4]
% My own derivation 3.40
%P = gn ./ V(end) .* abs(transpose(1./drdz)) .* A./A(end) .*dx_d; % [1/m^3]
% (S18) with dr
%P = gn ./ V(end) .* abs(transpose(drdz.^(-1))) .* A./A(end).* abs(dr); % [1/m^3]

figure(11)
semilogx(r_e*10^6, P) % Radius in [um]
title('Particle Size Distribution')
ylabel('Particle concentration [1/m^4] ')
xlabel('R(\mu m)_{exit}')