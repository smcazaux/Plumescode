% Retrieve the Schmidt 2007 data %
% Z/D % Pi [S/D]^2 % RHO [kg/m^3] % T [K] % v [m/s] % c [m/s] % R [m] % f % gamma_nuc [1/m^3/s] % supersat
function [x, A, rho, T, V, c, R, f, gn, S, data, r, L, n, dx, Length, M, p, x_d, dx_d] = getData_Verification(channelNr)
if channelNr == 100
    data = load('data/Example_0_numeric.txt');
elseif channelNr == 101
    data = load('data/Example_1_numeric.txt');
end
data = transpose(data);
x   = data(1,:); %[Z/D]
A   = data(2,:); % Pi [S/D]^2
rho = data(3,:); % [kg/m^3]
T   = data(4,:); % [K]
V   = data(5,:); % [m/s]
c   = data(6,:); % [m/s]
R   = data(7,:); % [m]
f   = data(8,:); % [-]
gn  = data(9,:); % [1/m^3/s]
S   = data(10,:);% [-]
M   = V ./ c; % Mach number [-]
m0 = 2.992*10^(-26);
kb = 1.38065*10^(-23);
Rg = kb/m0;
p   = rho .* Rg .* T; % Pressure
L = x(end);
% D = 1; %Common scaling factor
% x = x * D;
% A = A * D^2;
x = x ./ x(end);
r = sqrt(A/pi); % [m]
Length = x(end);
n = length(x);
dx = Length/(n-1);
x_d = x.*L;
dx_d = dx.*L;