% Retrieve the Schmidt 2007 data %
% Z/D % Pi [S/D]^2 % RHO [kg/m^3] % T [K] % v [m/s] % c [m/s] % R [m] % f % gamma_nuc [1/m^3/s] % supersat
function [x, A, rho, T, V, c, R, f, gn, S, data, r, L, n, dx, Length, M, p] = getDataWithOffset(channelNr, Rg)
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
M   = V ./ c; % Mach number
for i = 1:length(rho)% create a maximal of +- 10% offset for all data points
    rho(i) = rho(i)*(1+(rand()-0.5)*0.2); 
    T(i) = T(i)*(1+(rand()-0.5)*0.2);
    M(i) = M(i)* (1+(rand()-0.5)*0.2);
end
p   = rho .* Rg .* T; % Pressure 
D = 90/x(end); % Scaling-Factor
%x = x*D;
x = x ./ x(end);
A = A*D^2;
r = sqrt(A/pi); % [m]
%r = 0.5*( r / max(r));
%A = pi*r.^2;
Length = x(end);
n = length(x);
dx = Length/(n-1);
L = 90; % Scaling Factor for dimensionlaize.m 


