% Retrieve the Schmidt 2007 data %
% Z/D % Pi [S/D]^2 % RHO [kg/m^3] % T [K] % v [m/s] % c [m/s] % R [m] % f % gamma_nuc [1/m^3/s] % supersat
function [x, A, rho, T, V, c, R, f, gn, S, data, r, L, dx, Length, M, p, n] = getData(channelNr, Rg, n, n1, n2, grid_tuning)
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
p   = rho .* Rg .* T; % Pressure 
%D = 90/x(end); % Scaling-Factor
D = 1;
L = x(end);
x = x ./ x(end);


if grid_tuning == true
    A_throat = min(A);
    throat = find(A_throat == A); %n-loaction of throat
    if length(throat)>1
        display("Warning: 2 locations for throat found")
    end
    throat = throat(1);
    x1 = x(1:throat-1);
    A1 = A(1:throat-1);
    A1 = interp1(x1,A1,x1(1):(x1(end)-x1(1))/(n1-1):x1(end));
    x1 = x1(1):(x1(end)-x1(1))/(n1-1):x1(end);
    x2 = x(throat:end);
    A2 = A(throat:end);
    A2 = interp1(x2,A2,x2(1):(x2(end)-x2(1))/(n2-1):x2(end));
    x2 = x2(1):(x2(end)-x2(1))/(n2-1):x2(end);
    %x = 0:1/(n-1):1;
    A = [A1 A2];
    %A = A./A(throat); % [non-dimensional]
    x = [x1 x2];
    n = n1+n2;
    for i = 1:length(x)-1
        dx(i) = x(i+1)-x(i);
    end
else
    A = interp1(x,A,0:1/(n-1):x(end));
    A_throat = min(A);
    throat = find(A_throat == A); %n-loaction of throat
    if length(throat)>1
        display("Warning: 2 locations for throat found")
    end
    throat = throat(1);
    %A = A./A(throat); % [non-dimensional]
    x = 0:1/(n-1):1;
    Length = x(end);
    dx = Length/(n-1); 
end
r = sqrt(A/pi); % [non-dimensional]
Length = x(end);
%L = 150; % Scaling Factor for dimensionalize.m 