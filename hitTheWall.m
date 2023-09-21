% Chance of particle hitting a wall %
function E = getAccretion(V_d, Rg, T_d, p_d) 
k = 1.38064852*10^(-23); % Boltzmann constant
Dmol_H2O = 2.75*10^(-10); % Molecular diameter H2O
c = k/(pi*Dmol_H2O^2); % Constant
Lcoll = c.*T_d./p_d; % mean free path
crackLength = 500e3; % Lenght of crack = 500 km
m_sides = c.*crackLength.*V_d./Rg; % Mass flow at each side of the wall
for i = 1:length(Lcoll)
    massflow{i} = linspace(0,m_sides(i),1000); % kg/s
    h{i} = linspace(Lcoll(i),0,1000);
    Vx{i} = pi.*h{i}.^2.*(Lcoll(i)-h{i}./3);
    V = 4/3*pi*Lcoll.^3;
    %d{i} = r-h; % distance to wall [m] 
end

chance = Vx./V;
figure(1)
plot(d,chance)
ylabel('% chance to hit the wall ')
xlabel('distance to the wall [m]')
figure(32)
plot(massflow,chance)
ylabel('% chance to hit the wall ')
xlabel('mass flow in the channel (accumulating from left to right) [kg/s]')
totalE = trapz(massflow,chance)
averageChance = sum(chance(:))/numel(chance);
E = averageChance*massflow(end)