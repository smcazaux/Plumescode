% Chance of particle hitting a wall %
function [wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, run_name, k, nt, predictor, m_stick, T_w, c_d, r_sub) % [kg/s]
global approach1 mm_sigma mmm_sigma wallInteractions
if r_sub == 1 % (Cazaux)
    kkk = 4.96*10^15.*exp(-5705./T_w); 
elseif r_sub == 2 % (Speedy) 
    E = 120*48.25;
    kkk = 3.99*10^15.*exp(-E./T_w);
elseif r_sub == 3 % (Hasegawa & Herbst)
    kkk = 2.29*10^12.*exp(-5600./T_w); 
elseif r_sub == 4 % (Smith 2015)
    E = 120*48.25;
    kkk = 10^15.*exp(-E./T_w);
end

Dmol_H2O = 2.75*10^(-10); % Molecular diameter H2O
Angstrom = 1*10^-10; % Angstrom in [m]
n_molecules = 1/(3*Angstrom)^2; % number of molecules in [n/m2]
Area = c_d.*dx_d; % Circumference*dx: Actual area [m2]
rate_m = kkk.*n_molecules.*Area; % number of molecules that sublimate per s [n/[s]
N_A = 6.022140857*10^23; % number of Avogadro
m_H2O = 0.018; % Molar weight H2O [kg/mol] 
sublimation_d = rate_m./N_A.*m_H2O; % [kg/s]
mdot = rho_d.*A_d.*V_d;
kk = 1.38064852*10^(-23); % Boltzmann constant
sigmaa = pi*Dmol_H2O^2 ;
c = kk/(sigmaa); % Constant
Lcoll = c.*T_d./p_d; % mean free path
Vth = 1.5.*sqrt(T_d./100).*(1./sqrt(18)).*1000; % [m/s]
dt = kk.*T_d./(sigmaa.*Vth.*p_d);
t = dx_d./V_d; % 1 cell to 2 cell
n = round(t./dt); % nr. of coll
average_chance_sides = 0.187562562562563; % Chance to hit the wall
chance_channel = average_chance_sides.*2.*Lcoll./(2.*r_d); % Average chhance over whole channel to hit the wall
if m_stick == 1 % Cazaux 2011
    c_stick = (1+4.2*10^-2*sqrt(T_d+T_w) + 2.3*10^-3*T_d - 1.3*10^-7*T_d.^2).^-1; 
elseif m_stick == 2 % Vijay 2014
    y = 244./T_d;  
    A = 0.51; % over-estimated due to prefactor A (based on Tice = 70 K)
    c_stick = A.*(y.^2 + 0.8.*y.^3)./(1 + 2.4.*y + y.^2 + 0.8.*y.^3) ;% A = 0.51 (pre-factor for ice at 70 K) 
elseif m_stick == 3 % Buch 1991
    c_stick = ((T_d./102) + 1).^-2; 
elseif m_stick == 4
    c_stick = 1;
end
             
x = chance_channel.*c_stick; % Chance at every location of the channel 
% Approach when nt > 30.000
if k == 1 && nt > 1400000 && predictor == true% takes ~7h 
    approach1 = true;
    for i = 1:length(x)
        nn{i} = linspace(1,round(n(i).*1.4),round(n(i).*1.4));
        mm_sigma{i} = x(i).*(1-x(i)).^(nn{i}-1);
        for j = 1:length(mm_sigma{i})
            mmm_sigma{i}(j) = sum(mm_sigma{i}(1:j));
        end
    end
%     for i = 1:length(x)
%         mm_sigma{i} = x(i).*(1-x(i)).^(nn-1);
%     end
elseif k == 1 && nt <= 800000 && predictor == true% NEW APPROACH %
    approach1 = false;
    for i = 1:length(x)
        nn{i} = linspace(1,n(i)*10,n(i)*10);
        mm_sigma{i} = x(i).*(1-x(i)).^(nn{i}-1);
    end
end
if approach1 == true
    for i=1:length(n)
        if  n(i) > length(mmm_sigma{i})
            nn = linspace(1,n(i),n(i));
            m_sigma(i) = sum(x(i).*(1-x(i)).^(nn-1));
        else
            m_sigma(i) = mmm_sigma{i}(n(i));
        end
    end
elseif approach1 == false 
    for i = 1:length(x)
        if  n(i) > length(mm_sigma{i})
            nn = linspace(1,n(i),n(i));
            m_sigma(i) = sum(x(i).*(1-x(i)).^(nn-1));
        else
            m_sigma(i) = sum(mm_sigma{i}(1:n(i)));
        end
    end
end
% OLD APPROACH %
if wallInteractions == false
    for i = 1:length(n)
        nn = linspace(1,n(i),n(i));
        m_sigma(i) = sum(x(i).*(1-x(i)).^(nn-1));
    end
end
accretion_d = m_sigma.*mdot; % [kg/s]
wallMassFlux_d = accretion_d - sublimation_d; % [kg/s]
E = wallMassFlux_d ./ (crackLength.*dx_d); % [kg/m2s]
%%
if k==nt
    accretion = figure(991);
    subplot(6,1,1)
    semilogy(x_d,accretion_d)
    ylabel('ACC')
    subplot(6,1,2)
    plot(x_d,sublimation_d)
    ylabel('Sub')
    subplot(6,1,3)
    plot(x_d,n)
    ylabel('n')
    subplot(6,1,4)
    plot(x_d, wallMassFlux_d )
    ylabel('Wall mass flux [kg/s]')
    subplot(6,1,5)
    plot(x_d,E)
    ylabel('E [kg/m^2s]')
    subplot(6,1,6)
    semilogy(x_d,m_sigma);
    xlabel('x')
    ylabel('\sigma')
    name = ['accretion.png'];
    saveas(accretion, ['Figures/' run_name '/' name]);

end
    %%
% for i = 1:length(Lcoll)
%     h{i} = linspace(Lcoll(i),0,1000);
%     Vx{i} = pi.*h{i}.^2.*(Lcoll(i)-h{i}./3);
%     V(i) = 4/3*pi*Lcoll(i).^3;
%     chance{i} = Vx{i}./V(i);
%     chancee(i) = sum(chance{i})/numel(chance{i}); %average chance near walls
%     chanceee(i) = chancee(i)*2*Lcoll(i)./(2.*r_d); % Average chhance over whole channel
%     wallMassFlux_dt(i) = chanceee(i).*mdot    
%     wallMassFlux_d(i) = chancee(i)*m_sides(i)
%     %wallMassFlux_d(i) = n(i)*chancee(i)*m_sides(i);
% end  
% %%
% for i = 1:length(Lcoll)
%     massflow{i} = linspace(0,m_sides(i),1000); % kg/s
%     h{i} = linspace(Lcoll(i),0,1000);
%     Vx{i} = pi.*h{i}.^2.*(Lcoll(i)-h{i}./3);
%     V(i) = 4/3*pi*Lcoll(i).^3;
%     chance{i} = Vx{i}./V(i);
%     d{i} = Lcoll(i)-h{i}; % distance to wall [m] 
%     wallMassFlux_d(i) = trapz(massflow{i},chance{i}); % [kg/s]
% end
%%
% 
% 
% ylabel('% chance to hit the wall ')
% xlabel('distance to the wall [m]')
% figure(32)
% plot(massflow,chance)
% ylabel('% chance to hit the wall ')
% xlabel('mass flow in the channel (accumulating from left to right) [kg/s]')
% totalE = trapz(massflow,chance)
% averageChance = sum(chance(:))/numel(chance);
% E = averageChance*massflow(end)



% %%
% xdx = linspace(0,1,1000);
% figure(321)
% yline(average_chance_sides, 'k--') 
% hold on 
% plot(xdx,chance{1},'linewidth', 1.5)
% legend('Average \approx 0.19')
% grid on
% ylabel('Probability [-]')
% xlabel('\lambda')
end