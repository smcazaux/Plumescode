function [M Tstat pstat rho Tres Pres rhores Vel] = mini_moon(n, channelNr, x, A, dx, R)  
Tres = 273.16;               % reservoir Temperature
Pvap = 10^(8.07131+2.124903-(1730.63)/(233.426-273.15+Tres)); %Antoine equation (https://en.wikipedia.org/wiki/Antoine_equation) vapor pressure
Pres = 771.2; % set to generate sub-sonic -> super-sonic flow, does only influence pstat by multiplication

Pres = 611.2;
rhores = 4.85e-3;
%Pres/(R*Tres);
Pvacz = [100];% vacuum pressure (490 K)
P_exit_ratio = Pvacz/Pres; % pressure ratio: vacuum/reservoir
Legend=cell(length(Pvacz),1); % Create empty cell 
%% % NOZZLE 1 Anderson % 
% n = 375;
% Rg = 461,52;
% channelNr = 101;
% [x, dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9, dum10, r, dum11, dx, dum13, dum14, dum15] = getData(channelNr, Rg, n);

% we need: x, r
%[x, A, rho_initial, T_initial, V_initial, dum4, dum5, dum6, dum7, dum8, dum9, r, L, dx, Length, M_initial, p_initial] = getData(channelNr, Rg, n);
        
% Length = 1;
% L = 150;    %[m]
% n = 1500;
% x = linspace(0,Length,n);
% dx = Length/(n-1);
% % % Schmidt area crevasse %
% % r = -33.99 .* x.^5 + 104.1 .* x.^4 - 114.3 .* x.^3 + 52.68.*x.^2 - 8.869 .* x + 1.002;
% 
% % Schmidt Channel 1 %
% p1 = 2425.9223;
% p2 = -12291.7689;
% p3 = 26015.6998;
% p4 = -29617.7449;
% p5 = 19409.4629;
% p6 = -7217.5659;
% p7 = 1363.4469;
% p8 = -80.2421;
% p9 = -7.6893;
% p10 = 0.42572;
% p11 = 0.49841;
% r = p1*x.^10 + p2*x.^9 + p3*x.^8 + p4*x.^7 + p5*x.^6 + p6*x.^5 + p7*x.^4 + p8*x.^3 + p9*x.^2 + p10*x + p11 ;
r = sqrt(A./pi);
r_throat = min(r);
throat = find(r_throat == r); %n-loaction of throat
throat = throat(1); % Take first one if more than 1
x_throat = x(throat);
%x_throat = (throat(1)-1)*dx
%%
% A = pi .* r.^2;
% figure(7)
% plot([x_throat x_throat],[r_throat -r_throat],'b--','linewidth',2)
% hold on;
% plot(x,-r,'k','linewidth',2)
% plot(x,r,'k','linewidth',2)
% ylabel('Scaled Nozzle Radius [-]')
% xlabel('x/L [-]') 
% legend('Throat')
% set(gca,'fontsize',14)
% grid on
Ares = A(1);                % Might not be correct
%Ares = pi*(90e-3/2)^2;      % reservoir surface
Vres = Ares * 100e-3;       % reservoir volume    
%Athr = 3.9e-3*90e-3;        % throat area    
%Aex = 1.5*Athr;             % exit area
Aex = A(end);
Athr = A(throat);
alpha = 0.7;                % sticking coefficient
g = 4/3;                    % specific heat ratio Cp/Cv, close approximation at T =275 K
                            % H = Cp*T
                            % U = Cv*T

%Equilibrium condition (reservoir pressure)
Pres_eq = Pvap/(Athr/Ares * sqrt(2*pi*g*(2/(g+1))^((g+1)/(g-1)))/alpha+1);
% Isentropic 1D-flow, only change in cross-sectional area, perfect gas, 
% T0/T1 = temperature at location of interest divided by temperature at
% stagnation point. (Same for P0/P1)

% Determine flow conditions at throat for SUB and SUP %
% This assumes M = 1 in throat %
[dum1,dum2,pe3,dum4,dum5] = flowisentropic(g,Aex/Athr,'sub'); %flowisentropic(gamma, flow, mtype)
% pe3 = pressure ratio (subsonic to sonic, M<1 to M=1)
[M6,dum2,pe6,dum4,dum5] = flowisentropic(g,Aex/Athr,'sup'); % is this in equillibrium?
pe5 = pe6*(1+2*g/(g+1)*(M6^2-1)); % shockWave Equation
P_exit_ratio = pe5 - 0.001;
%pe5 = pe6*((2*g*M6^2-(g-1))/(g+1))
Press = 3000;
dt = 0.00002;               % time step
iend = 1000;                % 1000*0.00002 = 0.02 s 

for z = 1:length(Pvacz)
    Pvac = Pvacz(z);
    for i = 1:iend
        p(i) = Pres;
        if P_exit_ratio < pe3
        %if Pvac/Pres < pe3 %Pressure is below pe3, CHOKED FLOW (critical pressure ratio)(ALWAYS)
            Astr(i) = Athr; % = A* = 3.5100e-04 m^2]
            if P_exit_ratio > pe5
            %if Pvac/Pres > pe5 %NSW in nozzle
                % Compute Mach number downstream of shock, 
                Mach(i) = sqrt(1/(g-1) * ( sqrt(1+2*(g-1)*(2/(g+1))^((g+1)/(g-1))*(Astr(i)/Aex*Pres/Pvac)^2 )-1)); %Mach number in the exit
                % Compute A/A* that matches Mach number downstream of shock
                % (subsonic flow)
                [dum1,dum2,dum3,dum4,Aex_Astr2(i)] = flowisentropic(g,Mach(i),'mach'); %Virtual second throat needed for the shock
                % Compute pressure loss, Psup/Psub at exit
                p02_p01(i) = Athr/Aex*Aex_Astr2(i); %Total pressure loss over the shock = A*1/A*2
                [Ms1(i), dum2, dum3, dum4, Ms2(i), dum6, dum7] = flownormalshock(g, p02_p01(i), 'totalp');
                [dum1,dum2,dum3,dum4,Aex_Astr2_afterShock(i)] = flowisentropic(g,Ms2(i),'mach'); %Virtual second throat needed for the shock
                [dum1,dum2,dum3,dum4,As_Astr(i)] = flowisentropic(g,Ms1(i),'mach'); %Area cross section at the location of the shock M>1
                [dum1,dum2,dum3,dum4,As_Astr_forMst1(i)] = flowisentropic(g,Ms2(i),'mach'); %Area cross section at the location of the shock M<1
            elseif P_exit_ratio <= pe5
            %elseif Pvac/Pres <= pe5
                Mach(i) = M6;   % from super sonic relations
                                % If Pvac/Pres > Pb/P0|sup && Pvac/Pres < Pb/P0|sub
                                % -> Oblique shockwave.
                                % If Pvac/Pres < Pb/P0|sup
                                % -> Expansion fan.
                                % DOESN'T AFFECT FLOW IN CD NOZZLE
                %As_Astr = Aex/Athr+10.000001;
                As_Astr = Aex/Athr+inf; % Set shock wave location outside the channel length
            end
            
        elseif Pvac/Pres > pe3 && Pvac/Pres <1
        %elseif Pvac/Pres > pe3 && Pvac/Pres <1 %fully subsonic
            [Mach(i),dum2,dum3,dum4,Astr(i)] = flowisentropic(g,Pvac/Pres,'pres');
            As_Astr = 0;
            %Astr(i) = Aex/A;
            % subsonic flow %
    %         Astr(i) = 
    %         Aex_Astr2 = 1;
    %         As_Astr = 0;
            %As_Astr = Astr(i);
        elseif P_exit_ratio >=1    
        %elseif Pvac/Pres >= 1
            Astr(i) = 0;
            As_Astr(i) = 0;
            Aex_Astr2(i) = (Aex/Athr);
            p02_p01(i) = 1;
            Mach(i) = Mach(i-1);
        end

        if Pvap > Pres % liquid turns into gas
            m_dot_evap(i) = alpha*(Pvap-Pres)*sqrt(1/(2*R*pi*Tres))*Ares;
        else
            m_dot_evap(i) = 0;
        end
        % Mass flow after nozzle (choked flow) % it never stops since
        % equillibrium is reached, try different start state
        m_dot_flow(i) = Pres*Astr(i)/sqrt(Tres)*sqrt(g/R*(2/(g+1))^((g+1)/(g-1)));
        pdot(i) = R*Tres/Vres*(m_dot_evap(i)-m_dot_flow(i));
    %     Pres = Pres+pdot(i)*dt;
    %     Press(i+1) = Press(i)+pdot(i)*dt;
    end

    Temp = Tres./(1+(g-1)/2*Mach.^2); % static T
    U = Mach.*sqrt(g*R*Temp); % Velocity

    k = 1;

    %%
    %for Ar = 5:-0.1:1 %section from -0.5 till 0
    for Ar = 1:throat-1 % until throat
        x_L(k) = x(k);
        Ar = A(k)/Athr;
        %x_L(k) = 1*(Ar-1)/length(A);%0.5/200*Ar;
        %Ar = A(Ar)/Athr;
        [M(k),Tr,pr,dum4,dum5] = flowisentropic(g,Ar,'sub'); % assumes subsonic
        Area(k) = Ar; 

        %x_L(k) = -0.5*(Ar - 1)/(5 - 1);
        Tstat(k) = Tr * Tres;
        pstat(k) = pr * Pres;
        k = k + 1;
    end
    % shock is included
    %for Ar = 1:0.001:Aex/Athr % determine M, T, P, along nozzle, until 1, where Ar = 1.5 
%     A_shock = As_Astr(end);
%     if A > A_shock
%         A_shockk = A;
%     end
%     shock = find(A_shockk(1) == A); %n-loaction of shock
%     x_shock = (shock(1)-1)*dx;
    for Ar = throat:length(A) 
        
        x_L(k) = x(k);
        Ar = A(k)/Athr;
        %         x_L(k) =1*Ar/length(A);%0.5 + 0.5*(Ar-200)/(length(A)-200);
%         Ar = A(Ar)/Athr;
        As_Astr;
%         if As_Astr == inf
%             [M(k),Tr,pr,dum4,dum5] = flowisentropic(g,Ar,'sub');
%             Tstat(k) = Tr * Tres;
%             pstat(k) = pr * Pres;
%         end
        if Ar < As_Astr(end) % stop at shock @eq.
            [M(k),Tr,pr,dum4,dum5] = flowisentropic(g,Ar,'sup');
            Tstat(k) = Tr * Tres;
            pstat(k) = pr * Pres;
        %elseif x_L(k) > x_shock    
        else
            [M(k),Tr,pr,dum4,dum5] = flowisentropic(g,Ar*Aex_Astr2(end)*Athr/Aex,'sub'); %From virtual second throat / 1.5
            Tstat(k) = Tr * Tres;
            pstat(k) = pr * Pres * p02_p01(end); % pressure compensated with pressure loss ratio
        end
        Area(k) = Ar;
        %x_L(k) = (Ar - 1)/(Aex/Athr - 1); % Ar 1-1.5, x_L = 0-1
        k = k + 1;
    end
    rho = pstat ./ (R .* Tstat);
    Vel = M.*sqrt(g*R*Tstat);
    %%
%     figure(11), clf
%     plot(dt*[1:iend],m_dot_evap*1e3,'linewidth',2), hold on
%     plot(dt*[1:iend],m_dot_flow*1e3,'linewidth',2)
%     ylabel(['$\bf{\dot{m}}$ $\bf{[g/s]}$'],'Interpreter','latex')
%     xlabel('\bf{t} [s]')
%     leg1 = legend('$\dot{m}_{evap}$','$\dot{m}_{flow}$');
%     set(leg1,'Interpreter','latex');
%     title(['p_v_a_c = ' num2str(round(Pvac)) ' [Pa], T_r_e_s = ' num2str(round(Tres*10)/10) ' [K], p_v_a_p = ' num2str(round(Pvap)) ' [Pa]'])
%     set(gca,'fontsize',14)
%     grid on
%     %export_fig mdot.pdf -transparent
%     %%
%     figure(22), clf
%     title('P_{res} and Exit conditions over time (Non Steady-State)')
%     %set(gcf,'Position',[100,400,1200,500]);
%     subplot(2,2,1)
%     plot(dt*[1:iend],p,'linewidth',2)
%     ylabel('p_r_e_s [Pa]')
%     xlabel('t [s]')
%     set(gca,'fontsize',14)
%     grid on
%     subplot(2,2,2)
%     plot(dt*[1:iend],Mach,'linewidth',2)
%     ylabel('M_e_x [-]')
%     xlabel('t [s]')
%     set(gca,'fontsize',14)
%     grid on
%     subplot(2,2,3)
%     plot(dt*[1:iend],U,'linewidth',2)
%     ylabel('V_e_x [m/s]')
%     xlabel('t [s]')
%     set(gca,'fontsize',14)
%     grid on
%     subplot(2,2,4)
%     plot(dt*[1:iend],Temp,'linewidth',2)
%     ylabel('T_e_x [K]')
%     xlabel('t [s]')
%     set(gca,'fontsize',14)
%     grid on
%     %export_fig exit_time.pdf -transparent
% 
%     figure(33), clf
%     plot(dt*[1:iend],As_Astr)
%     ylabel('A_s / A* [-]')
%     xlabel('t [s]')
%     grid on
%     legend_text = ['P_{exit} = ', num2str(round(Pvac,0)) ,' Pa']
%     Legend{z}= legend_text;
%     figure(444)
%     %set(gcf,'Position',[100,400,1200,500]);
%     subplot(2,2,1)
%     hold on
%     plot(x_L,pstat,'linewidth',2);
%     ylabel('p_s_t_a_t [Pa]')
%     xlabel('x/L [-]')
%     set(gca,'fontsize',18)
%     grid on
%     if z == length(Pvacz)
%         legend(Legend,'Location','best');
%         legend('boxoff');  
%     end
%     subplot(2,2,2)
%     hold on
%     plot(x_L,Tstat,'linewidth',2);
%     ylabel('T_s_t_a_t [K]')
%     xlabel('x/L [-]')
%     set(gca,'fontsize',18)
%     grid on
%     subplot(2,2,3)
%     hold on
%     plot(x_L,M,'linewidth',2);
%     ylabel('Mach [-]')
%     xlabel('x/L [-]')
%     set(gca,'fontsize',18)
%     grid on
%     subplot(2,2,4)
%     hold on
%     plot(x_L,rho,'linewidth',2);
%     ylabel('\rho [kg/m^3]')
%     xlabel('x/L [-]')
%     set(gca,'fontsize',18)
%     grid on
%     %export_fig flow_channel.pdf -transparent
% 
%     figure(55), clf
%     semilogy(Tstat, pstat,'linewidth',2), hold on
%     semilogy([200:300],10.^(8.07131+2.124903-(1730.63)./(233.426-273.15+[200:300])),'k','linewidth',2);
%     semilogy([273.16 273.16],[612 1e4],'k:','linewidth',2);
%     ylabel('p_s_t_a_t [Pa]')
%     xlabel('T_s_t_a_t [K]')
%     title(['p_r_e_s = ' num2str(round(Pres)) ' [Pa], T_r_e_s = ' num2str(round(Tres*10)/10) ' [K], p_v_a_c = ' num2str(round(Pvac)) ' [Pa]'])
%     set(gca,'fontsize',14)
%     grid on
% 
%     figure(66)
%     hold on
%     plot([x_throat,x_throat], [0.5*Athr, -0.5*Athr],'k--')
%     plot(x_L, Area*Athr/2, x_L, -(Area*Athr/2),'k','LineWidth',1.5)
%     plot(x_L, Area*Athr/2, x_L, (Area*Athr/2),'k','LineWidth',1.5)
%     legend('Throat')
%     ylabel('Area CD Nozzle')
%     xlabel('x_L')
%     set(gca,'fontsize',14)
%     %export_fig pres_temp_diagram.pdf -transparent
%     pstat(1), pstat(end), M(end)
end
end