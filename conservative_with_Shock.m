%conservative - Solves conservative form of govering equations 

function [m, M, V, T, p, rho, m_d, x_d, rho_d, T_d, r_d] = conservative_with_Shock(x, dx, n, C, A, nt, y, throat, Cx, Rg, exit, pe, run_name, r_res_d, L, cv, r_throat, T_ref, rho_ref, p_ref, nucleation, channelNr, verification, b, M_ref, r, c, c_d, r_sub)
global convectionFriction cp wallInteractions m_stick approach1 mm_sigma mmm_sigma crackLength circShape gn E M U1 U2 U3 J2 Edot_L p rho T M b simp offset V0_fix convergenceCriteria critRes rho_m p_m T_m M_m f R fp x_d dR_dx saveAllData rho_reduction wallInteractions testrun As A_d r_d stop
variation = C;                      % Set parameter to vary (sensitivity study)
Legend=cell(length(variation),1);   % Create empty cell for parameter that is varied 
% Constants:
rho_ice = 920;                      % Ice density [kg/m3]
H2O_mol_mass = 0.01801528;          % Molar mass H2O [kg/mol]
avogadro = 6.02214076*10^(23);      % nr. of molecules in one mole
L_h = 2.8e6;                        % Latent Heat [J/kg] 
% Throat location (dimensional)
x_d = x .* L;
x_d_throat = [x_d(throat)];
% from Mini_moon.m (non-dimensional) %
rho = rho_m;
T = T_m;
p = p_m;
M = M_m;
V = M./ sqrt(T);
figure(2) % Plot of initial flow profile  
subplot(5,1,1)
plot(x_d,rho)
hold on
xline(x_d_throat,'k--')
ylabel('\rho')
subplot(5,1,2)
plot(x_d,T)
hold on
xline(x_d_throat,'k--')
ylabel('T')
subplot(5,1,3)
plot(x_d,V)
hold on
xline(x_d_throat,'k--')
ylabel('V')
subplot(5,1,4)
plot(x_d,p)
hold on
xline(x_d_throat,'k--')
ylabel('p')
subplot(5,1,5)
plot(x_d,M)
hold on
xline(x_d_throat,'k--')
ylabel('M')
a0 = sqrt(y*Rg*T_ref);
T_w = linspace(273,236,n); % Ingersoll (& Bouqet)
% Solution vectors
U1 = rho .* A;
U2 = rho .* A .* V;
U3 = rho .* A .* (T / (y - 1) + y / 2 * V.^2); 
if length(dx)>1 
    dx = [dx(1), dx];
end
for l = 1:length(variation) % Define parameter to vary
    C_ = variation(l);
    %C_ = C;
    p_e = pe;
    %p_e = variation(l);
    %C_x = variation(l); 
    C_x = Cx;
    ii=16;
    
    
% Iteration loop
    for k = 1 : nt   

        if mod(k,100)==0
            fprintf('firstrd3=%f',r_d(150))
            fprintf('dt=%f',dt)
            fprintf('k=%f',k)
            channelRadius_d = figure(ii);
%            plot([x_d_throat x_d_throat],[r_d(throat) -r_d(throat)],'k--')
            plot([x_d_throat x_d_throat],[8 -8],'k--')
            hold on
            plot(x_d,r_d,'k','linewidth',2)
            plot(x_d,-r_d,'k','linewidth',2)
            hold off
            r_max=max(r_d)
            D_min_max = r_throat/r_max;
            D_exit = r_d(end)*2;
            D_throat = 2*r_d(throat);
            legend(['D_{min}/D_{max} = ', num2str(round(D_min_max,2))],['D_{exit} = ', num2str(round(D_exit,2)), ' m'],['D_{throat} = ', num2str(round(D_throat,2)),' m'] );
            ylabel('Channel Radius [m]')
            xlabel('z [m]')
            ii=ii+1
            accretion = figure(ii+50);

    subplot(711)
    plot(x_d,V_d,'b','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,V_verification,'--r','linewidth',2)
        if offset == true
            plot(x_d,V_offset,'--y','linewidth',1.5) 
        end
    end
    ylabel('V [m/s]')
    xlabel('z [m]')
    grid minor
    hold off
    subplot(712)
    plot(x_d,T_d,'g','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,T_verification,'--r','linewidth',2)
        if offset == true
            plot(x_d,T_offset,'--y','linewidth',1.5) 
        end
    end
    hold off
    ylabel('T [K]')
    xlabel('z [m]')
    grid minor
    subplot(713)
    plot(x_d,rho_d,'r','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,rho_verification,'--b','linewidth',2)
        if offset == true
            plot(x_d,rho_offset,'--y','linewidth',1.5) 
        end
    end
    hold off
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid minor
    xlabel('z [m]')
    subplot(714)
    plot(x_d,p_d,'m','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,p_verification,'--b','linewidth',2)
        if offset == true
            plot(x_d,p_offset,'--y','linewidth',1.5) 
        end
    end
    ylabel('p [pa]')
    grid minor
    hold off
    xlabel('z [m]')
    subplot(7,1,5)
    plot(x_d,E)
    hold on
    xline(x_d_throat,'k--')
    ylabel('E [kg/m^2s]')
    xlabel('z [m]')
    
    subplot(7,1,6)
    plot(x_d,gn)
    hold on
    xline(x_d_throat,'k--')
    ylabel('nuc. [1/m^3s]')
    xlabel('z [m]')
    subplot(7,1,7)
    plot(x_d,f)
    hold on
    xline(x_d_throat,'k--')
    ylabel('f ')
    xlabel('z [m]')
    
        end
        %     fprintf('rd10=%f',r_d(150))       
        % Copies of solution vectors
        U1_old = U1; 
        U2_old = U2; 
        U3_old = U3;

        % dt calculations
        dxx = min(dx);
        dt = min(C_ .* dxx ./ (sqrt(T) + V));
        %dt=0.01 
        % Flux vectors     
        F1 = U2;
        F2 = U2.^2 ./ U1 + (y - 1) / y * (U3 - (y * U2.^2) ./ (2 * U1));
        F3 = y * U2 .* U3 ./ U1 - y * (y - 1) / 2 * U2.^3 ./ U1.^2; 
        r=r_d;
        [rho_d T_d p_d V_d A_d x_d dx_d r_d] = dimensionalize(rho, T, p, M, y, Rg, A, x, dx, rho_ref, T_ref, p_ref, r_res_d, L, throat, r, As, circShape, crackLength);
        if nucleation == true
            %fprintf('rd11=%f',r_d(150))       
            % Dimensionalize for computation: f', dE, dR. ('_d' = dimensionalized) %
            Q = mean(rho_d .* A_d .* V_d); % [kg/s], average because; mass flow = C  
            dR_dx = dR_dx_MC(T_d,rho_d,M, b); % [-]
            [R dR_dx] = getParticleRadius(dR_dx, x_d, dx_d); % [m], [-]
            [fp gn S] = solid_frac_MC(Q,dx_d,T_d,rho_d,R,dR_dx,A_d,x_d,V_d, n, Rg); % uses rho_d(i), T_d(i), M(i) 
            % Updating solid fraction, Change to Simpson's rule  %
            [f fp] = getSolidFraction(fp, x_d, dx_d);
            % Non-dimensionalize
            [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
            % Calculate energy of latent heat in all segments
            Edot_L = Edot_source(rho, f , fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, y); % Latent heat as Source term
        end
        if rho_reduction == true
            fprintf('rd12=%f',r_d(150))       
            mdot_d  = mdot_sink(f, rho, A, V); % Mass reduction ~Phase change vapor -> ice
        end
        if wallInteractions == true
            %fprintf('rd13=%f',r_d(150))       
            [tau tau_d] = wallStress(rho, V, rho_d, V_d, T_d, r_d, rho_ref, a0); % [-] [kg/ms2] 
            [friction] = getFriction(tau, r, L, A_d, throat, c, dx); % Non-dimensional
            [wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, run_name, k, nt, true, m_stick, T_w, c_d, r_sub); % [kg/s], [kg/m2s]
            [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
            % Non-dimensionalize for predictor step
            wallMassFlux = wallMassFlux_d./(rho_ref.*As.*a0);
            accretionFlux = accretion_d./(rho_ref.*As.*a0);
            sublimationFlux = sublimation_d./(rho_ref.*As.*a0);
            wallMomFlux = accretionFlux.*V_d./a0;%.*L./dx_d; % [non-dimensional]
            wallLatentHeatFlux = wallMassFlux .* Lh_nd; % [non-dimensional]
            wallEnergyLoss = accretionFlux.*(T./(y-1) + y.*V.^2./2 + p./rho); % [non-dimensional]
            wallEnergyLoss2 = accretion_d.*(cv.*T_d + V_d.^2./2 + p_d./rho_d)./(rho_ref.*As.*a0.*Rg.*T_ref);%./L); % [non-dimensional]
            wallEnergyWin = sublimationFlux.*((T_w./T_ref)./(y-1) + p_eq_sg(T_w)./p_ref .* rho_ref./(p_eq_sg(T_w)./(Rg.*T_w))); % Energy from sublimation from wall at Tw, and saturated vapor.
            % Convective heat transfer coefficient and Reynolds, Nussel &
            % Prandlt numbers
            [h Re Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d);
            % Convection (non-dimensional:
            wallHeatFlux = (h.*(T_d-T_w).*c_d.*dx_d)./(T_ref.*Rg.*rho_ref.*a0.*As); 
            h = h./(Rg.*rho_ref.*a0);
            wallHeatFlux2 = (h.*(T_d-T_w)./T_ref.*c.*dx);% convection
            %fprintf('rd14=%f',r_d(150))
        end
        % Predictor step:
        % Compute derivatives of the solution vectors
        for i = 2 : n - 1
            J2(i) = 1 / y * rho(i) .* T(i) * fd_dx(A, dx, i);
            if rho_reduction == true
                dU1dt_p(i) = -fd_dx(F1, dx, i) - fd_dx(mdot_d, dx, i); % mass is reduced ~Deposition
            elseif wallInteractions == true
                dU1dt_p(i) = -fd_dx(F1, dx, i) - wallMassFlux(i)/dx; % mass is reduced ~Wall condensation 
            else
                dU1dt_p(i) = -fd_dx(F1, dx, i);
            end
            if wallInteractions == true
                if convectionFriction == true
                    dU2dt_p(i) = -fd_dx(F2, dx, i) + J2(i) - wallMomFlux(i)/dx - friction(i)/dx;% - friction(i)/dx; 
                else
                    dU2dt_p(i) = -fd_dx(F2, dx, i) + J2(i) - wallMomFlux(i)/dx; 
                end
            else
                dU2dt_p(i) = -fd_dx(F2, dx, i) + J2(i);
            end
            if nucleation == true
                if wallInteractions ==  true
                    if convectionFriction == true
                        dU3dt_p(i) = -fd_dx(F3, dx, i) + fd_dx(Edot_L, dx, i) - wallHeatFlux(i)/dx - wallEnergyLoss(i)/dx + wallEnergyWin(i)/dx;% - wallHeatFlux(i);% + wallLatentHeatFlux(i);% + wallHeatFlux(i); % forward difference
                        latentheat_fd(i) = fd_dx(Edot_L, dx, i);
                    else
                        dU3dt_p(i) = -fd_dx(F3, dx, i) + fd_dx(Edot_L, dx, i) - wallEnergyLoss(i)/dx + wallEnergyWin(i)/dx;%- wallHeatFlux(i)/dx  - wallHeatFlux(i);% + wallLatentHeatFlux(i);% + wallHeatFlux(i); % forward difference
                    	latentheat_fd(i) = fd_dx(Edot_L, dx, i);
                    end
                else
                    dU3dt_p(i) = -fd_dx(F3, dx, i) + fd_dx(Edot_L, dx, i);
                    latentheat_fd(i) = fd_dx(Edot_L, dx, i);
                end
            else
                if wallInteractions ==  true
                    if convectionFriction == true
                        dU3dt_p(i) = -fd_dx(F3, dx, i) - wallHeatFlux(i)/dx - wallEnergyLoss(i)/dx + wallenergyWin(i)/dx;% wallLatentHeatFlux(i)/dx;% + wallHeatFlux(i); % forward difference
                    else
                        dU3dt_p(i) = -fd_dx(F3, dx, i) - wallEnergyLoss(i)/dx + wallEnergyWin(i)/dx;%- wallHeatFlux(i)/dx  wallLatentHeatFlux(i)/dx;% + wallHeatFlux(i); % forward difference
                    end
                else
                    dU3dt_p(i) = -fd_dx(F3, dx, i); % forward difference
                end
            end
            
            % Shock Capturing equations, based on 'old' (current) U1-3 
            cc = C_x * abs(cdiff(p,i))/(p(i+1) + 2 * p(i) + p(i-1));
            S1(i) = cc * cdiff(U1,i);
            S2(i) = cc * cdiff(U2,i);
            S3(i) = cc * cdiff(U3,i);

            % PREDICTED solution vectors (including S-terms)
            U1(i)= U1(i) + dU1dt_p(i) * dt + S1(i);
            U2(i)= U2(i) + dU2dt_p(i) * dt + S2(i);
            U3(i)= U3(i) + dU3dt_p(i) * dt + S3(i);
        end
        
        % PREDICTED primitive variables and flux vectors
        rho = U1 ./ A;
        T = (y - 1) * (U3 ./ U1 - (y / 2) * (U2 ./ U1).^2);
        p = rho .* T;
        % Compute new M %
        M = V ./ sqrt(T);
        [rho_d T_d p_d V_d A_d x_d dx_d r_d] = dimensionalize(rho, T, p, M, y, Rg, A, x, dx, rho_ref, T_ref, p_ref, r_res_d, L, throat, r, As, circShape, crackLength);

        if nucleation == true
            % Dimensionalize for computation: f', dE, dR. ('_d' = dimensionalized) %
            Q = mean(rho_d .* A_d .* V_d); % average because; mass flow = C  
            dR_dx = dR_dx_MC(T_d,rho_d,M, b);
            [R dR_dx] = getParticleRadius(dR_dx, x_d, dx_d);
            [fp gn S] = solid_frac_MC(Q,dx_d,T_d,rho_d,R,dR_dx,A_d,x_d,V_d, n, Rg); % uses rho_d(i), T_d(i), M(i) 
            % Updating solid fraction, Change to Simpson's rule  %
            [f fp] = getSolidFraction(fp, x_d, dx_d);
            % Non-dimensionalize
            [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
            % Calculate energy of latent heat in all segments
            Edot_L = Edot_source(rho, f , fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, y); % Latent heat as Source term
        end
        if wallInteractions == true
            [wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, run_name,k, nt, false, m_stick, T_w, c_d, r_sub);
            wallMassFlux = wallMassFlux_d./(rho_ref.*As.*a0);
            [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
            [tau tau_d] = wallStress(rho, V, rho_d, V_d, T_d, r_d, rho_ref, a0); % Tangential stress at walls non-dimensional
            [friction] = getFriction(tau, r, L, A_d, throat, c, dx); % non-dimensional friction
            accretionFlux = accretion_d./(rho_ref.*As.*a0);
            sublimationFlux = sublimation_d./(rho_ref.*As.*a0);
            wallMomFlux = accretionFlux.*V_d./a0;%.*L./dx_d; % [non-dimensional]
            wallLatentHeatFlux = wallMassFlux .* Lh_nd; % [non-dimensional]
            wallEnergyLoss = accretionFlux.*(T./(y-1) + y.*V.^2./2 + p./rho); % [non-dimensional]
            wallEnergyLoss2 = accretion_d.*(cv.*T_d + V_d.^2./2 + p_d./rho_d)./(rho_ref.*As.*a0.*Rg.*T_ref);%./L); % [non-dimensional]
            wallEnergyWin = sublimationFlux.*((T_w./T_ref)./(y-1) + p_eq_sg(T_w)./p_ref .* rho_ref./(p_eq_sg(T_w)./(Rg.*T_w))); % Energy from sublimation from wall at Tw, and saturated vapor.
            [h Re Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d);
            wallHeatFlux = (h.*(T_d-T_w).*c_d.*dx_d)./(T_ref.*Rg.*rho_ref.*a0.*As); 
            h = h./(Rg.*rho_ref.*a0);
            wallHeatFlux2 = (h.*(T_d-T_w)./T_ref.*c.*dx);% convection
        end        
        if rho_reduction == true
            mdot_d  = mdot_sink(f, rho, A, V); % Mass reduction ~Deposition
        end
        F1 = U2;
        F2 =(U2.^2 ./ U1) + ((y - 1) / y) * (U3 - (y * U2.^2) ./ (2 * U1));
        F3 =(y * U2 .* U3 ./ U1) - (y * (y - 1) / 2) * (U2.^3 ./ U1.^2);
        % Corrector step
        for j = 2 : n - 1
            J2(j) = 1 / y * rho(j) .* T(j) * rd_dx(A, dx, j); 
            if rho_reduction == true
                dU1dt_c(j) = -rd_dx(F1, dx, j) - rd_dx(mdot_d, dx, j); % mass is reduced ~Deposition
            elseif wallInteractions == true
                dU1dt_c(j) = -rd_dx(F1, dx, j) - wallMassFlux(j)/dx; % mass is reduced ~Wall condensation 
            else
                dU1dt_c(j) = -rd_dx(F1, dx, j);
            end
            if wallInteractions == true
                if convectionFriction == true
                    dU2dt_c(j) = -rd_dx(F2, dx, j) + J2(j) - wallMomFlux(j)/dx- friction(j)/dx; 
                else
                    dU2dt_c(j) = -rd_dx(F2, dx, j) + J2(j) - wallMomFlux(j)/dx;
                end
            else
                dU2dt_c(j) = -rd_dx(F2, dx, j) + J2(j); 
            end
            if nucleation == true
                if wallInteractions ==  true
                    if convectionFriction == true
                        dU3dt_c(j) = -rd_dx(F3, dx, j) + rd_dx(Edot_L, dx, j) - wallEnergyLoss(j)/dx + wallEnergyWin(j)/dx - wallHeatFlux(j)/dx;% + wallLatentHeatFlux(j);% + wallHeatFlux(j); % forward difference
                        latentheat_rd(j) = rd_dx(Edot_L, dx, j);
                    else
                        dU3dt_c(j) = -rd_dx(F3, dx, j) + rd_dx(Edot_L, dx, j) - wallEnergyLoss(j)/dx + wallEnergyWin(j)/dx;% + wallLatentHeatFlux(j);% + wallHeatFlux(j); % forward difference
                        latentheat_rd(j) = rd_dx(Edot_L, dx, j);
                    end
                elseif wallInteractions ==  false
                    dU3dt_c(j) = -rd_dx(F3, dx, j) + rd_dx(Edot_L, dx, j);
                    latentheat_rd(j) = rd_dx(Edot_L, dx, j);
                end
            else
                if wallInteractions ==  true
                    if convectionFriction == true
                        dU3dt_c(j) = -rd_dx(F3, dx, j) - wallEnergyLoss(j)/dx + wallEnergyWin(j)/dx - wallHeatFlux(j)/dx;%- wallLatentHeatFlux(j)/dx;% + wallHeatFlux(j); % forward difference
                    else
                        dU3dt_c(j) = -rd_dx(F3, dx, j) - wallEnergyLoss(j)/dx + wallEnergyWin(j)/dx;%- wallLatentHeatFlux(j)/dx;% + wallHeatFlux(j); % forward difference
                    end
                else
                    dU3dt_c(j) = -rd_dx(F3, dx, j); % rearward difference
                end
            end
                        
            % Shock Capturing equations
            cc = C_x * abs(cdiff(p,j))/(p(j+1) + 2 * p(j) + p(j-1));
            S1(j) = cc * cdiff(U1,j);
            S2(j) = cc * cdiff(U2,j);
            S3(j) = cc * cdiff(U3,j);
        end
        
        % Averaging corrector step and predictor step
        dU1dt_a = 0.5 * (dU1dt_p + dU1dt_c);
        dU2dt_a = 0.5 * (dU2dt_p + dU2dt_c);
        dU3dt_a = 0.5 * (dU3dt_p + dU3dt_c);

        % Updating solution vectors after averaging
        for m = 2 : n - 1 
           U1(m) = U1_old(m) + dU1dt_a(m) * dt + S1(m);
           U2(m) = U2_old(m) + dU2dt_a(m) * dt + S2(m);
           U3(m) = U3_old(m) + dU3dt_a(m) * dt + S3(m);
        end   
        resU1(k) = sum(abs(U1 - U1_old))/sum(U1_old);
        resU2(k) = sum(abs(U2 - U2_old))/sum(U2_old);
        resU3(k) = sum(abs(U3 - U3_old))/sum(U3_old);

        % Boundary conditions 
        % Inlet, fix rho & T 
        U1(1) = rho(1)*A(1);
        if V0_fix == true
            U2(1) = rho(1)*A(1)*V(1);
        else
            U2(1) = 2*U2(2) - U2(3);
            V(1) = U2(1)/U1(1);
        end
       
        T(1) = 1; 
        U3(1) = U1(1) * (T(1) / (y - 1) + y / 2 * V(1)^2);

        % Outlet
        U1(n) = 2*U1(n-1) - U1(n-2);
        U2(n) = 2*U2(n-1) - U2(n-2);
        if p_e > 0 % fix for shock-capturing
          U3(n) = p_e * A(n)/(y - 1) + y / 2 * U2(n) * (U2(n) / U1(n));%V(n);   
          p(n) = pe;
        else % if pe = undefined [subsonic-supersonic],  
          U3(n) = 2*U3(n-1) - U3(n-2); %Only when pe/p0 is NOT set 
        end
        % Updation of flow field variables
        rho = U1 ./ A;
        V = U2 ./ U1;
        T = (y - 1) * (U3 ./ U1 - (y / 2) * (U2 ./ U1).^2);
        T(1) = 1;
        p = rho .* T;
        
        % Defining mass flow and mach number
        m = rho .* A .* V;
        m_d = rho_d .* A_d .* V_d;
        M = V ./ sqrt(T);
        % Dimensional for computation: f', dE, dR. ('_d' = dimensionalized) %
        [rho_d T_d p_d V_d A_d x_d dx_d r_d] = dimensionalize(rho, T, p, M, y, Rg, A, x, dx, rho_ref, T_ref, p_ref, r_res_d, L, throat, r, As, circShape, crackLength);

        % Capturing the values at inlet, throat and exit for each iterations
        mach_in(k) = M(1);
        pressure_in(k) = p(1);
        density_in(k) = rho(1);
        temperature_in(k) = T(1);
        
        mach_thr(k) = M(throat);
        pressure_thr(k) = p(throat);
        density_thr(k) = rho(throat);
        temperature_thr(k) = T(throat);
        
        mach_ex(k) = M(n);
        pressure_ex(k) = p(n);
        density_ex(k) = rho(n);
        temperature_ex(k) = T(n);
        if nucleation == true
            gn_max = max(gn);
            gn_max_n = find(gn_max == gn);
            gn_max_z(k) = x_d(gn_max_n);
            f_exit(k) = f(end);
        end
        if convergenceCriteria == true && resU1(k) < critRes && resU2(k) < critRes  && resU3(k) < critRes
            break % Exit for-loop when convergence criteria is met
        end
        if isreal(M) == 0 || any(T<0) == true || sum(isnan(M)) > 0 || sum(isnan(T)) > 0
            display(['Imaginary numbers at ', num2str(k), ' iterations.'])
            break
        end
        if stop == true % type in command window "stop = true" to exit loop
            break
        end
        
        %fprintf('rd1=%f',r_d(150))
        r_old =r_d;
        r_d = r_old-E/1000*100000000*dt;
        r=r_d;    
        %fprintf('rd2=%f',r_d(150))
        %fprintf('El=%f',E(150))
        %fprintf('dt=%f',dt)
        %fprintf('k=%f',k)
        
        
    end
    [h Re Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d);
    numbers = figure(999);
    subplot(2,1,1)
    semilogy(x_d,Re)
    hold on
    semilogy(x_d,Nu)
    semilogy(x_d,Pr)
    legend('Re','Nu','Pr')
    grid on
    xlabel('x')
    ylabel('[-]')
    subplot(2,1,2)
    plot(x_d,h)
    ylabel('Convective heat transfer coefficient [J/(kgK)]')
    xlabel('x')
    name = ['numbers.png'];
    saveas(numbers, ['Figures/' run_name '/' name]);
    
    % Equilibrium density and super-saturation
    rho_eq = p_eq_lg(T_d)./(Rg*T_d);
    rho_eq_2 = p_eq_sg(T_d)./(Rg.*T_d);
    S = rho_d./rho_eq;
    S_2 = rho_d./rho_eq_2;
    SuperSaturation = figure(4);
    subplot(3,1,1)
    semilogy(x_d,rho_eq)
    hold on 
    semilogy(x_d, rho_eq_2)
    ylabel('\rho_{eq} [kg/m^3]')
    legend('solid/liquid','solid/gas')
    xline(x_d_throat,'k--')
    subplot(3,1,2)
    plot(x_d, rho_d)
    hold on
    xline(x_d_throat,'k--')
    ylabel('\rho [kg/m^3]')
    subplot(3,1,3)
    plot(x_d,S)
    hold on
    plot(x_d,S_2)
    xline(x_d_throat,'k--')
    ylabel('Super-Saturation [-]')
    xlabel('z [m]')
    name = ['SuperSaturation.png'];
    saveas(SuperSaturation, ['Figures/' run_name '/' name]);
    % Flux vectors and source vectors
    figure(5)
    subplot(4,1,1)
    plot(x_d,F1)
    hold on
    xline(x_d_throat,'k--')
    ylabel('F1')
    subplot(4,1,2)
    plot(x_d,F2)
    hold on
    xline(x_d_throat,'k--')
    ylabel('F2')
    subplot(4,1,3)
    plot(x_d,F3)
    hold on
    xline(x_d_throat,'k--')
    ylabel('F3')
    subplot(4,1,4)
    plot(x_d(2:length(J2)),J2(2:end))
    hold on
    xline(x_d_throat,'k--')
    ylabel('J2')
    xlabel('z [m]')
    % Solution vectors
    SolutionVectors = figure(6);
    subplot(3,1,1)
    plot(x_d,U1)
    hold on
    xline(x_d_throat,'k--')
    ylabel('U1')
    subplot(3,1,2)
    plot(x_d,U2)
    hold on
    xline(x_d_throat,'k--')
    ylabel('U2')
    subplot(3,1,3)
    plot(x_d,U3)
    hold on
    xline(x_d_throat,'k--')
    ylabel('U3')
    xlabel('z [m]')
    name = ['SolutionVectors.png'];
    saveas(SolutionVectors, ['Figures/' run_name '/' name]);
    
    if verification == true
        [dum1, dum2, rho_verification, T_verification, V_verification, c_verification, R_verification, f_verification, gn_verification, S_verification, dum3, dum4, dum5, dum6, dum7, dum8, M_verification, p_verification, x_d_Verification, dx_d_Verification] = getData_Verification(channelNr);
        [dum1, dum2, rho_offset, T_offset, V_offset, dum4, R_offset, f_offset, gn_offset, S_offset, dum9, dum10, dum11, dum12, dum13, dum14, M_offset, p_offset] = getData_withoffset(channelNr, Rg);
    end
    % Artificial viscosity magnitude
    ArtificialViscosity = figure(7);
    subplot(3,1,1)
    plot(x_d(1:n-1),S1)
    hold on
    xline(x_d_throat,'k--')
    ylabel('S1')
    subplot(3,1,2)
    plot(x_d(1:n-1),S2)
    hold on
    xline(x_d_throat,'k--')
    ylabel('S2')
    subplot(3,1,3)
    plot(x_d(1:n-1),S3)
    hold on
    xline(x_d_throat,'k--')
    ylabel('S3')
    xlabel('z [m]')
    % Dimensional results plot
    dim_Res = figure(9);
    subplot(411)
    plot(x_d,V_d,'b','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,V_verification,'--r','linewidth',2)
        if offset == true
            plot(x_d,V_offset,'--y','linewidth',1.5) 
        end
    end
    ylabel('V [m/s]')
    xlabel('z [m]')
    grid minor
    hold off
    subplot(412)
    plot(x_d,T_d,'g','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,T_verification,'--r','linewidth',2)
        if offset == true
            plot(x_d,T_offset,'--y','linewidth',1.5) 
        end
    end
    hold off
    ylabel('T [K]')
    xlabel('z [m]')
    grid minor
    subplot(413)
    plot(x_d,rho_d,'r','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,rho_verification,'--b','linewidth',2)
        if offset == true
            plot(x_d,rho_offset,'--y','linewidth',1.5) 
        end
    end
    hold off
    ylabel('\rho [kg/m^3]')
    grid minor
    xlabel('z [m]')
    subplot(414)
    plot(x_d,p_d,'m','linewidth',2)
    hold on
    xline(x_d_throat,'k--')
    if verification == true
        hold on
        plot(x_d_Verification,p_verification,'--b','linewidth',2)
        if offset == true
            plot(x_d,p_offset,'--y','linewidth',1.5) 
        end
    end
    ylabel('p [pa]')
    grid minor
    hold off
    xlabel('z [m]')
    name = ['DimensionalResults.png'];
    saveas(dim_Res, ['Figures/' run_name '/' name]);
    w = figure(11);
    clf;
    semilogy(resU1,'linewidth',2), hold on
    semilogy(resU2,'linewidth',2), semilogy(resU3,'linewidth',2)
    grid minor
    ylabel('Residuals [-]')
    xlabel('Iterations')
    title('Residuals of Solution Vectors (Conservative)', 'FontAngle', 'italic', 'FontSize', 13, 'FontWeight', 'bold')
    legend('U1','U2','U3');
    name = ['residualsCon.png'];
    saveas(w, ['Figures/' run_name '/' name]);
        dUdt_a = figure(18);
    subplot(3,1,1)
    plot(x_d(1:n-1),dU1dt_a)
    hold on
    xline(x_d_throat,'k--')
    ylabel('dU1dt_a')
    subplot(3,1,2)
    plot(x_d(1:n-1),dU2dt_a)
    hold on
    xline(x_d_throat,'k--')
    ylabel('dU2dt_a')
    subplot(3,1,3)
    plot(x_d(1:n-1),dU3dt_a)
    hold on
    xline(x_d_throat,'k--')
    ylabel('dU3dt_a')
    name = ['dUdt_a.png'];
    saveas(dUdt_a, ['Figures/' run_name '/' name]);
    MassFluxVariables = figure(19);
    subplot(3,1,1)
    plot(x_d,T_w)
    hold on
    plot(x_d, T_d)
    xline(x_d_throat,'k--')
    legend('wall','flow')
    ylabel('T [K]')
    subplot(3,1,2)
    plot(x_d, p_d)
    xline(x_d_throat,'k--')
    ylabel('P [pa]')
    subplot(3,1,3)
    plot(x_d, (p_d./sqrt(2.*pi.*Rg.*T_d)))
    xline(x_d_throat,'k--')
    ylabel('E terms')
    xlabel('z [m]')
    name = ['MassFluxVariables.png'];
    saveas(MassFluxVariables, ['Figures/' run_name '/' name]);
    
    % nucleation computations for figures
    if nucleation == true
        % Dimensionalize for computation: f', dE, dR. ('_d' = dimensionalized) %
        [rho_d T_d p_d V_d A_d x_d dx_d] = dimensionalize(rho, T, p, M, y, Rg, A, x, dx, rho_ref, T_ref, p_ref, r_res_d, L, throat, r, As, circShape, crackLength);
        Q = mean(rho_d .* A_d .* V_d); % average because; mass flow = C  
        [dR_dx rho_eq] = dR_dx_MC(T_d,rho_d,M, b);
        [R dR_dx] = getParticleRadius(dR_dx, x_d, dx_d);
        [fp gn S] = solid_frac_MC(Q,dx_d,T_d,rho_d,R,dR_dx,A_d,x_d,V_d, n, Rg); % uses rho_d(i), T_d(i), M(i) 
        % Updating solid fraction, Change to Simpson's rule  %
        [f fp] = getSolidFraction(fp, x_d, dx_d);
        % Non-dimensionalize
        [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
        % Calculate energy of latent heat in all segments
        Edot_L = Edot_source(rho, f , fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, y); % Latent heat as Source term
        [r_e P r_peak] = getParticleSizeDistribution(R, x_d, gn, dx_d, V_d, A_d);

        u = figure(3);
        subplot(611)
        %hold on 
        plot(mach_in, 'b', 'linewidth', 1.5)
        hold on
        plot(mach_thr, 'g', 'linewidth', 1.2)
        plot(mach_ex, 'r', 'linewidth', 1.2)
        hold off
        title('Variation of Nozzle Inlet, Throat and Exit Conditions for Conservative Form', 'FontAngle', 'italic', 'FontSize', 13, 'FontWeight', 'bold')
        ylabel('M', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor
        legend('Inlet','Throat','Exit','Location','best')
        legend('boxoff')
        subplot(612)
        %hold on 
        plot(pressure_in, 'b', 'linewidth', 1.5)
        hold on
        plot(pressure_thr, 'g', 'linewidth', 1.2)
        plot(pressure_ex, 'r', 'linewidth', 1.2)
        hold off
        ylabel('p/p_{0}', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor

        subplot(613)
        plot(density_in, 'b', 'linewidth', 1.5)
        hold on
        plot(density_thr, 'g', 'linewidth', 1.2)
        plot(density_ex, 'r', 'linewidth', 1.2)
        hold off
        ylabel('\rho/\rho_{0}', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor

        subplot(614)
        plot(temperature_in, 'b', 'linewidth', 1.5)
        hold on
        plot(temperature_thr, 'g', 'linewidth', 1.2)
        plot(temperature_ex, 'r', 'linewidth', 1.2)
        hold off
        ylabel('T/T_{o}', 'FontSize', 11, 'FontWeight', 'bold')
        xlabel('Number of Iterations', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor

        subplot(615)
        plot(f_exit, 'b', 'linewidth', 1.5)
        hold off
        ylabel('f_{exit}', 'FontSize', 11, 'FontWeight', 'bold')
        xlabel('Number of Iterations', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor

        subplot(616)
        plot(gn_max_z, 'b', 'linewidth', 1.5)
        hold off
        ylabel('\gamma_{nuc} peak location', 'FontSize', 11, 'FontWeight', 'bold')
        xlabel('Number of Iterations', 'FontSize', 11, 'FontWeight', 'bold')
        grid minor
        name = ['nozzle_throat.png'];
        saveas(u, ['Figures/' run_name '/' name]);
        
        SuperSaturation_Nucleationrate = figure(8);
        subplot(5,1,1)
        plot(x_d,S)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(S) max(S)],'k--')
        if verification == true
            hold on
            plot(x_d_Verification,S_verification,'--r','linewidth',2)
            if offset == true
                plot(x_d,S_offset,'--y','linewidth',1.5)
            end
        end
        hold off
        ylabel('Super-saturation')
        subplot(5,1,2)
        semilogy(x_d, gn)
        hold on
        xline(x_d_throat,'k--')%semilogy(x_d_throat,[min(gn) max(gn)],'k--')
        if verification == true
            hold on
            plot(x_d_Verification,gn_verification,'--r','linewidth',2)
            if offset == true
                plot(x_d,gn_offset,'--y','linewidth',1.5) 
            end
        end
        ylabel('log10(y_{nuc}[(m^{-3}s^{-1})])')
        ylim([1 10^20])
        xlim([0 max(x_d)])
        hold off
        subplot(5,1,3)
        plot(x_d, T_d)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(T_d) max(T_d)],'k--')
        if verification == true
            hold on
            plot(x_d_Verification,T_verification,'--r','linewidth',2)
            if offset == true
                plot(x_d,T_offset,'--y','linewidth',1.5) 
            end
        end
        ylabel('T [K]')
        hold off
        subplot(5,1,4)
        plot(x_d(1:length(latentheat_fd)),latentheat_fd)
        hold on
        plot(x_d(1:length(latentheat_rd)),latentheat_rd)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(latentheat_fd) max(latentheat_fd)],'k-')
        legend('forward','rearward')
        ylabel('L_h [J/s]/dx')
        xlabel('x (dimensional)')
        subplot(5,1,5)
        plot(x_d,Edot_L)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(Edot_L) max(Edot_L)],'k--')
        ylabel('Edot_L')
        xlabel('x (dimensional)')
        name = ['SuperSaturation_Nucleationrate.png'];
        saveas(SuperSaturation_Nucleationrate, ['Figures/' run_name '/' name]);
        particleSizeDistribution = figure(20);
        semilogx(r_e*10^6, P)
        title('Particle Size Distribution')
        ylabel('Particle concentration [n/m^3]')
        xlabel('R(\mu m)_{exit}')
        name = ['particleSizeDistribution.png'];
        saveas(particleSizeDistribution, ['Figures/' run_name '/' name]);

        R_fp_f_gn = figure(10);
        subplot(5,1,1)
        semilogy(x_d,dR_dx)
        ylim([ 10^-9 max(dR_dx)*10]) 
        %ylim([ 10^-4 max(dR_dx)*10])
        hold on
        xline(x_d_throat,'k--')%semilogy(x_d_throat,[min(dR_dx) max(dR_dx)],'k--')
        xlabel('z [m]')
        ylabel('dR/dx') 
        grid on
        subplot(5,1,2)
        plot(x_d,R)
        hold on
        xline(x_d_throat,'k--')
        if verification == true
            hold on
            plot(x_d_Verification,R_verification,'--b','linewidth',2)
            if offset == true
                plot(x_d,R_offset,'--y','linewidth',1.5) 
            end
        end
        xlabel('z [m]')
        ylabel('R')
        hold off
        %ylim([0 max(R)*1.1]) 
        grid on
        subplot(5,1,3)
        semilogy(x_d,fp)
        ylim([10^-100 max(fp)*10])
        hold on
        xline(x_d_throat,'k--')%semilogy(x_d_throat,[min(fp) max(fp)],'k--')
        xlabel('z [m]')
        ylabel('fp')
        grid minor
        subplot(5,1,4)
        plot(x_d,f)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(f) max(f)],'k--')
        if verification == true
            hold on
            plot(x_d_Verification,f_verification,'--b','linewidth',2)
            if offset == true
                plot(x_d,f_offset,'--y','linewidth',1.5) 
            end
        end
        hold off
        xlabel('z [m]')
        ylabel('f')
        %ylim([0 1.1*max(f)])
        grid minor
        subplot(5,1,5)
        semilogy(x_d, gn)
        hold on
        xline(x_d_throat,'k--')%plot(x_d_throat,[min(gn) max(gn)],'k--')
        if verification == true
            hold on
            plot(x_d_Verification,gn_verification,'--b','linewidth',2)
            if offset == true
                plot(x_d,gn_offset,'--y','linewidth',1.5)
            end
        end
        xlabel('z [m]')
        ylabel('Nucleation Rate [m^{-3}s^{-1}]')
        ylim([10^0 max(gn)*10])
        grid minor
        hold off
        name = ['R_fp_f_gn.png'];
        saveas(R_fp_f_gn, ['Figures/' run_name '/' name]);
    end
    
    % Wall interactions computations for figures
    m_d = rho_d .* A_d .* V_d; % mass flow along z
    m_d_e = m_d(end); % Exit mass flow
    [wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, run_name, k, nt, false, m_stick, T_w, c_d, r_sub);
    wallMassFlux = wallMassFlux_d./(rho_ref.*As.*a0);%E = linspace(EE,EE,n);
    [fp_nd Lh_nd Rg_nd E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, y, As);   
    [tau tau_d] = wallStress(rho, V, rho_d, V_d, T_d, r_d, rho_ref, a0); % Tangential stress at walls non-dimensional
    [friction] = getFriction(tau, r, L, A_d, throat, c, dx); % non-dimensional friction
    accretionFlux = accretion_d./(rho_ref.*As.*a0);
    sublimationFlux = sublimation_d./(rho_ref.*As.*a0);
    wallMomFlux = accretionFlux.*V_d./a0;%.*L./dx_d; % [non-dimensional]
    wallLatentHeatFlux = wallMassFlux .* Lh_nd; % [non-dimensional]
    wallEnergyLoss = accretionFlux.*(T./(y-1) + y.*V.^2./2 + p./rho); % [non-dimensional]
    wallEnergyLoss2 = accretion_d.*(cv.*T_d + V_d.^2./2 + p_d./rho_d)./(rho_ref.*As.*a0.*Rg.*T_ref);%./L); % [non-dimensional]
    wallEnergyWin = sublimationFlux.*((T_w./T_ref)./(y-1) + p_eq_sg(T_w)./p_ref .* rho_ref./(p_eq_sg(T_w)./(Rg.*T_w))); % Energy from sublimation from wall at Tw, and saturated vapor.
    wallEnergyWin2 = sublimation_d.*(cv.*T_w + p_eq_sg(T_w)./(p_eq_sg(T_w)./(Rg.*T_w)))./(rho_ref.*As.*a0.*Rg.*T_ref);%./L); % [non-dimensional]
    
    % convective heat transfer coefficient, Reynolds, Nusselt % Prandtl
    % number
    [h Re Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d);
    wallHeatFlux = (h.*(T_d-T_w).*c_d.*dx_d)./(T_ref.*Rg.*rho_ref.*a0.*As); 
    wallHeatFlux_d = (h.*(T_d-T_w).*c_d.*dx_d);
    h = h./(Rg.*rho_ref.*a0);
    wallHeatFlux2 = (h.*(T_d-T_w)./T_ref.*c.*dx);% convection
    tau_d = tau.*rho_ref.*a0.^2;%.*As./L;
    [friction_d] = getFriction(tau_d, r_d, L, A_d, throat, c_d, dx_d); % non-dimensional friction
    %wallMassFlux_d = (E.*c_d);%.*(L/sqrt(A_d(throat)*y)); % [non-dimensional]
    wallMomFlux_d = accretion_d.*V_d; % [non-dimensional]
    wallLatentHeatFlux_d = wallMassFlux_d .* L_h; % [non-dimensional]
    q_d = 100; % 5.000 - 100.000 W/(m^2K)
    
    wallInteractions_nd = figure(145);
    title('Non-dimensional variables')
    subplot(2,5,1)
    plot(x_d,wallHeatFlux)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallHeatFlux')
    subplot(2,5,2)
    semilogy(x_d,E_nd)
    ylabel('E_{nd}')
    subplot(2,5,3)
    plot(x_d,wallMassFlux)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallMassFlux')
    subplot(2,5,4)
    plot(x_d,wallMomFlux)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallMomFlux')
    subplot(2,5,5)
    plot(x_d,wallLatentHeatFlux)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallLatentHeatFlux')
    subplot(2,5,6)
    plot(x_d,wallHeatFlux2)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallHeatFlux 2')
    subplot(2,5,7)
    plot(x_d, wallEnergyLoss)
    hold on
    xline(x_d_throat,'k--')
    ylabel('Wall Energy Loss')
    subplot(2,5,8)
    plot(x_d, wallEnergyWin)
    hold on
    xline(x_d_throat,'k--')
    ylabel('Wall Energy Win')
    subplot(2,5,9)
    plot(x_d, friction)
    hold on
    xline(x_d_throat,'k--')
    ylabel('friction')
    subplot(2,5,10)
    plot(x_d, c_stick)
    hold on
    xline(x_d_throat,'k--')
    ylabel('Sticking coefficient')
    name = ['wallInteractions_nd.png'];
    saveas(wallInteractions_nd, ['Figures/' run_name '/' name]);
    wallInteractions_d = figure(1045);
    title('Dimensional variables')
    subplot(2,4,1)
    semilogy(x_d,E)
    hold on
    xline(x_d_throat,'k--')
    ylabel('E [kg/(m^2s)]')
    subplot(2,4,2)
%     plot(x_d,wallMassFlux_d)
%     hold on
%     plot(x_d, accretion_d)
%     plot(x_d, sublimation_d)
    semilogy(x_d,wallMassFlux_d)
    hold on
    semilogy(x_d, accretion_d)
    semilogy(x_d, sublimation_d)    
    xline(x_d_throat,'k--')
    ylabel('[kg/s]')
    legend('Wall mass flux','accretion','sublimation') 
    wallMassFlux_pm_d = wallMassFlux_d./dx_d;
    accretedMass(1) = 0;
    accretedFraction(1) = 0;
    for i = 2:length(wallMassFlux_d)
        accretedMass(i) = trapz(x_d(1:i),wallMassFlux_pm_d(1:i));
        %accretedFraction(i) = trapz(x_d(1:i),m_sigma(1:i));
        
    end
    accretedFraction = accretedMass./(accretedMass + m_d); % Works well when flow is converged 
    subplot(2,4,3)
    plot(x_d, accretedMass)
    hold on
    xline(x_d_throat, 'k--')
    ylabel('Accreted mass [kg/s]')
    xlabel('x [m]')
    subplot(2,4,4)
    plot(x_d, accretedFraction)
    hold on
    xline(x_d_throat, 'k--')
    ylabel('Accreted fraction [-]')
    xlabel('x [m]')
    subplot(2,4,5) 
    plot(x_d,wallHeatFlux_d)
    hold on
    xline(x_d_throat,'k--')
    ylabel('Convective Heat Flux')
    subplot(2,4,6)
    plot(x_d,wallMomFlux_d)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallMomFlux')
    subplot(2,4,7)
    plot(x_d,wallLatentHeatFlux_d)
    hold on
    xline(x_d_throat,'k--')
    ylabel('wallLatentHeatFlux')
    subplot(2,4,8)
    plot(x_d, friction_d)
    hold on
    xline(x_d_throat,'k--')
    ylabel('friction')
    name = ['wallInteractions_d.png'];
    saveas(wallInteractions_d, ['Figures/' run_name '/' name]);
    
    viscosity = 0.00013; % water vapor
    Re = rho_d.*2.*r_d.*V_d./viscosity; % Reynolds Number
    f_f = 16./Re; % Friction coefficient
     
    
    if rho_reduction == true
        mdot_d  = mdot_sink(f, rho, A, V); % Mass reduction ~Deposition
    end
    % Show some results %
    display('Final Results')
    display(['n = ', num2str(n)])
    display(['\rho = ', num2str(rho(end))])
    display(['V = ', num2str(V(end))])
    display(['V_{initial} = ', num2str(V(1))])
    display(['T = ', num2str(T(end))])
    display(['p = ', num2str(p(end))])
    display(['M0 = ', num2str(M(1))])
    display(['M_throat = ', num2str(mach_thr(end))])
    display(['M1 = ', num2str(M(end))])
    display(['massflow = ', num2str(U2(end))])
    display(['U1 = ', num2str(U1(end))])   
    display(['U2 = ', num2str(U2(end))])
    display(['U3 = ', num2str(U3(end))])
    if nucleation == true
        display(['f = ', num2str(f(end))])
        display(['y_{nuc} max = ', num2str(round(max(gn)))])
    end
    display(['T = ', num2str(T_d(end))])
    display(['p_e/p_0 = ', num2str(round(p_e,3))])
    
    % Knudsen Number %
    kb = 1.38064852*10^(-23); %Boltzmann
    dH2O = 3.533095*10^(-10);
    Lmax = r_res_d; % characteristic length (reservoir diameter)
    Lmin = 0.008; % throat diameter
    Kn = (kb .* T_d) ./ (sqrt(2) .* pi .* dH2O .^2 .* p_d .* Lmax);
    Kn2 = (kb .* T_d) ./ (sqrt(2) .* pi .* dH2O .^2 .* p_d .* Lmin);
    
    % Save all data
    if testrun == false 
        % Generate Excel File %
        %fileName = 'runData.xlsx';
        fileName = 'runData.txt';
        folder = ['Figures/' run_name];                  % Folder Name (for figures)
        fullFileName = fullfile(folder, fileName); 
        if nucleation == true
            T_result = T_d(end); T_min = min(T_d); T_max = max(T_d); rho_result = rho_d(end); p_result = p_d(end);V_result = V_d(end);
            M_result = M(end); f_result = f(end); R_result = R(end); S_max = max(S); gn_max = max(gn); Edot_L_max = max(Edot_L);
            matrix = [T_result T_min T_max rho_result p_result V_result M_result f_result R_result S_max gn_max Edot_L_max];        
        else
            T_result = T_d(end); T_min = min(T_d); T_max = max(T_d); rho_result = rho_d(end); p_result = p_d(end);V_result = V_d(end);
            M_result = M(end); S_max = max(S); 
            matrix = [T_result T_min T_max rho_result p_result V_result M_result  S_max ];        
        end
        writematrix(matrix, fullFileName);

%         saveTable = table(T_result, T_min, T_max, rho_result, p_result, V_result, M_result, f_result, R_result, S_max, gn_max, Edot_L_max, m_d_e, Q, r_peak);
%         %saveTable.Properties.VariableUnits = {'-' ,'m/s','-'}; 
%         writetable(saveTable, fullFileName, 'Sheet',1);
    end
    % Save All Data %
    if saveAllData == true 
        %fileName = 'runFullData.xlsx';
        fileName = 'runFullData.txt';
        folder = ['Figures/' run_name];                  % Folder Name (for figures)
        fullFileName = fullfile(folder, fileName); 
        %matrix = ['f' str2num(f);'V' str2num(V_d); 'M' str2num(M)];
        if nucleation == true
            matrix = [T_d; rho_d; p_d; V_d; M; f; R; S; gn; fp; dR_dx; Edot_L; E; A_d; x_d; r_d; T_w; r_e; P; m_d;];
            matrix = transpose(matrix);
            matrix2 = ['T',{'rho'},'p','V','M','f','R','S','y_nuc','fp','dR_dx','Edot_L', 'E','A', 'x', 'r','T_w','r_e','P','m_d'];
            matrix2 = ["T","rho","p","V","M","f","R","S","y_nuc","fp","dR_dx","Edot_L", "E","A", "x", "r","T_w","r_e","P","m_d"];
       else
            matrix = [T_d; rho_d; p_d; V_d; M; f; R; S; gn; fp; dR_dx; E; A_d; x_d; r_d; T_w; m_d;];
            matrix = transpose(matrix);
            matrix2 = ['T',{'rho'},'p','V','M','f','R','S','y_nuc','fp','dR_dx', 'E','A', 'x', 'r','T_w','m_d'];
            matrix2 = ["T","rho","p","V","M","f","R","S","y_nuc","fp","dR_dx", "E","A", "x", "r","T_w","m_d"];
       end
        writematrix(matrix, fullFileName);
        fileName = 'constantsData.txt';
        folder = ['Figures/' run_name];                  % Folder Name (for figures)
        fullFileName = fullfile(folder, fileName); 
        %matrix = ['f' str2num(f);'V' str2num(V_d); 'M' str2num(M)];
        matrix = [L; L_h; cv; T_ref; Rg; rho_ref; y; As; crackLength; m_stick; r_sub; x_d_throat];
        matrix = transpose(matrix);
        matrix2 = ['L', {'L_h'},'cv','T_ref','Rg','rho_ref','y','As','crackLength','m_stick','r_sub','x_d_throat'];
        writematrix(matrix, fullFileName);
        fileName = 'Residuals.txt';
        folder = ['Figures/' run_name];                  % Folder Name (for figures)
        fullFileName = fullfile(folder, fileName); 
        matrix = [resU1; resU2; resU3];
        matrix = transpose(matrix);
        matrix3 = ["resU1","ResU2","ResU3"];
        writematrix(matrix, fullFileName); % Did for 'the cluster' can't read xlsx.
    end
    

    if length(Cx) > 1
        ccc = figure(21);
        subplot(211)
        hold on
        if rem(l,2)==0
            plot(x, M,'--','linewidth',1.5)
        else
            plot(x, M,'linewidth',1.5)
        end
        if l == length(Cx)
            %title('Parameters of Quasi - 1D Nozzle Flow', 'FontSize', 13, 'FontWeight', 'bold', 'FontAngle', 'italic')
            ylabel('M', 'FontWeight', 'bold', 'FontSize', 11)
            grid minor
        end

        subplot(212)
        hold on
        if rem(l,2)==0
            plot(x, p,'--','linewidth',1.5)
        else
            plot(x, p,'linewidth',1.5)
        end
        legend_text = ['Cx = ', num2str(C_x)];
        Legend{l}= legend_text;
        if l == length(Cx)
            ylabel('P/P_0', 'FontWeight', 'bold', 'FontSize', 11)
            xlabel('Non-Dimensional Length of the Nozzle', 'FontWeight', 'bold', 'FontSize', 11)
            grid minor
            name = ['P_M_Behaviour_Cx_variation.png'];
            legend(Legend,'Location','best');
            legend('boxoff');
            %magnifyOnFigure;
            saveas(ccc, ['Figures/' run_name '/' name]);
        end
    end
    
    if length(C) > 1
        CCC = figure(22);
        subplot(1,2,1)
        hold on
        plot(x_d, V_d, 'linewidth', 1.3)
        grid minor
        legend_text = ['C = ', num2str(C)];
        Legend{l}= legend_text;
        if l == length(C)
            %ylabel('P/P_0', 'FontWeight', 'bold', 'FontSize', 11)
            ylabel('V [m/s]')
            xlabel('z [m]');%, 'FontWeight', 'bold', 'FontSize', 11)
            name = ['Res_V_Behaviour_C_variation.png'];
            legend(Legend,'Location','best');
            legend('boxoff');
        end 
        subplot(1,2,2)
        hold on
        semilogy(resU1,'linewidth',2), hold on
        semilogy(resU2,'linewidth',2), semilogy(resU3,'linewidth',2)
        grid minor
        ylabel('Residuals')
        xlabel('Iterations')
        if l == length(C)
            saveas(CCC, ['Figures/' run_name '/' name]);
        end 
    end

    
    if length(pe) > 1
        ppp = figure(23);
        subplot(221)
        hold on
        plot(x, p_d, 'linewidth', 2)
        grid minor
        legend_text = ['p_e = ', num2str(round(p_e*p_ref)), ' Pa'];
        Legend{l}= legend_text;


        if l == length(pe)
            %ylabel('P/P_0', 'FontWeight', 'bold', 'FontSize', 11)
            ylabel('p [Pa]')
            xlabel('x/L');%, 'FontWeight', 'bold', 'FontSize', 11)
            name = ['P_M_Behaviour_pe_variation.png'];
            legend(Legend,'Location','best');
            legend('boxoff');
        end        
        subplot(223)
        hold on
        plot(x, M, 'linewidth', 2)
        grid minor
        ylabel('M')
        subplot(222)
        hold on
        plot(x, T_d, 'linewidth', 2)
        grid minor
        ylabel('T [K]')
        subplot(224)
        hold on
        plot(x, rho_d, 'linewidth', 2)
    
        grid minor
        ylabel('\rho [kg/m^3]')
        if l == length(pe)
            saveas(ppp, ['Figures/' run_name '/' name]);
        end    
    end 
    

end


 