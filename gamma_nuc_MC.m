function [g_nuc S] = gamma_nuc(T,rho)
% Computes the nucleation rate for icy grains (see Schmidt et al)
R = 461.52; %Gas constant water vapour [J/(kg K)]
T_range = [259 253.7 248.5 244.1 238.8 233.5 228.2 222.6 217.1];
ds = 0;
s_max = 10.^[(0.89 + ds)  (0.96 + ds) (1.0 + ds) (1.05 + ds) (1.11 + ds) (1.17 + ds) (1.24 + ds) (1.31 + ds) (1.375 + ds)] - 1;
for i = 1:length(rho) % To avoid NaN %
    p_eq_lg(i) = 610.8.*exp(-5.1421.*log(T(i)./273.15)-6828.77*(1./T(i)-1/273.15));
    rho_eq(i) = p_eq_lg(i)./(R*T(i));
    if rho(i)<= rho_eq(i)
        g_nuc(i) = 0;
        S(i) = rho(i)/rho_eq(i);
        continue 
    end
    % If T > 259 or T < 217 : Boundaries % 
    S(i) = (rho(i)/rho_eq(i)); 
%     if T(i) < T_range(end)   
%         T(i) = 217.1; % The valid range is from 259.0 - 217.1 K
%         p_eq_lg(i) = 610.8.*exp(-5.1421.*log(T(i)./273.15)-6828.77*(1./T(i)-1/273.15));
%         rho_eq(i) = p_eq_lg(i)./(R*T(i));
%         S(i) = (rho(i)/rho_eq(i));
%         if S(i) > s_max(end)
%             S(i) = s_max(end);
%         end
%     elseif T(i) > T_range(1)   
%         T(i) = 259; % The valid range is from 259.0 - 217.1 K
%         p_eq_lg(i) = 610.8.*exp(-5.1421.*log(T(i)./273.15)-6828.77*(1./T(i)-1/273.15));
%         rho_eq(i) = p_eq_lg(i)./(R*T(i));
%         S(i) = (rho(i)/rho_eq(i));
%         if S(i) > s_max(1)
%             S(i) = s_max(1);
%         end
%     end
    
    B(i) = 10.^(-32.025880 - 0.070074889 * T(i) + 0.00040240774 * T(i).^2);
    n(i) = 1.1743362 + 0.011626596 .* T(i) + 0.00047597152 .* T(i).^2;
    g_nuc(i) = (B(i).*(S(i)-1).^n(i));
    %g_nuc(i) = (B(i).*(S(i)).^n(i));
end
end  

