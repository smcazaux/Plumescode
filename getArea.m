% getArea
function [A, r] = getArea(n, L, x, dx, Length, channelNr)
if channelNr == 47
        %% Channel 47 %%
    p1 = 4082.5323;
    p2 = -22474.3605;
    p3 = 52477.3896;
    p4 = -67628.5973;
    p5 = 52410.7194;
    p6 = -24930.4765;
    p7 = 7102.3397;
    p8 = -1118.4622;
    p9 = 80.3391;
    p10 = -1.4341;
    p11 = 0.38115;
    r = p1*x.^10 + p2*x.^9 + p3*x.^8 + p4*x.^7 + p5*x.^6 + p6*x.^5 + p7*x.^4 + p8*x.^3 + p9*x.^2 + p10*x + p11;
    A = pi * r .^2;
elseif channelNr == 1
        % Channel 1 %
    p1 = 2425.9223;
    p2 = -12291.7689;
    p3 = 26015.6998;
    p4 = -29617.7449;
    p5 = 19409.4629;
    p6 = -7217.5659;
    p7 = 1363.4469;
    p8 = -80.2421;
    p9 = -7.6893;
    p10 = 0.42572;
    p11 = 0.49841;
    r = p1*x.^10 + p2*x.^9 + p3*x.^8 + p4*x.^7 + p5*x.^6 + p6*x.^5 + p7*x.^4 + p8*x.^3 + p9*x.^2 + p10*x + p11 ;
    A = pi * r.^2;
elseif channelNr == 2
    % quadratic %
    p1 = 1.4855;
    p2 = -2.6113;
    p3 = 1.2562;
    r = p1*x.^2 + p2*x + p3;
% 7th degree polynomial %
%     p1 = -452.5698;
%     p2 = 1779.9012;
%     p3 = -2763.1799;
%     p4 = 2122.8659;
%     p5 = -819.8566;
%     p6 = 140.3533;
%     p7 = -8.4192;
%     p8 = 1.0766;
%     r = p1*x.^7 + p2*x.^6 + p3*x.^5 + p4*x.^4 + p5*x.^3 + p6*x.^2 + p7*x + p8; 
    A = pi*r.^2;    
elseif channelNr == 11
    A = 1 + 2.2 * (x - 1.5).^2; % Sub-sonic to super-sonic [Anderson]
    r = sqrt(A./pi);
elseif channelNr == 12 % SHOCK CAPTURING
    for i = 1:length(x)
        if x(i) <= 3
            A(i) = 1 + 2.2 * (x(i) - 1.5).^2; % Sub-sonic to super-sonic [Anderson]
        else
            A(i) = 1 + 0.2223*(x(i) - 1.5).^2;
        end
    end
    r = sqrt(A./pi);    
elseif channelNr == 12
    for i = 1:length(x)
        if x(i) <= 1.5
            A(i) = 1.0 + 2.2*(x(i)-1.5).^2; %Nozzle geometry
        end
        if x(i) >= 1.5
            A(i) = 1.0 + 0.2223*(x(i) - 1.5).^2;
        end
    end

    r = sqrt(A./pi);
elseif channelNr == 101
    p1 = 3.4236e-18;
    p2 = -1.5302e-15;
    p3 = 2.9007e-13;
    p4 = -3.0336e-11;
    p5 = 1.9041e-09;
    p6 = -7.2944e-08;
    p7 = 1.6557e-06;
    p8 = -2.0574e-05;
    p9 = 0.00012081;
    p10 = -0.00022353;
    p11 = 0.026483;
    r = p1*x.^10 + p2*x.^9 + p3*x.^8 + p4*x.^7 + p5*x.^6 + p6*x.^5 + p7*x.^4 + p8*x.^3 + p9*x.^2 + p10*x + p11; 
    %[x, A, rho, T, V, c, R, f, gn, S, data, r] = getData(channelNr)
    A = pi*r.^2;
end
%dx = 0.05;
%  % NOZZLE 2 %
%x = (0:dx:1);
% n = length(x);
% A = pi*(-0.6072.*x.^4 + 2.753.*x.^2 - 3.046.*x + 1).^2;

% % normale Nozzle %
% dx = 0.05;
% x = (0:dx:1);
% A = 1.0 + 0.1*(x-0.5).^2;
% r = sqrt(A/pi);% raduis 
%% Initial conditions for the Anderson Book 
% for i = 1 : length(x)
%     if x(i) >= 0 && x(i) <= 0.5
%         rho(i) = 1;
%         T(i) = 1;
%     elseif x(i) >= 0.5 && x(i) <= 1.5
%         rho(i) = 1 - 0.366 * (x(i) - 0.5);
%         T(i) = 1 - 0.167 * (x(i) - 0.5);
% %      elseif x(i) >= 1.5 && x(i) <= 2.1
% %         rho(i) = 0.634 - 0.702 * (x(i) - 1.5);
% %         T(i) = 0.833 - 0.4908 * (x(i) - 1.5);
% %      elseif x(i) >= 2.1 && x(i) <= 3.0
% % %         rho(i) = 0.5892 - 0.10228 * (x(i) - 2.1);
% %         T(i) = 0.93968 - 0.0622 * (x(i) - 2.1);
%     elseif x(i) >= 1.5 && x(i) <= 3.5 
%         rho(i) = 0.634 - 0.3879 * (x(i) - 1.5);
%         %T(i) = 0.833 - 0.3507 * (x(i) - 1.5);
%         T(i) = 0.833 - 0.1 * (x(i) - 1.5);
%      end
% end 