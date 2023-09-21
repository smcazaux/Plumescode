% Particle Size Distribution
function [r_e P r_peak] = getParticleSizeDistribution(R, x_d, gn, dx_d, V_d, A_d)
r_e = R(end) - R; % Gives particle size at exit, nucleated from location z
% figure(117)
% semilogy(x_d,r_e);
% ylabel('Particle Radius exit, generated from nucleation at z [m]')
% xlabel('z')
n_z(1) = 0;
if length(dx_d)>1
    for i = 2:length(gn)
        n_z(i) = gn(i-1).*dx_d(i-1)./V_d(i-1);
    end
else
    n_z = gn.*dx_d./V_d;
end
% figure(119)
% semilogy(x_d, n_z)
% ylabel('Particle concentration at the point of their origin [n/m^3]')
% xlabel('z')
drdz = gradient(r_e(:)) ./ gradient(x_d(:));
for i = 1:length(r_e)-1
    dr(i) = r_e(i+1)-r_e(i);
end
dr = [dr 0];
P = gn ./ V_d(end) .* abs(transpose(drdz.^(-1))) .* A_d./A_d(end);% .* abs(dr);%.*dx_d;
%P = gn ./ V_d(end) .* abs(transpose(1./drdz)) .* A_d./A_d(end) .*dx_d;
% figure(111)
% loglog(r_e*10^6, P)
% %ylim([10^(-2) 1.3*max(P)])
% %ylim([10^19 3*10^20])
% %ylim([10^9 2*10^10])
% %xlim([0.001 1000])
% title('Particle Size Distribution')
% ylabel('Particle concentration [n/m^3]')
% xlabel('R(\mu m)_{exit}')
x = find(max(P) == P); % find location of max particle density
r_peak = r_e(x); % find particle radius at max. density
end