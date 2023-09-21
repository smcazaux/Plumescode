% Energy release
QQ = rho_d.*A_d.*V_d.*f.*L_h./dx_d;
hh = rho_d.*A_d.*V_d.*cv.*T_d./dx_d;
vv = rho_d.*A_d.*V_d.*0.5.*V_d.^2./dx_d;
%vs isentropic d(CvT), d(1/2V^2) 
figure(438)
plot(x_d, QQ,'linewidth',2)
hold on
plot(x_d, hh,'linewidth',2)
plot(x_d, vv,'linewidth',2)
legend('Q','CvT','1/2V^2')
grid on
ylabel('[J/ms]')
xlabel('z [m]')