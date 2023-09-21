% Paper figures !!!
% Run twice, once with baseline data, once with new simulation data to
% generate plots with 2 lines.
clear all
baseline = true;
if baseline == true % Put in file-location of baseline 
    % Baseline
    oldResults_xl = 'C:\Users\smcazaux\Desktop\Channels\Isentropic\J_DB2.0_V0fix_Channel51_rho0.004626_T268.9086_fp(i)_pe0_Cx0.3_i10000_n201_C0.1_b0.2_De9';
else                % Put in file-location of simulation 
    % Wall interactions
    oldResults_xl = 'C:\Users\nicky\Documents\Thesis\CompleteModel\Figures\J_DB2.0_Walls_V0fix_ChannelNr51_rho0.0042743_T261.914_pe0_Cx0.3_i11000_n201_C0.1_b0.2_L150_De9_m_stick3_r_sub3_crackL200';
end
% Load data from txt. files
results = load([oldResults_xl, '\runData.txt']);
constants = load([oldResults_xl, '\constantsData.txt']);
fullResults = load([oldResults_xl, '\runFullData.txt']);
residuals = load([oldResults_xl, '\Residuals.txt']);
% Dimensional Results %
T_d = fullResults(:,1);
rho_d = fullResults(:,2);
p_d = fullResults(:,3);
V_d = fullResults(:,4);
M = fullResults(:,5);
f = fullResults(:,6);
R = fullResults(:,7);
S = fullResults(:,8);
%S = fullResults(:,6);
gn = fullResults(:,9);
fp = fullResults(:,10);
dRdx = fullResults(:,11);
E_dotL = fullResults(:,12);
E = fullResults(:,13);
A_d = fullResults(:,14);
%E = fullResults(:,7);
%A_d = fullResults(:,8);
x_d = fullResults(:,15);
r = fullResults(:,16);
T_w = fullResults(:,17);
% x_d = fullResults(:,9);
% r = fullResults(:,10);
% T_w = fullResults(:,11);
r_e = fullResults(:,18);
P = fullResults(:,19);
m_d = fullResults(:,20);
%m_d = fullResults(:,12);
dx_d = x_d(end)/200;

dV(1) = 0; % Velocity increase per cell
for i = 1:length(V_d)-1
    dV(i+1) = V_d(i+1) - V_d(i);
end
% Transpose %
T_d = transpose(T_d); rho_d = transpose(rho_d); p_d = transpose(p_d); V_d = transpose(V_d); M = transpose(M);
f = transpose(f); R = transpose(R); 
S = transpose(S); 
gn = transpose(gn); fp = transpose(fp);dRdx = transpose(dRdx); E_dotL = transpose(E_dotL); 
A_d = transpose(A_d); x_d = transpose(x_d); 
r_e = transpose(r_e); P = transpose(P);
T_w = transpose(T_w);  r = transpose(r); m_d = transpose(m_d);
%x_d = A_d
% constants %
n = length(x_d);
L = constants(1);
L_h = constants(2);
cv = constants(3);
T_ref = constants(4);
Rg = constants(5);
rho_ref = constants(6);
y = constants(7);
As = constants(8);
crackLength = constants(9); 
m_stick = constants(10); 
r_sub = constants(11); 

%Q = mean(m_d);
% Results
% R_peak = results(15);
%hours = results(16)/3600
% Residuals
resU1 = residuals(:,1);
resU2 = residuals(:,2);
resU3 = residuals(:,3);
x_d_throat = constants(12);
throat = find(x_d == x_d_throat);

%x_d_throat = 103.5; % change
p_ref = p_d(1);
V_ref = V_d(1);
M_ref = M(1);
circShape = false; %% change if necesarry (C100 & C101)
if circShape == true
    r_d = sqrt(A_d./pi);
    c_d = 2.*pi.*r_d; % Circumference of the channel [m]    
else
    r_d = A_d./(2.*crackLength);
    c_d = 4.*r_d + 2.*crackLength; % Circumference of the channel [m]
end
n = length(x_d);
dx_d = L/(n-1);
dx = dx_d./L;
%x_d = linspace(0,150,201);
p_eq = 610.8.*exp(-5.1421.*log(T_d./273.15)-6828.77*(1./T_d-1/273.15));
rho_eq = p_eq./(Rg.*T_d); 
drho = rho_d - rho_eq;
%[r_e P r_peak] = getParticleSizeDistribution(R, x_d, gn, dx_d, V_d, A_d);
% figure(127)
% semilogx(r_e*10^6, P./max(P), '-','linewidth',1.2)
% hold on
% ylabel('Non-dimensional number density [-]')
% xlabel('R [\mu m]')
% xlim([0.1 75])
% set(gca,'FontSize',9)
% grid on
% %legend('Multi-phase','Throat','D_{exit}/D_{throat} = 2.4','Length','Irregular','Wall interactions')
% %
% dV
figure(259)
if baseline == true
    i = 1;
    j = 1.25;
    k = 9;
    subplot(i,1,1)
    plot(x_d, dV, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dV [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
else
        i = 1;
        j = 1.5;
        k = 9;
    subplot(i,1,1)
    plot(x_d, dV, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dV [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
end
% %% MTrhoS %%
figure(123)
if baseline == true
    i = 4;
    j = 1.25;
    k = 9;
    subplot(i,1,1)
    plot(x_d, M, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    set(gca,'ytick',linspace(0,2,11))
    ylim([min(M) 1.401])
    ylabel('M [-]')
    xlabel('z [m]')
    xlim([0 x_d(end)])
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'g--','linewidth',j)
    set(gca,'ytick',linspace(200,300,6))
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,3)
    plot(x_d, rho_d, 'm--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,4)
    plot(x_d, S, 'b--','linewidth',j)
    hold on
    ylim([0 30])
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
else
        i = 4;
        j = 1.5;
        k = 9;
    subplot(i,1,1)
    plot(x_d, M, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    set(gca,'ytick',linspace(0.4,1.6,4))
    ylim([min(M) max(M)])
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(a)', 'FontSize', k);
    subplot(i,1,2)
    plot(x_d, T_d, 'g','linewidth',j)
    hold on
    set(gca,'ytick',linspace(200,300,6))%set(gca,'ytick',linspace(100,300,6))
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(b)', 'FontSize', k);
    subplot(i,1,3)
    plot(x_d, rho_d, 'm','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(c)', 'FontSize', k);
    subplot(i,1,4)
    plot(x_d, S, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(d)', 'FontSize', k);
end
% drho, S, gn, dRdx
i=4; % number of rows
%
figure(333)
if baseline == true
    j = 1.25; % line width
    subplot(i,1,4)
    hold on
    plot(x_d, drho,'--','Color',[0.085, 0.325, 0.098],'linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,2)
    plot(x_d, f, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,1)
    semilogy(x_d, gn,'--', 'Color',[0.7, 0.680, 0.184],'linewidth',j)
    hold on
    ylim([10^0 10^12]) 
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    xlim([0 150])
    subplot(i,1,3)
    semilogy(x_d, dRdx, '--','color',[0.9100    0.4100    0.1700],'linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
else
    j=1.5;
    subplot(i,1,4)
    plot(x_d, drho,'Color',[0.085, 0.325, 0.098],'linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(d)', 'FontSize', k);
    set(gca,'ytick',linspace(-1e-3,9e-3,11))
    subplot(i,1,2)
    plot(x_d, f, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    set(gca,'ytick',linspace(0,0.2,11))
    title('(b)', 'FontSize', k);
    subplot(i,1,1)
    %plot(x_d, gn, 'Color',[0.7, 0.680, 0.184],'linewidth',j)
    semilogy(x_d, gn, 'Color',[0.7, 0.680, 0.184],'linewidth',j)
    %ylim([max(gn)/10^2 max(gn)])
    %ylim([10^0 10^12])
    ylim([10^0 10^12])
    set(gca,'ytick',logspace(0,12,4))
    hold on
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(a)', 'FontSize', k);
    %set(gca,'ytick',linspace(round(max(gn))/10^4,round(max(gn)),5))
    xlim([0 max(x_d)])
    subplot(i,1,3)
    semilogy(x_d, dRdx, 'color',[0.9100    0.4100    0.1700],'linewidth',j)
    hold on
    set(gca,'ytick',logspace(-10,-5,6))
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    title('(c)', 'FontSize', k);
end
figure(172)
plot(x_d,m_d,'linewidth',2)
ylabel('Mass flow [kg/s]')
xlabel('z [m]')
grid on
hold on
set(gca,'FontSize',k)
xline(x_d_throat,'--k')


% Nozzle
figure(16);
plot([x_d_throat x_d_throat],[r_d(throat) -r_d(throat)],'k--')
hold on
plot(x_d,r_d,'k','linewidth',2)
plot(x_d,-r_d,'k','linewidth',2)
hold off
D_exit = r_d(end)*2;
D_throat = 2*r_d(throat);
legend(['D_{throat} = ', num2str(round(D_throat,1)),' m'],['D_{exit} = ', num2str(D_exit), ' m']);
ylabel('Channel Radius [m]')
xlabel('z [m]')
% Residuals
figure(321)
hold off
semilogy(resU1,'linewidth',1.5)
hold on
semilogy(resU2,'linewidth',1.5)
semilogy(resU3,'linewidth',1.5)
legend('U1','U2','U3')
% Plots %%
%% dependencies (T, S, gn)
figure(146)
gnnn = [];
iii = [];
T_ddd = [];
SSS = [];
x_ddd = [];
for i = 1:length(gn)
    if gn(i) > 10^9.6
        iii = [iii,i];
        gnnn = [gnnn, gn(i)];
        T_ddd = [T_ddd, T_d(i)];
        SSS = [SSS, S(i)];
        x_ddd = [x_ddd, x_d(i)];
    end
end
T_dd = interp1(x_ddd,T_ddd,x_ddd(1):dx_d/2000:x_ddd(end));
x_dddd = x_ddd(1):(x_ddd(end)-x_ddd(1))./(length(T_dd)-1):x_ddd(end); 
p = polyfit(x_dddd, T_dd, 2);
T_dd = p(1).*x_dddd.^2+p(2).*x_dddd + p(3);
gnn = interp1(x_ddd,gnnn,x_ddd(1):dx_d/2000:x_ddd(end));
p = polyfit(x_dddd, gnn, 2);
gnn = p(1).*x_dddd.^2+p(2).*x_dddd.^1 + p(3);
% gnn = p(1).*x_dddd.^3+p(2).*x_dddd.^2 + p(3).*x_dddd.^1 + p(4);
% gnn = p(1).*x_dddd.^4+p(2).*x_dddd.^3 + p(3).*x_dddd.^2 + p(4).*x_dddd.^1 + p(5);
% gnn = p(1).*x_dddd.^5+p(2).*x_dddd.^4 + p(3).*x_dddd.^3 + p(4).*x_dddd.^2 + p(5).*x_dddd + p(6);
SS = interp1(x_ddd,SSS,x_ddd(1):dx_d/2000:x_ddd(end));
p = polyfit(x_dddd, SS, 2);
SS = p(1).*x_dddd.^2+p(2).*x_dddd + p(3);
scatter(T_dd,gnn,500,SS,'.')
hold on
set(gca,'yscale','log')
%scatter(fitnessValue, verticalTime, 500, objectiveVariables(:,5).*180/pi, '.');
ylabel('\gamma_{nuc} [1/m^3s]')
xlabel('T [K]')
grid on
% hc = colorbar;
% cb = linspace(min(p_d),max(p_d),20);
% set(hc, 'ylim', [min(p_d) max(p_d)]);
h = colorbar;
set(get(h,'label'),'string','S [-]');
%set(gca,'ColorScale','log')
set(gca, 'fontsize',11)
%ylim([min(gnnn) max(gnnn)*1.1])
%set(hc, 'YTick',cb, 'YTickLabel',cb)

%% phase diagram plot
T_range1 = linspace(240,273.16,300);
T_range2 = linspace(273.16,280,100); 
T_range3 = linspace(273.16,273.15,300);
p_range3 = linspace(611.657,101325,300);
%T_range = linspace(
p_eq_lg1 = p_eq_lg(T_range2);
p_eq_lg2 = p_eq_lg(T_range1);
figure(147)
semilogy(T_range2,p_eq_lg1,'b-','linewidth',2)
hold on
semilogy(T_range1,p_eq_lg2,'b-','linewidth',2)
semilogy(T_range3, p_range3, 'b-', 'linewidth',2)
scatter(T_d3,p_d3,500,S3,'.')
%scatter(fitnessValue, verticalTime, 500, objectiveVariables(:,5).*180/pi, '.');
ylim([  10^2 1000])
xlim([175 280])
ylabel('p [Pa]')
xlabel('T [K]')
grid on
% hc = colorbar;
% cb = linspace(min(p_d),max(p_d),20);
% set(hc, 'ylim', [min(p_d) max(p_d)]);
h = colorbar;
set(get(h,'label'),'string','S [-]');
set(gca, 'fontsize',11)
%set(hc, 'YTick',cb, 'YTickLabel',cb)
% T S
%% %% MTrhoS %%
figure(123)
if baseline == true
    i = 2;
    j = 1.25;
    k = 9;
    subplot(i,1,1)
    plot(x_d, T_d, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,2)
    plot(x_d, S, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
else
        i = 2;
        j = 1.5;
        k = 9;
    subplot(i,1,1)
    plot(x_d, T_d, 'g','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,2)
    plot(x_d, S, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
end


%% %% MTrhoS %%
figure(123)
if baseline == true
    i = 4;
    j = 1.25;
    
else
        i = 4;
        j = 1.5;
    subplot(i,1,1)
    plot(x_d, M, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,2)
    plot(x_d, T_d, 'g','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,3)
    plot(x_d, rho_d, 'm','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
    subplot(i,1,4)
    plot(x_d, S, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    set(gca,'FontSize',k)
    grid on
end
%%
clear all
% Schmidt c100
baseline = true;
Rg = 461.4472;
x_d_throat = 139.2186;
throat = 349;
[x, A, rho_d, T_d, V_d, c, R, f, gn, S, data, r, L, n, dx, Length, M, p_d, x_d, dx_d] = getData_Verification(100);
A_d = [57.2500747137866,56.5878877354116,55.9017566377783,55.2256510903809,54.5911868323818,54.0261579965598,53.5533666606631,53.1897676986551,52.9459247091292,52.8257608336307,52.8265861186831,52.9393878395277,53.1493774454339,53.4367926204614,53.7779512615350,54.1465436354819,54.5151296884929,54.8567830524000,55.1467963549685,55.3643397038636,55.4939512364824,55.5267395871598,55.4611948264173,55.3035357911766,55.0675638441366,54.7740398981456,54.4496458867014,54.1256269291149,53.8362310351426,53.6170667809381,53.5034865879271,53.5290777071127,53.7243107131859,54.1153633017799,54.7231121757737,55.5622728015608,56.6406680547289,57.9586212524506,59.5084927079737,61.2744053889127,63.2322271709964,65.3498877235175,67.5881023665633,69.9015514628266,72.2405237500532,74.5529805541885,76.7869426596502,78.8930515782761,80.8271205184071,82.5524740906017,84.0418832026847,85.2789324946751,86.2587081756940,86.9877574028996,87.4833377602670,87.7720380247312,87.8879014952206,87.8702153296646,87.7611413063187,87.6033561437139,87.4378466920791,87.3019725108325,87.2278718027520,87.2412520880487,87.3605785644075,87.5966529112863,87.9525632597253,88.4239802300116,88.9997712933755,89.6629029826173,90.3915950639206,91.1606814118478,91.9431193548277,92.7115745703190,93.4399951806496,94.1050798835918,94.6875436867884,95.1730929368913,95.5530391201703,95.8245069669624,95.9902238898785,96.0579109936946,96.0393268764866,95.9490407653973,95.8030288411099,95.6171959377914,95.4059245891029,95.1807461594996,94.9492166928540,94.7140652940895,94.4726670136586,94.2168761964714,93.9332400806331,93.6035954232102,93.2060321986183,92.7161873727002,92.1088086086327,91.3595038189068,90.4465701807394,89.3527788599409,88.0669827881967,86.5854174826758,84.9125809206164,83.0616078811223,81.0540948013959,78.9193788918619,76.6933243268041,74.4167124326332,72.1333659912782,69.8881555637102,67.7250359029266,65.6852435866212,63.8057561771275,62.1180738975044,60.6473436100909,59.4118084525636,58.4225403248147,57.6833998817876,57.1911703983899,56.9358256654450,56.9009134601404,57.0640592133715,57.3976131509098,57.8694732998896,58.4441132410228,59.0838269227884,59.7501755761181,60.4055884520836,61.0150358761743,61.5476663800296,61.9782849075473,62.2885497951129,62.4677832602989,62.5133215749938,62.4303725997540,62.2313939562892,61.9350483323037,61.5648273971895,61.1474583350100,60.7112151272200,60.2842509671719,59.8930512756562,59.5610829310346,59.3076893239637,59.1472570861480,59.0886619826402,59.1349899060824,59.2835237163790,59.5259858158155,59.8490268878254,60.2349500809806,60.6626546922369,61.1087730080227,61.5489589357957,61.9592695003937,62.3175634227606,62.6048285211618,62.8063448910521,62.9125959319455,62.9198548611467,62.8303991144056,62.6523360588916,62.3990566763338,62.0883648685753,61.7413547821970,61.3811241879286,61.0314172642898,60.7152856781588,60.4538447358850,60.2651847306262,60.1634798785358,60.1583214272646,60.2542895961856,60.4507715319879,60.7420285861198,61.1175139980849,61.5624390394205,62.0585795602659,62.5853042566999,63.1207907397957,63.6433770274221,64.1329770970577,64.5724731330742,64.9489876472547,65.2549385940916,65.4887913759581,65.6554428316484,65.7662016229534,65.8383630202065,65.8944092081501,65.9608941774924,66.0670912695182,66.2434894419087,66.5202214027997,66.9254951332092,67.4840838865268,68.2159132785800,69.1347721421062,70.2471697655696,71.5513672563878,73.0366238559794,74.6827164671515,76.4598070041152,78.3287412463849,80.2418589984722,82.1443746749308,83.9763489518516,85.6752182817136,87.1787856575794,88.4285114982488,89.3728877733830,89.9706411443423,90.1934996926295,90.0282770803461,89.4780777047632,88.5625018957301,87.3168226392351,85.7902029492744,84.0431129449596,82.1441758577361,80.1667133039163,78.1852672057158,76.2723488483561,74.4956099721423,72.9155560650602,71.5838404434862,70.5421023320956,69.8212546377377,69.4410958221925,69.4101189309148,69.7254179102757,70.3726403801210,71.3259967255574,72.5483953570332,73.9918209088402,75.5980957749440,77.3001593678046,79.0239624793061,80.6910100732903,82.2215034368984,83.5379439609795,84.5689795274326,85.2532135577359,85.5426665665577,85.4055866207933,84.8283492251845,83.8162642273241,82.3932080020889,80.6001107213520,78.4924367286416,76.1368872805537,73.6076180627654,70.9822918642163,68.3382774800123,65.7492623073480,63.2824757171169,60.9966339368501,58.9406273038385,57.1528895603801,55.6613266305616,54.4836459647022,53.6279198091232,53.0932350288501,52.8703228376316,52.9421155296688,53.2842341130201,53.8654607258255,54.6482845303276,55.5896238310476,56.6418183519854,57.7539555281045,58.8735482531760,59.9485262144449,60.9294474867198,61.7717901801673,62.4381530344692,62.9001839230320,63.1400681754100,63.1514430243899,62.9396558593906,62.5213454867877,61.9233890975600,61.1813148271835,60.3373234228715,59.4380874397963,58.5325001644544,57.6695297463870,56.8963003224310,56.2564770503444,55.7889831328001,55.5270316093681,55.4974196034647,55.7200127090359,56.2073446538877,56.9642718136758,57.9876503677892,59.2660404544708,60.7794797822648,62.4994016404037,64.3887927467130,66.4026902707863,68.4891025711411,70.5904054615390,72.6452187658242,74.5907124661554,76.3652353412614,77.9111095224588,79.1773990665719,80.1224449155596,80.7159653958792,80.9405506357275,80.7924278833721,80.2814370659382,79.4302246318334,78.2727306345426,76.8521014038449,75.2182017317733,73.4249222725735,71.5274785640934,69.5798792686612,67.6327068897983,65.7313100618871,63.9144590407621,62.2134715339852,60.6517796505219,59.2448838415492,58.0006273839041,56.9197242011598,55.9964808216981,55.2196661165135,54.5734959128187,54.0387099225761,53.5937231211482,53.2158318249540,52.8824470077977,52.5723160804933,52.2666825980271,51.9503245696193,51.6124091775710,51.2471066160519,50.8539188013150,50.4376986404574,50.0083597489583,49.5803014494745,49.1715958059993,48.8029990783086,48.4968571740784,48.2759728258858,48.1624923905546,48.1768548886007,48.3368286767398,48.6566457656238,49.1462335713779,49.8105409778809,50.6489605125526,51.6548599354092,52.8152517369998,54.1106439651837,55.5151261285499,56.9967458206245,58.5182226191503,60.0380250505525,61.5118054170755,62.8941495322635,64.1405589126908,65.2095475148366,66.0647093065436,66.6766013419535,67.0242922085942,67.0964480191599,66.8918653561680,66.4194084174817,65.6973602203382,64.7522486454033,63.6172512351933];
r_d = sqrt(A_d./pi);    
p_eq = 610.8.*exp(-5.1421.*log(T_d./273.15)-6828.77*(1./T_d-1/273.15));
rho_eq = p_eq./(Rg.*T_d); 
drho = rho_d - rho_eq;
%% Schmidt c101
clear all

baseline = false;
Rg = 461.4472;
x_d_throat = 103.6138;
throat = 260;
[x, A, rho_d, T_d, V_d, c, R, f, gn, S, data, r, L, n, dx, Length, M, p_d, x_d, dx_d] = getData_Verification(101);
A_d = [57.2500747137866,56.5878877354116,55.9017566377783,55.2256510903809,54.5911868323818,54.0261579965598,53.5533666606631,53.1897676986551,52.9459247091292,52.8257608336307,52.8265861186831,52.9393878395277,53.1493774454339,53.4367926204614,53.7779512615350,54.1465436354819,54.5151296884929,54.8567830524000,55.1467963549685,55.3643397038636,55.4939512364824,55.5267395871598,55.4611948264173,55.3035357911766,55.0675638441366,54.7740398981456,54.4496458867014,54.1256269291149,53.8362310351426,53.6170667809381,53.5034865879271,53.5290777071127,53.7243107131859,54.1153633017799,54.7231121757737,55.5622728015608,56.6406680547289,57.9586212524506,59.5084927079737,61.2744053889127,63.2322271709964,65.3498877235175,67.5881023665633,69.9015514628266,72.2405237500532,74.5529805541885,76.7869426596502,78.8930515782761,80.8271205184071,82.5524740906017,84.0418832026847,85.2789324946751,86.2587081756940,86.9877574028996,87.4833377602670,87.7720380247312,87.8879014952206,87.8702153296646,87.7611413063187,87.6033561437139,87.4378466920791,87.3019725108325,87.2278718027520,87.2412520880487,87.3605785644075,87.5966529112863,87.9525632597253,88.4239802300116,88.9997712933755,89.6629029826173,90.3915950639206,91.1606814118478,91.9431193548277,92.7115745703190,93.4399951806496,94.1050798835918,94.6875436867884,95.1730929368913,95.5530391201703,95.8245069669624,95.9902238898785,96.0579109936946,96.0393268764866,95.9490407653973,95.8030288411099,95.6171959377914,95.4059245891029,95.1807461594996,94.9492166928540,94.7140652940895,94.4726670136586,94.2168761964714,93.9332400806331,93.6035954232102,93.2060321986183,92.7161873727002,92.1088086086327,91.3595038189068,90.4465701807394,89.3527788599409,88.0669827881967,86.5854174826758,84.9125809206164,83.0616078811223,81.0540948013959,78.9193788918619,76.6933243268041,74.4167124326332,72.1333659912782,69.8881555637102,67.7250359029266,65.6852435866212,63.8057561771275,62.1180738975044,60.6473436100909,59.4118084525636,58.4225403248147,57.6833998817876,57.1911703983899,56.9358256654450,56.9009134601404,57.0640592133715,57.3976131509098,57.8694732998896,58.4441132410228,59.0838269227884,59.7501755761181,60.4055884520836,61.0150358761743,61.5476663800296,61.9782849075473,62.2885497951129,62.4677832602989,62.5133215749938,62.4303725997540,62.2313939562892,61.9350483323037,61.5648273971895,61.1474583350100,60.7112151272200,60.2842509671719,59.8930512756562,59.5610829310346,59.3076893239637,59.1472570861480,59.0886619826402,59.1349899060824,59.2835237163790,59.5259858158155,59.8490268878254,60.2349500809806,60.6626546922369,61.1087730080227,61.5489589357957,61.9592695003937,62.3175634227606,62.6048285211618,62.8063448910521,62.9125959319455,62.9198548611467,62.8303991144056,62.6523360588916,62.3990566763338,62.0883648685753,61.7413547821970,61.3811241879286,61.0314172642898,60.7152856781588,60.4538447358850,60.2651847306262,60.1634798785358,60.1583214272646,60.2542895961856,60.4507715319879,60.7420285861198,61.1175139980849,61.5624390394205,62.0585795602659,62.5853042566999,63.1207907397957,63.6433770274221,64.1329770970577,64.5724731330742,64.9489876472547,65.2549385940916,65.4887913759581,65.6554428316484,65.7662016229534,65.8383630202065,65.8944092081501,65.9608941774924,66.0670912695182,66.2434894419087,66.5202214027997,66.9254951332092,67.4840838865268,68.2159132785800,69.1347721421062,70.2471697655696,71.5513672563878,73.0366238559794,74.6827164671515,76.4598070041152,78.3287412463849,80.2418589984722,82.1443746749308,83.9763489518516,85.6752182817136,87.1787856575794,88.4285114982488,89.3728877733830,89.9706411443423,90.1934996926295,90.0282770803461,89.4780777047632,88.5625018957301,87.3168226392351,85.7902029492744,84.0431129449596,82.1441758577361,80.1667133039163,78.1852672057158,76.2723488483561,74.4956099721423,72.9155560650602,71.5838404434862,70.5421023320956,69.8212546377377,69.4410958221925,69.4101189309148,69.7254179102757,70.3726403801210,71.3259967255574,72.5483953570332,73.9918209088402,75.5980957749440,77.3001593678046,79.0239624793061,80.6910100732903,82.2215034368984,83.5379439609795,84.5689795274326,85.2532135577359,85.5426665665577,85.4055866207933,84.8283492251845,83.8162642273241,82.3932080020889,80.6001107213520,78.4924367286416,76.1368872805537,73.6076180627654,70.9822918642163,68.3382774800123,65.7492623073480,63.2824757171169,60.9966339368501,58.9406273038385,57.1528895603801,55.6613266305616,54.4836459647022,53.6279198091232,53.0932350288501,52.8703228376316,52.9421155296688,53.2842341130201,53.8654607258255,54.6482845303276,55.5896238310476,56.6418183519854,57.7539555281045,58.8735482531760,59.9485262144449,60.9294474867198,61.7717901801673,62.4381530344692,62.9001839230320,63.1400681754100,63.1514430243899,62.9396558593906,62.5213454867877,61.9233890975600,61.1813148271835,60.3373234228715,59.4380874397963,58.5325001644544,57.6695297463870,56.8963003224310,56.2564770503444,55.7889831328001,55.5270316093681,55.4974196034647,55.7200127090359,56.2073446538877,56.9642718136758,57.9876503677892,59.2660404544708,60.7794797822648,62.4994016404037,64.3887927467130,66.4026902707863,68.4891025711411,70.5904054615390,72.6452187658242,74.5907124661554,76.3652353412614,77.9111095224588,79.1773990665719,80.1224449155596,80.7159653958792,80.9405506357275,80.7924278833721,80.2814370659382,79.4302246318334,78.2727306345426,76.8521014038449,75.2182017317733,73.4249222725735,71.5274785640934,69.5798792686612,67.6327068897983,65.7313100618871,63.9144590407621,62.2134715339852,60.6517796505219,59.2448838415492,58.0006273839041,56.9197242011598,55.9964808216981,55.2196661165135,54.5734959128187,54.0387099225761,53.5937231211482,53.2158318249540,52.8824470077977,52.5723160804933,52.2666825980271,51.9503245696193,51.6124091775710,51.2471066160519,50.8539188013150,50.4376986404574,50.0083597489583,49.5803014494745,49.1715958059993,48.8029990783086,48.4968571740784,48.2759728258858,48.1624923905546,48.1768548886007,48.3368286767398,48.6566457656238,49.1462335713779,49.8105409778809,50.6489605125526,51.6548599354092,52.8152517369998,54.1106439651837,55.5151261285499,56.9967458206245,58.5182226191503,60.0380250505525,61.5118054170755,62.8941495322635,64.1405589126908,65.2095475148366,66.0647093065436,66.6766013419535,67.0242922085942,67.0964480191599,66.8918653561680,66.4194084174817,65.6973602203382,64.7522486454033,63.6172512351933];
r_d = sqrt(A_d./pi);    
p_eq = 610.8.*exp(-5.1421.*log(T_d./273.15)-6828.77*(1./T_d-1/273.15));
rho_eq = p_eq./(Rg.*T_d); 
drho = rho_d - rho_eq;

%% Schmidt V&V
% drho, S, gn, dRdx
i=3; % number of rows
k = 9;
figure(333)
if baseline == true
    j = 0.5; % line width
    subplot(i,1,3)
    plot(x_d, f, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,1)
    semilogy(x_d, gn,'b','linewidth',j)
    hold on
    ylim([1 max(gn)]) 
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    grid on
    xlim([x_d(1) x_d(end)])
    set(gca,'FontSize',k)
    subplot(i,1,2)
    semilogy(x_d, R, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('R_{max} [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
else
    j = 2; % line width
    subplot(i,1,3)
    plot(x_d, f, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,1)
    semilogy(x_d, gn,'--', 'Color',[0.7, 0.680, 0.184],'linewidth',j)
    hold on
    ylim([1 max(gn)]) 
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    xlim([x_d(1) x_d(end)])
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,2)
    semilogy(x_d, R,'--','color',[0.9100    0.4100    0.1700],'linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('R_{max} [m]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    
end
i = 4;
figure(334)
if baseline == true
    j = 0.5; % line width
    subplot(i,1,1)
    plot(x_d, M, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,2)
    plot(x_d, T_d, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,3)
    plot(x_d, rho_d, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,4)
    plot(x_d, S, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    
else
    j = 2; % line width
    subplot(i,1,1)
    plot(x_d, M, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,2)
    plot(x_d, T_d, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,3)
    plot(x_d, rho_d, 'm--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    subplot(i,1,4)
    plot(x_d, S, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    set(gca,'FontSize',k)
    
end
%% Width (double x-axis),  MTrhoS %%
figure(123)
if baseline == true
    i = 4;
    j = 1.25;
    subplot(i,1,1)
    plot(x_d,M,'r--','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('M [-]')
    x2 = x_d*2;
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,M,'r--','linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'g--','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('T [K]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,T_d, 'g--','linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,3)
    plot(x_d, rho_d, 'm--','linewidth',j)
    hold on
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,rho_d, 'm--','linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,4)
    plot(x_d, S, 'b--','linewidth',j)
    ylabel('S [-]')
    xlabel('z [m]')
    hold on
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,S, 'b--','linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
else
        i = 4;
        j = 1.5;
        %x_d = x_d./2;
    subplot(i,1,1)
    plot(x_d./2,M,'r','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('M [-]')
    x2 = x_d;
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,M)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,2)
    plot(x_d./2, T_d, 'g','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('T [K]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,T_d)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,3)
    plot(x_d./2, rho_d, 'm','linewidth',j)
    hold on
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,rho_d)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,4)
    plot(x_d./2, S, 'b','linewidth',j)
    ylabel('S [-]')
    xlabel('z [m]')
    hold on
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,S)
    xline(x_d_throat,'--k')
    grid on
end
%% Length (double x-axis),  ynuc_f_dRdz_drho %%
figure(123)
if baseline == true
    i = 4;
    j = 1.25;
    subplot(i,1,1)
    semilogy(x_d,gn,'--','Color',[0.7, 0.680, 0.184],'linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('\gamma_{nuc} [1/m^3s]')
    x2 = x_d*2;
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    semilogy(hAx(2),x2,gn,'Color',[0.7, 0.680, 0.184],'linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,2)
    plot(x_d, f, 'k--','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('f [-]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,f, 'k--','linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,3)
    plot(x_d, dRdx, '--','color',[0.9100    0.4100    0.1700],'linewidth',j)
    hold on
    ylabel('dR/dz [-]')
    xlabel('z [m]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,dRdx, '--','color',[0.9100    0.4100    0.1700],'linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
    subplot(i,1,4)
    plot(x_d, drho,'--','Color',[0.085, 0.325, 0.098],'linewidth',j)
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    hold on
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,drho,'--','Color',[0.085, 0.325, 0.098],'linewidth',j)
    xline(x_d_throat*2,'--k')
    grid on
else
        i = 4;
        j = 1.5;
        %x_d = x_d./2;
        subplot(i,1,1)
    semilogy(x_d./2,gn,'Color',[0.7, 0.680, 0.184],'linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('\gamma_{nuc} [1/m^3s]')
    x2 = x_d;
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    semilogy(hAx(2),x2,gn,'Color',[0.7, 0.680, 0.184],'linewidth',j)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,2)
    plot(x_d./2, f, 'k','linewidth',j)
    hold on
    xlabel('z [m]')
    ylabel('f [-]')
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,f, 'k','linewidth',j)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,3)
    plot(x_d./2, dRdx, 'color',[0.9100    0.4100    0.1700],'linewidth',j)
    ylabel('dR/dz [-]')
    xlabel('z [m]')
    hold on
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,dRdx, 'color',[0.9100    0.4100    0.1700],'linewidth',j)
    xline(x_d_throat,'--k')
    grid on
    subplot(i,1,4)
    plot(x_d./2, drho,'Color',[0.085, 0.325, 0.098],'linewidth',j)
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    hold on
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','origin','color','none');
    hold(hAx(2),'on')
    plot(hAx(2),x2,drho,'Color',[0.085, 0.325, 0.098],'linewidth',j)
    xline(x_d_throat,'--k')
    grid on
end
%% E vs E (C51 nuc + wall)
% c_stick = 1 (nr. 4), r_sub = 3 (H&H)
[wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, 'xxx', 1, 2000, 1, 4, T_w, c_d, r_sub); % [kg/s]
A = 3.63*10^(12); % [Pa]
B =  6147; % [K]
pw = A*exp(-B./T_w);
EEE = -2.*(p_d./sqrt(2.*pi.*Rg.*T_d) - pw/sqrt(2*pi*Rg*T_w));
figure(444)
hold on
plot(x_d, -EEE,'k--','linewidth',1.3)
hold on
plot(x_d, E,'k-','linewidth',2)
grid on
ylabel('E [kg/m^2s]')
xlabel('z [m]')
legend('E_{model} (c_{stick} = 1)','E_{Ingersoll}')
%% E vs E (C51 nuc + wall)
% c_stick = Buch (nr. 3), r_sub = 3 (H&H)
[wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, 'xxx', 1, 2000, 1, 3, T_w, c_d, r_sub); % [kg/s]
%%
figure(444)
plot(x_d, wallMassFlux_d./dx_d,'b','linewidth',1.4)
hold on
plot(x_d, accretion_d./dx_d,'r','linewidth',1.4)
plot(x_d, sublimation_d./dx_d,'g','linewidth',1.4)
grid on
ylabel('Wall mass flow [kg/ms]')
xlabel('z [m]')
legend('Total','Accretion','Sublimation')
%% friction & convection
[wallMassFlux_d E m_sigma c_stick accretion_d sublimation_d] = getAccretion(V_d, Rg, T_d, p_d, r_d,rho_d, A_d, dx_d, crackLength, x_d, 'xxx', 1, 2000, 1, 3, T_w, c_d, r_sub); % [kg/s]
[tau tau_d] = wallStress(1, 1, rho_d, V_d, T_d, r_d, rho_ref, sqrt(Rg*y*T_ref)); % T
%[friction] = getFriction(tau, r, L, A_d, throat, c, dx); % non-dimensional friction
[friction_d] = getFriction(tau_d, r_d, L, A_d, throat, c_d, dx_d); % non-dimensional friction
[h Re, Nu Pr] = getReNuPr(T_d,rho_d,V_d,r_d);    
%wallHeatFlux_d = tau_d./V_d.*cp.*(T_d-T_w).*c_d.*dx_d;% convection
wallHeatFlux_d = h.*(T_d-T_w).*c_d.*dx_d;% convection

wallEnergyLoss2 = accretion_d.*(cv.*T_d + V_d.^2./2 + p_d./rho_d);%./L); % [non-dimensional]
wallMomFlux_d = accretion_d.*V_d; % [non-dimensional]
mom_tot = rho_d.*V_d.^2.*A_d;
e_tot = rho_d.*(cv.*T_d+V_d.^2./2).*A_d.*V_d;
   %%
figure(444)
subplot(2,1,1)
semilogy(x_d, abs(friction_d),'r--','linewidth',2.5)
hold on
semilogy(x_d, wallMomFlux_d,'r','linewidth',0.8)
semilogy(x_d, mom_tot,'r','linewidth',1.5)
grid on
ylabel('Momentum flux')
xlabel('z [m]')
legend('Friction','Wall mass flow','Total momentum')
subplot(2,1,2)
semilogy(x_d, abs(wallHeatFlux_d),'b--','linewidth',2.5)
hold on
semilogy(x_d, wallEnergyLoss2, 'b', 'linewidth',0.8) 
semilogy(x_d, e_tot, 'b', 'linewidth', 1.5) 
grid on
%ylim([10*min(wallHeatFlux_d) max(wallEnergyLoss2) ]) 
ylabel('Energy flux')
xlabel('z [m]')
legend('Convection','Wall mass flow','Total energy')
%%
Vth = 1.5.*sqrt(T_d./100).*(1./sqrt(18)).*1000; % [m/s]
figure(444)
plot(x_d, Vth,'b','linewidth',1.4)
hold on
plot(x_d, V_d,'r','linewidth',1.4)
grid on
ylabel('Velocity [m/s]')
xlabel('z [m]')
legend('Thermal','Flow')
%% f, T, M
    figure(123)
if baseline == true
    
    i = 3;
    j = 1.25;
    subplot(i,1,1)
    hold on
    plot(x_d, f, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    hold on
    plot(x_d, M, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    
else
    i = 3;
    j = 1.5;
    subplot(i,1,1)
    hold on
    plot(x_d, f, 'b','linewidth',j)
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    hold on
    plot(x_d, M, 'g','linewidth',j)
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    
end

%% M,T
    figure(123)
if baseline == true
    
    i = 2;
    j = 1.25;
    subplot(i,1,1)
    hold on
    plot(x_d, M, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
else
        i = 2;
    j = 1.5;
    subplot(i,1,1)
    hold on
    plot(x_d, M, 'g','linewidth',j)
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
end
%% T, rho, drho, dRdx, re
i=5; % number of rows

figure(321)
if baseline == true
    j = 1.25; % line width
    subplot(i,1,1)
    plot(x_d, T_d, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, rho_d, 'm--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    hold on
    plot(x_d, drho, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,4)
    semilogy(x_d, dRdx, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,5)
    semilogx(r_e*10^6, P, 'k--','linewidth',j)
    hold on
    %xline(x_d_throat,'--k')
    ylabel('Particle density [m^{-3}]')
    xlabel('Particle radius [\mu m]')
    grid on
else
    j=1.5;
    subplot(i,1,1)
    plot(x_d, T_d, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, rho_d, 'm','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    hold on
    plot(x_d, drho, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,4)
    semilogy(x_d, dRdx, 'g','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,5)
    semilogx(r_e*10^6, P, 'k','linewidth',j)
    hold on
    %xline(x_d_throat,'--k')
    ylabel('Particle density [m^{-3}]')
    xlabel('Particle radius [\mu m]')
    grid on
end
%% S ynuc f %
    figure(123)
if baseline == true
    
    i = 3;
    j = 1.25;
    subplot(i,1,1)
    plot(x_d, S, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, gn, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    plot(x_d, f, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
else
        i = 3;
    j = 1.5;
    subplot(i,1,1)
    plot(x_d, S, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
    subplot(i,1,2)
    plot(x_d, gn, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
    subplot(i,1,3)
    plot(x_d, f, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
end
%% dRho,S %%
figure(1234)
if baseline == true
    
    i = 4;
    j = 1.25;
    subplot(i,1,1)
    hold on
    plot(x_d, drho, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    hold on
    semilogy(x_d, dRdx, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    plot(x_d, f, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,4)
    semilogx(r_e*10^6, P, 'g--','linewidth',j)
    hold on
    %xline(x_d_throat,'--k')
    ylabel('Particle density [m^{-3}]')
    xlabel('Particle radius [\mu m]')
    grid on
else
    i = 4;
    j = 1.5;
    subplot(i,1,1)
    hold on
    plot(x_d, drho, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho - \rho_{eq} [kg/m^3]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
    subplot(i,1,2)
    hold on
    semilogy(x_d, dRdx, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('dR/dx [-]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
    subplot(i,1,3)
    plot(x_d, f, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('f [-]')
    xlabel('z [m]')
    grid on
    xlim([0 150])
    subplot(i,1,4)
    
    semilogx(r_e*10^6, P, 'g','linewidth',j)
    hold on
    %xline(x_d_throat,'--k')
    ylabel('Particle density [m^{-3}]')
    xlabel('Particle radius [\mu m]')
    grid on
    xlim([0 150])
end
    
    
%% MT rho S ynuc %%

    figure(123)
if baseline == true
    
    i = 5;
    j = 1.25;
    subplot(i,1,1)
    hold on
    plot(x_d, M, 'g--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    plot(x_d, rho_d, 'm--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,4)
    plot(x_d, S, 'r--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,5)
    plot(x_d, gn, 'b--','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    grid on
else
        i = 5;
    j = 1.5;
    subplot(i,1,1)
    hold on
    plot(x_d, M, 'g','linewidth',j)
    xline(x_d_throat,'--k')
    ylabel('M [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,2)
    plot(x_d, T_d, 'k','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('T [K]')
    xlabel('z [m]')
    grid on
    subplot(i,1,3)
    plot(x_d, rho_d, 'm','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('\rho [kg/m^3]')
    xlabel('z [m]')
    grid on
    subplot(i,1,4)
    plot(x_d, S, 'r','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('S [-]')
    xlabel('z [m]')
    grid on
    subplot(i,1,5)
    plot(x_d, gn, 'b','linewidth',j)
    hold on
    xline(x_d_throat,'--k')
    ylabel('y_{nuc} [1/m^3s]')
    xlabel('z [m]')
    grid on
end
%% R',P(r),f %%
close all
i=3; % number of rows
j = 1.5; % line width
figure(321)
subplot(i,1,1)
semilogy(x_d, dRdx, 'g','linewidth',j)
hold on
xline(x_d_throat,'--k')
ylabel('dR/dx [-]')
xlabel('z [m]')
grid on
subplot(i,1,2)
semilogx(r_e*10^6, P, 'k','linewidth',j)
hold on
%xline(x_d_throat,'--k')
ylabel('Particle density [m^{-3}]')
xlabel('Particle radius [\mu m]')
grid on
subplot(i,1,3)
plot(x_d, f, 'm','linewidth',j)
hold on
xline(x_d_throat,'--k')
ylabel('f [-]')
xlabel('z [m]')
grid on
%% TSfV %%
close all
figure(123)
subplot(5,1,1)
plot(x_d, T_d, 'g','linewidth',2)
hold on
xline(111,'--k')
ylabel('T [K]')
xlabel('z [m]')
grid on
subplot(5,1,2)
plot(x_d, S, 'k','linewidth',2)
hold on
xline(111,'--k')
ylabel('S [-]')
xlabel('z [m]')
grid on
subplot(5,1,3)
plot(x_d, gn, 'm','linewidth',2)
hold on
xline(111,'--k')
ylabel('\gamma_{nuc} [1/m^3s]')
xlabel('z [m]')
grid on
subplot(5,1,4)
plot(x_d, f, 'r','linewidth',2)
hold on
xline(111,'--k')
ylabel('f [-]')
xlabel('z [m]')
grid on
subplot(5,1,5)
plot(x_d, V_d, 'b','linewidth',2)
hold on
xline(111,'--k')
ylabel('V [m/s]')
xlabel('z [m]')
grid on
%%
figure(856)
j = 1.5;
z = throat;
z = 1;
yline(1,'--','linewidth',1)
%plot(x_d,m_d./m_d(z),'linewidth',j)
hold on
plot(x_d, rho_d./rho_d(z),'linewidth',j)
plot(x_d, V_d./V_d(z),'linewidth',j)
plot(x_d, A_d./A_d(z),'linewidth',j)
plot(x_d, M,'linewidth',j)
xline(x_d_throat,'k--')

legend('Mass flow','\rho','V','A','M','throat')
xlim([(x_d_throat -20) (x_d_throat + 20)])
grid on