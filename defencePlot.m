% Defence bars
close all
x = [-1, 0,1,2]./2;
V = [562; 685; 558; 623; 910];
R = [2.8 1.5 4.5 10 8.5];
RR = [0.1 75];
f = [4.1 7 4.4 8.2 12];
ff = [6.9 20];
T = [257 249 258 246 209];
TT = [170 210];
j = 2;
a = x(1);
b = x(end);
figure(1)
X = categorical({'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
X = reordercats(X,{'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
hold on
for i = 1:length(V)
    h=bar(X(i),V(i));
    if V(i) < 950 && V(i) > 350
        set(h,'FaceColor','g');
    elseif V(i) > 950 || V(i) < 350
        set(h,'FaceColor','r');
    end
end
% hb(1,1).FaceColor = 'r';
% hb(2,1).FaceColor = 'm';
% hb(3).FaceColor = 'c';
% set(hb(1),'FaceColor','blue');
% set(hb(2),'FaceColor','red');
hold on
yline(350,'k--','linewidth', j)
hold on
grid on
yline(950, 'k--','linewidth', j)

ylim([350 950])
% xlim([a b])
ylabel('V [m/s]')
figure(2)
X = categorical({'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
X = reordercats(X,{'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
hold on
for i = 1:length(V)
    h=bar(X(i),R(i));
    if R(i) < RR(2) && R(i) > RR(1)
        set(h,'FaceColor','g');
    elseif R(i) > RR(2) || R(i) < RR(1)
        set(h,'FaceColor','r');
    end
end
yline(RR(1),'k--','linewidth', j)
hold on
grid on
yline(RR(2), 'k--','linewidth', j)
ylim([0 RR(2)])
ylabel('R [\mu m]')

figure(3)
X = categorical({'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
X = reordercats(X,{'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
hold on
for i = 1:length(V)
    h=bar(X(i),f(i));
    if f(i) < ff(2) && f(i) > ff(1)
        set(h,'FaceColor','g');
    elseif f(i) > ff(2) || f(i) < ff(1)
        set(h,'FaceColor','r');
    end
end
yline(ff(1),'k--','linewidth', j)
hold on
grid on
yline(ff(2), 'k--','linewidth', j)
ylim([0 ff(2)])
ylabel('f [%]')

figure(4)
X = categorical({'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
X = reordercats(X,{'Baseline','Expansion ratio \uparrow','Throat \leftarrow','Wall interactions','Expansion ratio = 4'});
hold on
for i = 1:length(V)
    h=bar(X(i),T(i));
    if T(i) < TT(2) && T(i) > TT(1)
        set(h,'FaceColor','g');
    elseif T(i) > TT(2) || T(i) < TT(1)
        set(h,'FaceColor','r');
    end
end
yline(TT(1),'k--','linewidth', j)
hold on
grid on
yline(TT(2), 'k--','linewidth', j)
ylim([170 260])
ylabel('T [K]')
