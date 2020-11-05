%Yoshida Method
clear;clc;close all;

[r0,v0,m] = solarSystemData;
%Remove Pluto from data
r0 = r0(:,1:end-1);
v0 = v0(:,1:end-1);
m = m(1:end-1);

%Initial data
nBodies = length(m);
t0 = 0;     % Days
tEnd = 200*365.25; % 200 years --> Days
dt = 0.1;     % Days

%Integrator
[rVec, vVec] = yoshida4Integrator(nBodies, r0, v0, m, t0, tEnd, dt);

%Plot Integrator data
figure(1)
for index = 1:9
    plot3(rVec{1,index}(1,:),rVec{1,index}(2,:),rVec{1,index}(3,:))
    hold on
end
xlim([-32,32])
ylim([-32,32])
zlim([-32,32])
set(gca,'FontSize',12)
xlabel('[au]','FontSize',16,'Interpreter','latex')
ylabel('[au]','FontSize',16,'Interpreter','latex')
zlabel('[au]','FontSize',16,'Interpreter','latex')
legend({'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'},'Interpreter','latex','Location','northeast','FontSize',14)
hold off
%Save plot
savedir = '/Users/colincoane/Documents/College Work/Fall 2020/PHYS 304/Project/Code';
saveas(figure(1),fullfile(savedir,'yoshida_plot.png'));

%Plot Closeup data
figure(2)
for index = 1:6
    plot3(rVec{1,index}(1,:),rVec{1,index}(2,:),rVec{1,index}(3,:))
    hold on
end
xlim([-6.5,6.5])
ylim([-6.5,6.5])
zlim([-2,2])
set(gca,'FontSize',12)
xlabel('[au]','FontSize',16,'Interpreter','latex')
ylabel('[au]','FontSize',16,'Interpreter','latex')
zlabel('[au]','FontSize',16,'Interpreter','latex')
legend({'Sun','Mercury','Venus','Earth','Mars','Jupiter'},'Interpreter','latex','Location','northeast','FontSize',14)
hold off
%Save plot
savedir = '/Users/colincoane/Documents/College Work/Fall 2020/PHYS 304/Project/Code';
saveas(figure(2),fullfile(savedir,'yoshida_closeup.png'));

%Capture end of ephemeris data
[r0End,v0End] = solarSystemEndData;
r0End = r0End(:,1:end-1);
v0End = v0End(:,1:end-1);

%Capture end of ephemeris data from integrator
rEnd = zeros(3,nBodies);
vEnd = zeros(3,nBodies);
for index = 1:nBodies
    rEnd(:,index) = rVec{1,index}(:,end);
    vEnd(:,index) = vVec{1,index}(:,end);
end

r0E = zeros(1,nBodies);
v0E = zeros(1,nBodies);
rE = zeros(1,nBodies);
vE = zeros(1,nBodies);

%Scalar positions, velocities
for index = 1:nBodies
    r0E(1,index) = sqrt(r0End(1,index)^2 + r0End(2,index)^2 + r0End(3,index)^2);
    v0E(1,index) = sqrt(v0End(1,index)^2 + v0End(2,index)^2 + v0End(3,index)^2);
    rE(1,index) = sqrt(rEnd(1,index)^2 + rEnd(2,index)^2 + rEnd(3,index)^2);
    vE(1,index) = sqrt(vEnd(1,index)^2 + vEnd(2,index)^2 + vEnd(3,index)^2);
end

%Error Analysis (of scalar positions velocities)
rErr = zeros(1,nBodies);
vErr = zeros(1,nBodies);

for idx1 = 1:nBodies
    rErr(idx1) = abs((rE(idx1) - r0E(idx1))/r0E(idx1))*100;
    vErr(idx1) = abs((vE(idx1) - v0E(idx1))/v0E(idx1))*100;
end

disp(rErr)
disp(vErr)