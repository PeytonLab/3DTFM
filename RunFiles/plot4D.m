%% 4D Plots
x = [-1 0 1]; %X = 0.5 to 2 mM RGD
y = [1 0 -1]; %Y = 5wt to 20 wt
z = [-1 0 1]; %Z = 0 to 50 deg
[X,Y,Z] = meshgrid(x,y,z);

%% Displacement Phase Space
mpm = 0.502709629;  %0.5 mM RGD, 20wt, 0% MMP
ppm = 0.309132485; %2.0 mM RGD, 20wt, 0% MMP
mpp = 0.588873026; %0.5 mM RGD, 20wt, 50% MMP
ppp = 0.160470821; %2.0 mM RGD, 20wt, 50% MMP
ctr = 0.46935435; %1.0 mM RGD, 10wt, 25% MMP
mmm = 0.873686098; %0.5 mM RGD, 5wt, 0% MMP
pmm = 0.858414628; %2.0 mM RGD, 5wt, 0% MMP
mmp = 3.995620409; %0.5 mM RGD, 5wt, 50% MMP
pmp = 1.832755374; %2.0 mM RGD, 5wt, 50% MMP

D(:,:,1) = [mpm, NaN, ppm;
            NaN, NaN, NaN;
            mmm, NaN, pmm];    
D(:,:,2) = [NaN, NaN, NaN;
            NaN, ctr, NaN;
            NaN, NaN, NaN];        
D(:,:,3) = [mpp, NaN, ppp;
            NaN, NaN, NaN;
            mmp, NaN, pmp];          
[XI,YI,ZI] = meshgrid(-1:0.05:1);

idx = ~(isnan(X) | isnan(Y) | isnan(Z) | isnan(D));
DI = griddata(X(idx),Y(idx),Z(idx),D(idx), XI, YI, ZI);
bottom = DI(:,:,1);
top = DI(:,:,end);
front = squeeze(DI(1,:,:))';
back = squeeze(DI(end,:,:))';
left = squeeze(DI(:,1,:))';
right = squeeze(DI(:,end,:))';
figure
surface([-1 1; -1 1],[-1 -1; 1 1],[-1 -1; -1 -1],'FaceColor','texturemap','CData',bottom);
surface([-1 1; -1 1],[-1 -1; 1 1],[1 1; 1 1],'FaceColor','texturemap','CData',top);
surface([-1 1; -1 1],[-1 -1; -1 -1],[-1 -1; 1 1],'FaceColor','texturemap','CData',front);
surface([-1 1; -1 1],[1 1; 1 1],[-1 -1; 1 1],'FaceColor','texturemap','CData',back);
surface([-1 -1; -1 -1],[-1 1; -1 1], [-1 -1; 1 1],'FaceColor','texturemap','CData',left);
surface([1 1; 1 1],[-1 1; -1 1],[-1 -1; 1 1],'FaceColor','texturemap','CData',right);
axis image
set(gcf,'color','w')
colormap(jet)
c = colorbar;
ylabel(c, 'Velocity, {\mu}m/hr','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})
zlabel('Degradability','FontWeight','Bold')
zticks([-1 1])
zticklabels({'0%','50%'})

%3D Cross sections of MMP axis
%0 MMP
surf(XI(:,:,1),YI(:,:,1),DI(:,:,1),'edgecolor','interp')
axis image
set(gcf,'color','w')
colormap(jet)
c = colorbar;
ylabel(c, 'Velocity, {\mu}m/hr','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})

%50 MMP
surf(XI(:,:,end),YI(:,:,end),DI(:,:,end),'edgecolor','interp')
axis image
set(gcf,'color','w')
colormap(jet)
c = colorbar;
ylabel(c, 'Velocity, {\mu}m/hr','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})


%% Traction Phase Space
mpm = 38.68832822;  %0.5 mM RGD, 20wt, 0% MMP
ppm = 13.72389451; %2.0 mM RGD, 20wt, 0% MMP
mpp = 682.0670671; %0.5 mM RGD, 20wt, 50% MMP
ppp = 11.08822647; %2.0 mM RGD, 20wt, 50% MMP
ctr = 11.70816081; %1.0 mM RGD, 10wt, 25% MMP
mmm = 1027.463402; %0.5 mM RGD, 5wt, 0% MMP
pmm = 112.6630259; %2.0 mM RGD, 5wt, 0% MMP
mmp = 16642.91742; %0.5 mM RGD, 5wt, 50% MMP
pmp = 14996.00699; %2.0 mM RGD, 5wt, 50% MMP

D(:,:,1) = [mpm, NaN, ppm;
            NaN, NaN, NaN;
            mmm, NaN, pmm];    
D(:,:,2) = [NaN, NaN, NaN;
            NaN, ctr, NaN;
            NaN, NaN, NaN];        
D(:,:,3) = [mpp, NaN, ppp;
            NaN, NaN, NaN;
            mmp, NaN, pmp];          
[XI,YI,ZI] = meshgrid(-1:0.05:1);
idx = ~(isnan(X) | isnan(Y) | isnan(Z) | isnan(D));
DI = griddata(X(idx),Y(idx),Z(idx),D(idx), XI, YI, ZI);
bottom = DI(:,:,1);
top = DI(:,:,end);
front = squeeze(DI(1,:,:))';
back = squeeze(DI(end,:,:))';
left = squeeze(DI(:,1,:))';
right = squeeze(DI(:,end,:))';

figure
surface([-1 1; -1 1],[-1 -1; 1 1],[-1 -1; -1 -1],'FaceColor','texturemap','CData',bottom);
surface([-1 1; -1 1],[-1 -1; 1 1],[1 1; 1 1],'FaceColor','texturemap','CData',top);
surface([-1 1; -1 1],[-1 -1; -1 -1],[-1 -1; 1 1],'FaceColor','texturemap','CData',front);
surface([-1 1; -1 1],[1 1; 1 1],[-1 -1; 1 1],'FaceColor','texturemap','CData',back);
surface([-1 -1; -1 -1],[-1 1; -1 1], [-1 -1; 1 1],'FaceColor','texturemap','CData',left);
surface([1 1; 1 1],[-1 1; -1 1],[-1 -1; 1 1],'FaceColor','texturemap','CData',right);
axis image
set(gcf,'color','w')
colormap(jet)
c = colorbar;
ylabel(c, 'Traction, Pa/s','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})
zlabel('Degradability','FontWeight','Bold')
zticks([-1 1])
zticklabels({'0%','50%'})

%3D Cross sections of MMP axis
%0 MMP
figure
surf(XI(:,:,1),YI(:,:,1),DI(:,:,1),'edgecolor','interp')
axis image
set(gcf,'color','w')
set(gca, 'ZScale', 'log')
colormap(jet)
c = colorbar;
ylabel(c, 'Traction, Pa/s','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})

%50 MMP
figure
surf(XI(:,:,end),YI(:,:,end),log(DI(:,:,end)),'edgecolor','interp')
axis image
set(gcf,'color','w')
colormap(jet)
c = colorbar;
ylabel(c, 'log(Traction), log(Pa/s)','FontWeight','Bold')
xlabel('Adhesivity','FontWeight','Bold')
set(gca, 'ZScale', 'log')
xticks([-1,1])
xticklabels({'0.5mM','2mM'})
ylabel('Stiffness','FontWeight','Bold')
yticks([-1 1])
yticklabels({'1.25 kPa','3.75 kPa'})