% SFTDATA -- data-assimilative flux transport model.
%
% - A.R. Yeates, Durham University 8/3/16

close all; clear all;

%% (0) User-defined parameters:
% Type of magnetogram and path to magnetogram files:
magType = 'kp';
magPath = '~/Desktop/kpmag/';
% Path to output (if any):
outPath = '~/Desktop/2016_03_08_sft1/';
% Number of cells in simulation grid in sin(latitude) and longitude:
ns = 90;
np = 180;
% Carrington rotations to start and stop at:
rot0=1910;
rot1=1912;
% Output plots of Br to screen, once per rotation?
showplots=0;
% Save aforementioned plots to outPath as pngs?
saveplots=0;

%% (1) Read initial synoptic map
[br0 pc sc flux] = readSynoptic(magType, magPath, rot0, ns, np);
ns=size(br0,1);
np=size(br0,2);

%% (2) Precompute coordinate arrays
% ( i.e. equally spaced cells in longitude and sin(latitude) )
latc=asin(sc);
[junk, sthc]=meshgrid(pc,cos(asin(sc)));
sthcth=cat(2,sthc,sthc(:,1));
sg=2./(ns)*((1:(ns+1))-1)-1;
[junk, sthg]=meshgrid(pc,cos(asin(sg)));
[junk, slat]=meshgrid(pc,sc);
slatcth=cat(2,slat,slat(:,1));
latg=asin(sg);
[junk, latph]=meshgrid(pc,latg);
dp=pc(2)-pc(1);
ds=sc(2)-sc(1);

%% (3) Compute vector potential on ribs from Bz at cell centres:
% Br = 1/(sin(th))*( d/dth(sin(th)Aph) - d/dph(Ath) )
% Aph(th,ph) = 0
% Ath(th,ph) = -sin(th)*int_0^ph Br(ph',th) dph'
aph=zeros(ns+1,np);
ath=-sthc.*cumsum(br0,2)*dp;
ath=cat(2,zeros(ns,1),ath);

%% (4) Initialise flows and diffusion
% (a) Supergranular diffusion coefficient:
eta=500./6.96e5^2;
% (b) Differential rotation (on th ribs):
omA=(13.7195 - 360/27.2753)/180*pi/86400.;
omB=-2.396/180*pi/86400.;
omC=-1.787/180*pi/86400.;
om=omA + omB*slatcth.^2 + omC*slatcth.^4;
vph_om = om.*sthcth;
% (c) Meridional flow (on ph ribs):
% -Schuessler-Baumann profile
v0=11./6.96e8;
vth_mf = -sin(2*latph).*exp(pi*(1.-2*abs(latph)/pi));
vth_mf = vth_mf*v0/max(abs(vth_mf(:)));

%% (5) Choose timestep
% Timesteps (for CFL condition):
hphmin=min(abs(sthc(:)*dp));
hthmin=min(latg(2:end)-latg(1:end-1));
dt_eta=min(hphmin^2/eta,hthmin^2/eta);
%
t_mf=abs((latph(2:end,:) - latph(1:end-1,:))./vth_mf(1:end-1,:));
dt_mf=min(t_mf(:));
%
dt_om=min(abs(sthcth(:)*dp./vph_om(:)));
%
dt=min([dt_eta,dt_mf,dt_om])*0.4;
% Modify to fit exactly into one Carrington rotation:
ndt=round(27.2753*86400/dt);
dtday=27.2753/ndt;
dt=dtday*86400;

%% (6) Main time integration
% Color table:
load('Bluered.mat');
% Set up plotting window if required:
if (showplots==1)
    f1=figure(1);
    set(gcf,'Units','centimeters','Position',[10 10 16 16], 'PaperPositionMode', 'auto');
    cmax=10;
end  
% Create output directory (if necessary):
system(['mkdir ' outPath]);
% Initialize arrays:
brc=zeros(ns+2,np+2);
nrot=rot1-rot0;
bfly=zeros(ns,nrot);
unflux=zeros(1,nrot);
%
% Precompute as much as possible to speed execution:
f1 = eta./sthcth/dp;
f2 = eta*sthg/ds;
f3 = sthg/ds;
f4 = dp*sthc;
vph_om = vph_om*0.5;
vth_mf = vth_mf*0.5;
%
% Evolve for one Carrington rotation at a time:
for n=1:nrot
    disp(sprintf('Rotation %g',rot0+n));
    %
    % (a) Generate list of emerging regions for this rotation, and their times:
    [map1 lm reglon]=getEmergingRegions(magType, magPath, rot0+n, ns, np);
    nlm=max(lm(:));
    % Convert Carrington longitude to time (in days):
    regt=(360-reglon)*27.2753/360.;
    % Round to nearest day:
    regt=round(regt);
    %
    % (b) Sort into day order:
    [regt isort]=sort(regt,'descend');
    lm0=lm;
    for i=1:nlm
       lm(lm0==i)=isort(i); 
    end
    day=0.;
    nxtreg=nlm;
    if (nxtreg > 0)
        nxtregt = regt(nxtreg);
    else
        nxtregt = -1;
    end
    %
    tic
    % (c) Loop over individual timesteps within the rotation:
    for step=1:ndt
        % (i) Compute B at cell centres from A
        % Note d/dth = d(cos(th))/dth*d/d(cos(th)) = -sin(th)d/d(cos(th))
        brc(2:end-1,2:end-1) = ...
            f3(1:end-1,:).*aph(1:end-1,:) - f3(2:end,:).*aph(2:end,:) ...
            +(ath(:,1:end-1) - ath(:,2:end))./f4;
        %
        % (ii) Add in new regions, if any:
        if (nxtregt==floor(day))
            while (nxtregt == floor(day))
                br1=brc(2:end-1,2:end-1);
                npts=sum(lm(:)==nxtreg);
                oldflux=sum(br1(lm==nxtreg)*dp*ds);
                br1(lm==nxtreg) = map1(lm==nxtreg) + oldflux/(npts*dp*ds);                
                brc(2:end-1,2:end-1)=br1;
                % Recompute vector potential:
                aph=zeros(ns+1,np);
                ath=-f4.*cumsum(brc(2:end-1,2:end-1),2);
                ath=cat(2,zeros(ns,1),ath);
                nxtreg=nxtreg-1;
                if (nxtreg > 0)
                    nxtregt = regt(nxtreg);
                else
                    nxtregt = -1;
                end
            end
        end
        %
        % (iii) Apply boundary conditions (periodic in ph):
        brc(2:end-1,1)=brc(2:end-1,end-1);
        brc(2:end-1,end)=brc(2:end-1,2);
        %
        % (v) Compute emf on ribs from differential rotation
        emfth = -vph_om.*(brc(2:end-1,1:end-1) + brc(2:end-1,2:end));
        %
        % (vi) Compute emf on ribs from meridional flow
        emfph = vth_mf.*(brc(1:end-1,2:end-1) + brc(2:end,2:end-1));
        %
        % (vii) Compute emf on ribs from supergranular diffusion
        % jth = 1/(sin(th))*d/dph(Br)
        % jph = -d/dth(Br) = sin(th)*d/d(cos(th))(Br)
        emfth = emfth + f1.*(brc(2:end-1,2:end)-brc(2:end-1,1:end-1));
        emfph = emfph + f2.*(brc(2:end,2:end-1)-brc(1:end-1,2:end-1));
        %
        % (viii) Update vector potential
        ath = ath - dt*emfth;
        aph = aph - dt*emfph;
        %
        day = day + dtday;
    end
    %
    % (d) Add to butterfly diagram array:
    bfly(:,n) = mean(brc(2:end-1,2:end-1),2);
    %
    % (e) Add to unsigned flux array:
    unflux(n) = sum(sum(abs(brc(2:end-1,2:end-1))*dp*ds));
    %
    if (showplots==1)
        figure(1);
        % (f) Plot Br in lat-long:
        % Simulated:
        subplot(2,1,1);
        colormap(cmap);
        h=pcolor(pc,latc,brc(2:end-1,2:end-1));
        set(h,'EdgeColor','none');
        caxis([-cmax, cmax]);
        colorbar;
        title(sprintf('End of CR %g',rot0+n));
        xlabel('Longitude'); ylabel('Latitude');
        % Observed:
        subplot(2,1,2);
        map1 = readSynoptic(magType, magPath, rot0+n, ns, np, 0);
        % [Correct flux imbalance:]
        fnet=sum(map1(:)*dp*ds);
        map1=map1-fnet/(ns*np*dp*ds);
        h=pcolor(pc,latc,map1);   
        set(h,'EdgeColor','none');
        caxis([-cmax, cmax]);
        colorbar;
        title(sprintf('Observed synoptic CR %g',rot0+n));
        xlabel('Longitude'); ylabel('Latitude');
        pause(0.001);
        % Save the plot to file:
        if (saveplots==1)
            saveas(gcf(),strcat(outPath,sprintf('br_cr%4.4d.png',rot0+n)));
        end
    end
end

%% (7) Read in observed synoptic data for comparison:
bfly_obs=zeros(ns,nrot);
unflux_obs=zeros(1,nrot);
for n=1:nrot
    map1 = readSynoptic(magType, magPath, rot0+n, ns, np, 0);
    %
    % Smooth observed magnetogram:    
    myfilt = fspecial('gaussian',[13 13], 5);
    map1sm=imfilter(map1, myfilt, 'replicate');
    % Add column to butterfly diagram:
    bfly_obs(:,n) = mean(map1sm,2);
    % Compute unsigned flux
    unflux_obs(n) = sum(sum(abs(map1)*dp*ds));
end

%% (8) Plot the butterfly diagrams:
figure();
colormap(cmap);
set(gcf,'Units','centimeters','Position',[10 10 16 16], 'PaperPositionMode', 'auto');
cmax=6;
% Simulated:
subplot(2,1,1);
h=pcolor((rot0+1):rot1,latc,bfly);
set(h,'EdgeColor','none');
caxis([-cmax, cmax]);
colorbar;
title('SIMULATION');
xlabel('Carrington Rotation');
ylabel('Latitude');
%Observed:
subplot(2,1,2);
h=pcolor((rot0+1):rot1,latc,bfly_obs);
set(h,'EdgeColor','none');
caxis([-cmax, cmax]);
colorbar;
title('INPUT MAGNETOGRAMS');
xlabel('Carrington Rotation');
ylabel('Latitude');
saveas(gcf(),strcat(outPath,'bfly.png'));

%% (9) Plot the total unsigned fluxes:
% (note: there is an offset between the two because of the smoothing out of
% small-scale flux in the simulation)
figure();
plot((rot0+1):rot1,unflux*6.96e10^2,'b');
hold on;
plot((rot0+1):rot1,unflux_obs*6.96e10^2,'r');
plot((rot0+1):rot1,unflux*6.96e10^2,'sb');
plot((rot0+1):rot1,unflux_obs*6.96e10^2,'ro');
grid on;
xlabel('Carrington Rotation');
ylabel('Total Unsigned Flux (Mx)');
legend('SIMULATED','MAG','Location','NW');
saveas(gcf(),strcat(outPath,'unflux.png'));
