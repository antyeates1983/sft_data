function [BR pc sc flux] = readSynoptic(type, magPath, rot0, ns, np, plot)
% READSYNOPTIC.m Read in observed synoptic map and remove monopole component by 
% a multiplicative correction.
%
% INPUTS:
% type -- name of observatory to use (affects filename and how data read)
%      -- possibilities are 'kp' for Kitt Peak
% magPath -- path and directory where magnetograms can be found
% rot0 -- required carrington rotation (int)
% ns, np -- number of grid cells in sin(lat) and long respectively (ints)
% [OPTIONAL] plot -- either 0 or 1; if 1 then plots original and remapped magnetogram
%
% OUTPUTS:
% BR -- 2d array of br at cell centres of simulation grid (sinlat, lon)
% pc -- 1d array of longitude at cell centres of simulation grid
% sc -- 1d array of sin(lat) at cell centres of simulation grid
% flux -- total unsigned flux of BR (Mx)
%
% - A.R. Yeates, Durham University 8/3/16

if (nargin < 6)
    plot = 0;
end

%% Form grid arrays (using only a single subgrid):
dp = 2*pi/np;
pc = 0.5*dp:dp:(2*pi - 0.5*dp);
ds = 2/ns;
sc = (-1+0.5*ds):ds:(1-0.5*ds);

% 2d (ph,th) array:
[PHC,SC]=meshgrid(pc,sc);


%% (1) Read in synoptic map:
if (type=='kp')
    fname = strcat(magPath,sprintf('m%4.4df.fits',rot0));
    if exist(fname,'file')
        BR0=fitsread(fname);
    else
        BR0=zeros(ns,np);
    end
    nlat=size(BR0,1);
    nlon=size(BR0,2);
    sc0=2./nlat*((1:nlat) - 0.5) - 1;
    phc0=2.*pi/nlon*((1:nlon) - 0.5);
end
    
% 2d (ph,th) array:
[PHC0,SC0]=meshgrid(phc0,sc0);

%% (2) Remap to computational grid:
BR=interp2(PHC0,SC0,BR0,PHC,SC,'spline');

%% (3) New method of monopole removal (multiplicative):
BR1 = BR;
BR1(BR1 < 0) = 0;
fluxp = abs(sum(BR1(:)))*dp*ds;
BR1 = BR;
BR1(BR1 > 0) = 0;
fluxn = abs(sum(BR1(:)))*dp*ds;
% Rescale either positive or negative field, to match polarity with lower
% flux.
if (fluxp < fluxn)
    BR(BR<0) = BR(BR<0)*fluxp/fluxn;
else
    BR(BR>0) = BR(BR>0)*fluxn/fluxp;
end

%% (4) Unsigned flux (on simulation grid):
flux = sum(abs(BR(:)))*dp*ds;
flux = flux*(6.96e10)^2;

if (plot == 1)
    %% (5) Plot the two side-by-side to double check:
    figure();
    subplot(1,2,1)
    h=pcolor(PHC0*180/pi,SC0,BR0);
    set(h,'EdgeColor','none');
    caxis([-25,25]);
    ylim([-1,1]);
    colorbar;
    ylabel('sin(Latitude)');
    xlabel('Longitude');
    title('Original Grid');
    %
    subplot(1,2,2);
    h=pcolor(PHC*180/pi,SC,BR);
    set(h,'EdgeColor','none');
    caxis([-25,25]);
    ylim([-1,1]);
    colorbar;
    ylabel('sin(Latitude)');
    xlabel('Longitude');
    title('Simulation Grid');
end

end