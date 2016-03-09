function [map1 lm reglon] = getEmergingRegions(magType, magPath, rot, ns, np, plots)
% GETEMERGINGREGIONS.m Determine strong flux regions in an observed synoptic map.
%
% INPUTS:
% magType -- name of observatory to use (affects filename and how data read)
%      -- possibilities are 'kp' for Kitt Peak
% magPath -- path and directory where magnetograms can be found
% rot0 -- required carrington rotation (int)
% ns, np -- number of grid cells in sin(lat) and long respectively (ints)
% [OPTIONAL] plots -- 0 or 1; if 1 then plot each new region
%
% OUTPUTS:
% map1 -- 
% lm -- 
% reglon -- 
%
% - A.R. Yeates, Durham University 8/3/16

if (nargin < 6)
    plots = 0;
end

% Grid cell sizes:
dp = 2*pi/np;
ds = 2/ns;

%% Read in observed synoptic map and correct flux imbalance:
[map1 pc sc flux] = readSynoptic(magType, magPath, rot, ns, np);
[pc2 sc2]=meshgrid(pc,sc);

%% Identify strong flux regions by smoothing absolute flux then contouring:
% Absolute value of synoptic map:
map=abs(map1);
% Smooth by convolving with gaussian filter:
myfilt = fspecial('gaussian',[12 12], 3);
map=imfilter(map, myfilt, 'replicate');
% Cut-off in flux:
maski=(map > 15);
mask=maski*1.;
% Identify and number separate regions:
cc=bwconncomp(mask,8);
lm=labelmatrix(cc);
nlm=max(lm(:));
% Correct each of the selected regions for flux balance:
unb=[];
for i=1:nlm
   npts=sum(lm(:)==i);
   freg=sum(map1(lm==i)*ds*dp);
   fabs=sum(abs(map1(lm==i))*ds*dp);
   if (abs(freg)/abs(fabs) > 0.5)
       unb=[unb,i];
   else
       map1(lm==i)=map1(lm==i) - freg/(npts*ds*dp);
   end
end
% Remove regions with very unbalanced flux from list:
if (length(unb)>0)
   for i=1:length(unb)
       irm=unb(i);
       lm(lm==irm)=0;
       lm(lm>irm)=lm(lm>irm)-1;
       nlm=nlm-1;
       unb=unb-1;
   end
end
% Get (area) centroid of each region:
reglon=zeros(1,nlm);
for i=1:nlm
   reglon(i)=mean(pc2(lm==i))*180/pi; 
end

if (plots==1)
    for i=1:nlm
        map2 = map1*0;
        map2(lm==i) = map1(lm==i);
        f4=figure(4);
        set(f4,'Units','centimeters','Position',[10 10 10 10], 'PaperPositionMode', 'auto');
        load('Bluered.mat');
        colormap(cmap);
        imagesc(map2);
        set(gca,'Ydir','normal');
        colorbar;
        title('REGION');
        cmax=50;%max(abs(map2(:)));
        caxis([-cmax,cmax]);
        pause(1)
    end
end

end