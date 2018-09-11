function [bem lem reglon leadpol] = getEmergingRegions(magType, magPath, rot, ns, np)
% GETEMERGINGREGIONS.m Determine strong flux regions in an observed synoptic map.
%
% INPUTS:
% magType -- name of observatory to use (affects filename and how data read)
%      -- possibilities are 'kp' for Kitt Peak
% magPath -- path and directory where magnetograms can be found
% rot0 -- required carrington rotation (int)
% ns, np -- number of grid cells in sin(lat) and long respectively (ints)
%
% OUTPUTS:
% map1(nrg,ns,np) -- Br for each emerging region (nrg of them)
% lm(nrg,ns,np) -- binary mask for each emerging region
% reglon(nrg) -- central longitude of emerging region
% leadpol(nrg) -- sign of leading polarity of each emerging region
%
% - A.R. Yeates, Durham University 30/8/18

if (nargin < 6)
    plots = 0;
end

% Grid cell sizes:
dp = 2*pi/np;
ds = 2/ns;

%% Read in observed synoptic map and neighbours, and correct flux imbalance:
[map1 pc sc flux] = readSynoptic(magType, magPath, rot, ns, np);
[mapl pc sc fluxl] = readSynoptic(magType, magPath, rot+1, ns, np);
[mapr pc sc fluxr] = readSynoptic(magType, magPath, rot-1, ns, np);

[pc2 sc2]=meshgrid(pc,sc);

%% Concatenate maps:
map1 = cat(2, mapl, map1, mapr);

%% Identify strong flux regions by smoothing absolute flux then contouring:
SIG = 3;  % width of smoothing filter
BPAR = 39.8; % cut-off
% Absolute value of synoptic map:
map=abs(map1);
% Smooth by convolving with gaussian filter:
myfilt = fspecial('gaussian',[4*SIG 4*SIG], SIG);
map=imfilter(map, myfilt, 'replicate');
% Cut-off in flux:
maski=(map > BPAR);
mask=maski*1.;
% Identify and number separate regions:
cc=bwconncomp(mask,8);
lm=labelmatrix(cc);
nlm=max(lm(:));

%% Correct each of the selected regions for flux balance:
unb=[];
for i=1:nlm
   npts=sum(lm(:)==i);
   freg=sum(map1(lm==i));
   fabs=sum(abs(map1(lm==i)));
   if (abs(freg)/abs(fabs) > 0.5)
       unb=[unb,i];
   else
       map1(lm==i)=map1(lm==i) - freg/npts;
   end
end

%% Remove regions with very unbalanced flux from list:
if (length(unb)>0)
   for i=1:length(unb)
       irm=unb(i);
       lm(lm==irm)=0;
       lm(lm>irm)=lm(lm>irm)-1;
       nlm=nlm-1;
       unb=unb-1;
   end
end

%% Get longitude centroid of each region:
np3 = np*3;
pc3 = 6*pi/np3*((1:np3) - 0.5);
[pc2 sc2]=meshgrid(pc3,sc);
reglon=zeros(1,nlm);
for i=1:nlm
   reglon(i)=mean(pc2(lm==i)); 
end

%% Remove regions with centroid outside original map:
kp=ones(1,nlm);
outsd=[];
for i=1:nlm
   if (reglon(i)<2*pi | reglon(i)>=4*pi)
       outsd=[outsd,i];
       kp(i)=0;
   end
end
if (length(outsd)>0)
   for i=1:length(outsd)
       irm=outsd(i);
       lm(lm==irm)=0;
       lm(lm>irm)=lm(lm>irm)-1;
       nlm=nlm-1;
       outsd=outsd-1;
   end
end
reglon = reglon(kp==1);
reglon = (reglon - 2*pi)*180/pi;

% Store bmap for each region in array (for single map - account for
% neighbouring maps by rotating):
bem=zeros(1,ns,np);
for j=1:nlm
   bem = cat(1, bem, zeros(1,ns,np));
   btmp = map1.*(lm==j);
   bem(end,:,:) =  btmp(:,(np+1):2*np) + btmp(:,1:np) + btmp(:,(2*np+1):3*np);
end
bem=bem(2:end,:,:);  

% Do same for label map:
lem=zeros(1,ns,np);
for j=1:nlm
   lem = cat(1, lem, zeros(1,ns,np));
   ltmp = lm*0 + 1;
   ltmp(lm ~= j) = 0;
   lem(end,:,:) =  ltmp(:,(np+1):2*np) + ltmp(:,1:np) + ltmp(:,(2*np+1):3*np);
end
lem=lem(2:end,:,:); 

% Determine sign of leading polarity, for each region:
leadpol=ones(1,nlm);
for j=1:nlm
    btmp = map1.*(lm==j);
    btmp(btmp < 0) = 0;
    cenpos=sum(sum(abs(btmp).*pc2))/sum(sum(abs(btmp)));
    btmp = map1.*(lm==j);
    btmp(btmp > 0) = 0;
    cenneg=sum(sum(abs(btmp).*pc2))/sum(sum(abs(btmp)));
    if (cenneg > cenpos)
        leadpol(j) = -1;
    end
end

end
