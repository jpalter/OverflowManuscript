

%Claude Sonnet 4.5 Code (11/24/25) for downloading monthly mean
%climatological SST and mapping. I had to do a bunch of tweaking to make it
%work, including manual download.  The only thing left from the original
%code was some of the mapping commands. NOT A TIME SAVER This time


%downloaded the 1981-2021 long term mean sst from:
% https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
%NOAA OISST v2 High Resolution page:
%https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
%this is the 0.25 degree climatology 1991-2020;

%% 
%ds = netcdf_load('sst.mon.ltm.1981-2010.nc');
ds = netcdf_load('sst.mon.ltm.1991-2020.nc');

% Read metadata
lon =double(ds.lon);
lat = double(ds.lat);
time = ds.time;
sst = double(ds.sst);
sst(sst==min(sst(:))) =NaN;

%sst_annual = mean(sst,3);
sst = mean(cat(3,sst(:,:,12),sst(:,:,1:2)),3);


[lati loni] = meshgrid(lat,lon);

% Convert longitude from 0-360 to -180-180
% Find where to split (at 180°)
split_idx = find(lon >= 180, 1);

% Rearrange: put 180-360 at beginning, 0-180 at end
lon = [lon(split_idx:end) - 360; lon(1:split_idx-1)];
%sst = [sst_annual(split_idx:end, :); sst_annual(1:split_idx-1, :)];
sst = [sst(split_idx:end, :); sst(1:split_idx-1, :)];

% Extract North Atlantic indices
lon_idx = find(lon >= -80 & lon <= 20);
lat_idx = find(lat >= 20 & lat <= 85);


% Read the SST data

% Calculate annual mean (average over all months)
%sst_annual = mean(sst, 3);

% Extract coordinates
lon_natl = lon(lon_idx);
lat_natl = lat(lat_idx);


%load etop02 bathymetry
bathy = netcdf_load('MATLAB_PATH_FUNCTIONS/ETOPO2v2g_f4.nc');
bathy_lon = double(bathy.x);
bathy_lat = double(bathy.y);
D = double(bathy.z);
% Extract North Atlantic indices
lon_idxb = find(bathy_lon >= -80 & bathy_lon <= 20);
lat_idxb = find(bathy_lat >= 20 & bathy_lat <= 85);
lonb_natl = bathy_lon(lon_idxb);
latb_natl = bathy_lat(lat_idxb);
D_natl = D(lon_idxb,lat_idxb);
[lab lob] = meshgrid(latb_natl,lonb_natl);


% Create map
figure('Position', [100, 100, 1000, 600]);


axesm('MapProjection', 'miller', 'MapLatLimit', [35 80], ...
    'MapLonLimit', [-70 20])
[c h] = contourfm(lat_natl, lon_natl, sst(lon_idx,lat_idx)',-2:2:22,'linestyle','none');
colorbar
cmocean('thermal',length(-2:2:22)-1);
hold on
contourm(lab,lob,D_natl,[-600 -600],'w','LineWidth',1)


geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8])
caxis([-2 22])
gridm on
title('DJF Mean SST (1991-2020)')
tightmap
mlabel('south');plabel

% Set figure properties for clean export
set(gcf, 'Renderer', 'painters')
set(gcf, 'Color', 'white')

% Export as EPS
print('SST_map.eps', '-depsc', '-painters', '-r300')

%% 

%Use this same m-file to make the surface air temperature anomalies from
%global mean

sat = netcdf_load('NOAA_20thcenturyreanalysis_SAT_DJFclim19812010.nc');
zonal_mean_sat = mean(sat.var,1);
sat_anom = double(sat.var-zonal_mean_sat);
[lati loni] = meshgrid(double(sat.lat),double(sat.lon));

figure;
levels= [-25:5:25];
axesm('miller', ...
      'MapLatLimit', [0 85], 'MapLonLimit', [-180 180],...
      'Frame','on','Grid','on');
[c h] =contourfm(lati,loni,sat_anom,levels,'linestyle','none'); 
tightmap; colorbar;
geoshow('landareas.shp', ...
'FaceColor', 'none', ...   % No fill
    'EdgeColor', 'k', ...      % Outline color (black)
    'LineWidth', 1);

cmocean('balance',numel(levels),'pivot',0,'color','w','BackgroundColor','none');
plabel
mlabel
set(gca,'fontsize',14)
[t] = clabelm(c,h,'manual');
set(t,'color','w','backgroundcolor','none')
cmocean('balance',numel(levels),'pivot',0,'color','w','BackgroundColor','none');