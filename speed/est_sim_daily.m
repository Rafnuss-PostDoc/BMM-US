

%% Load data
load('data/density/inference-trans.mat','g','radar')
load('data/speed/inference-trans.mat')
load('data/speed/inference-daily.mat')
addpath('./functions/'); 


%% ----- SERVER
% Simulation
gn=15; % number of date to simulation togeter
gnm=4; % margin around those date
sim_n=2;
plotit=false;

uv={'u','v'};

for i_uv = 1:2

    % unconditional realiation
    Resd.(uv{i_uv})=simulationDailyUnconditional(uv{i_uv}, g, gn, sim_n, dcov, plotit); % 15min
    % save('data/speed/simuncond-daily_v','ResdV','dcov','gn','-v7.3')
    
    % Simulation and Estimation
    condVal = victddt.(uv{i_uv}); % conditional value
    [Ed.(uv{i_uv}),Sd.(uv{i_uv}),ResdC.(uv{i_uv})] = estimationSimulationDailyConditional(uv{i_uv}, g, gn, gnm, radar, dcov, condVal, Resd.(uv{i_uv})); % 1h50 
    
    % Transform inverse
    victd_mean_doy_ll_t.(uv{i_uv})=nan(numel(doy_g),sum(~g.mask_water(:)));
    for i_doy = 1:numel(doy_g)
         victd_mean_doy_ll_t.(uv{i_uv})(i_doy,:)=victd_trend.(['sf_mean_' uv{i_uv}]){i_doy}([g.LON(~g.mask_water),g.LAT(~g.mask_water)]);
    end
    victd_mean_doy_ll_t_f.(uv{i_uv}) = interp1(doy_g,victd_mean_doy_ll_t.(uv{i_uv}),1:366,'pchip')';
    
    f_trans_inver = @(X)  ( X ... 
        .* victd_trend.(['t_std_' uv{i_uv}])(g.day_doy)' ...
        .* victd_trend.(['sf_std_' uv{i_uv}])([g.LON(~g.mask_water),g.LAT(~g.mask_water)]) )...
        + victd_mean_doy_ll_t_f.(uv{i_uv})(:,g.day_doy) ;

    Edm.(uv{i_uv}) = f_trans_inver(Ed.(uv{i_uv}));
    Sdm.(uv{i_uv}) = Sd.(uv{i_uv}) .* abs(victd_trend.(['sf_std_' uv{i_uv}])([g.LON(~g.mask_water),g.LAT(~g.mask_water)])) .* abs(victd_trend.(['t_std_' uv{i_uv}])(g.day_doy))';
    ResdCm.(uv{i_uv}) = f_trans_inver(ResdC.(uv{i_uv}));
    
    save(['data/speed/estsim-daily_' uv{i_uv}],'Edm','Sdm','ResdCm','-v7.3') % 5Go
end




%% ----

i_uv=1;
load(['data/speed/estsim-daily_' uv{i_uv}])

%% Figure
i_v1= find(year(g.day)==2010,1)+100;

% load('data/speed/inference-trans') % get origianl data
i_uv=1;
sfv = Edm.(uv{i_uv});
ptv = vicd.(uv{i_uv});

figure('position',[0 0 1600 900]);
tiledlayout('flow','TileSpacing','tight','Padding','tight');

for i=1:10
    nexttile; hold on; set(gca,'ydir','normal'); axis equal tight
    tmp=nan(size(g.mask_water));
    tmp(~g.mask_water) = sfv(:,i_v1+i,1);%EUm(:,i_v1+i);
    imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    borders('states','k');
    id = ~isnan(ptv(i_v1+i,:));
    scatter(radar.lon(id),radar.lat(id),100,ptv(i_v1+i,id),'filled','MarkerEdgeColor','k');
    axis([-125 -68 23 50]); 
    %caxis([-3 3]); 
    colorbar; 
    caxis([-10 10])
    title(datestr(g.day(i_v1+i)))
end

figure;
ha=tight_subplot(1,2);
axes(ha(1)); hold on; set(gca,'ydir','normal'); axis equal tight
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = EUm(:,i_v1);
imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
borders('states','k');
id = ~isnan(victd.u(i_v1,:));
scatter(radar.lon(id),radar.lat(id),100,vicd.u(i_v1,id),'filled','MarkerEdgeColor','k');
axis([-125 -68 23 50]); 
caxis([-10 10]); colorbar;
title(datestr(g.day(i_v1)))

axes(ha(2)); hold on; set(gca,'ydir','normal'); axis equal tight
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = SU(:,i_v1);
imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
borders('states','k');
axis([-125 -68 23 50]); 
% caxis([-10 10]); 
% colorbar;
title(datestr(g.day(i_v1)))
