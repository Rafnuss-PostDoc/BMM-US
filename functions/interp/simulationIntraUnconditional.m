function Res = simulationIntraUnconditional(duv, g, radar, idt_y, sim_n, icov, plotit)

% gn : group of days simulated together (15 or 30)
% sim_n: number of simulaiton
% plotit: display plot or not

dlat = lldistkm([g.lat((end-1)/2) mean(g.lon)],[g.lat((end-1)/2+1) mean(g.lon)]);
dlon = lldistkm([mean(g.lat) g.lon((end-1)/2) ],[ mean(g.lat) g.lon((end-1)/2+1)]);
[Ddistlat, Ddistlon] =  ndgrid(-numel(g.lat):numel(g.lat), -numel(g.lon):numel(g.lon));
Ddist = sqrt((dlat*Ddistlat).^2+(dlon*Ddistlon).^2);

g_nlm = sum(~g.mask_water(:));
g_nlat = numel(g.lat); 
g_nlon=numel(g.lon);

%% Perform simulation 
tic

% initial variable
Res = nan(g_nlm, numel(idt_y), sim_n, 'single');

 % find the day of year
day_id = g.day_id(idt_y);
day_id_unique = unique(day_id);

for i_d=1:5%numel(day_id_unique)

    % Find the index of the element to sim/est for the year (i.e., in idt_y)
    i_ey = find(day_id==day_id_unique(i_d));

    % Find the index of the element to estimate on the full time (absolute id)
    i_en = idt_y(i_ey);
    
    % mask water
    g_mw = repmat(~g.mask_water,1,1,numel(i_ey));

    % compute their delta time
    dt_max=days(g.time(i_en(end))-g.time(i_en(1)));
    % dt_max=4/24; % we can assume a max range of 4 hours
    dt = -dt_max:(1/24/4):dt_max; 
  
    % Define the parameter of the covariance
    if duv=='d'
        parm = icov.f_parm(g.day_doy(day_id_unique(i_d)));
    elseif duv=='u'
        parm = icov.f_parm_u(g.day_doy(day_id_unique(i_d)));
    else
        parm = icov.f_parm_v(g.day_doy(day_id_unique(i_d)));
    end
    Cdist = icov.Gneiting_dist(Ddist,parm(3),parm(6));
    Ctime = icov.Gneiting_time(abs(dt),parm(4),parm(5));
    K = parm(2) .* Cdist .* reshape(Ctime,1,1,[]);
    
    fftnK = fftn(K).^.5;

    for k=1:sim_n
        tmp = real(ifftn(fftn(randn(size(K))).*fftnK));
        tmp = tmp(1:g_nlat,1:g_nlon,1:numel(i_en));
        Res(:,i_ey,k) = reshape(tmp(g_mw),[],numel(i_ey));
    end
    ttoc = toc;
    disp(['estuncond - ' datestr(g.day(day_id_unique(i_d))) ' - ' num2str(ttoc/60,3) ' min'])  
end


%% Figure Illustration 
if plotit
    k=1;
    i_s = 50;
    i_step = 1; 

    figure; ha=tight_subplot(5,5);
    for u=1:25
        axes(ha(u)); set(gca,'ydir','normal'); hold on;
        tmp=nan(g_nlat,g_nlon);
        tmp(~g.mask_water) = Res(:,i_s+u*i_step,k);
        imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water)
        id = ~isnan(vidti(i_s+u*i_step,:));
        scatter(radar.lon(id),radar.lat(id),[], vidti(i_s+u*i_step,id),'filled','MarkerEdgeColor','k')
        title(datestr(g.time(idt_y(i_s)+u*i_step)))
        borders('states','k');
        axis equal;axis([-125 -68 23 50]);
    end
end

end