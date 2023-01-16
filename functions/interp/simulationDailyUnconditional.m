function Res=simulationDailyUnconditional(duv, g, gn, sim_n, dcov, plotit)

% gn : group of days simulated together (15 or 30)
% sim_n: number of simulaiton
% plotit: display plot or not

it=1:gn:numel(g.day);

dlat = lldistkm([g.lat((end-1)/2) mean(g.lon)],[g.lat((end-1)/2+1) mean(g.lon)]);
dlon = lldistkm([mean(g.lat) g.lon((end-1)/2) ],[ mean(g.lat) g.lon((end-1)/2+1)]);
[Ddistlat, Ddistlon] =  ndgrid(-numel(g.lat):numel(g.lat), -numel(g.lon):numel(g.lon));
Ddist = sqrt((dlat*Ddistlat).^2+(dlon*Ddistlon).^2);

Dtime = abs(-gn:gn);

g_nlm = sum(~g.mask_water(:));
res = nan(g_nlm,gn*2+1,numel(it),sim_n,'single');
g_nlat = numel(g.lat); 
g_nlon=numel(g.lon);
g_mw = repmat(~g.mask_water,1,1,gn*2+1);

%% Perform simulation 
tic
for iit=1:numel(it)

    if duv=='d'
        parm = dcov.f_parm(g.day_doy(it(iit)));
        Cdist = dcov.Gneiting_dist(Ddist,parm(4),parm(7));
        Ctime = dcov.Gneiting_time(Dtime,parm(5),parm(6));
        K = parm(3) .* Cdist.*reshape(Ctime,1,1,[]);
        % K(end/2+.5,end/2+.5,:) = K(end/2+.5,end/2+.5,:)+parm(2);
        % K(end/2+.5,end/2+.5,end/2+.5) = K(end/2+.5,end/2+.5,end/2+.5)+parm(1);
    else
        if duv=='u'
            parm = dcov.f_parm_u(g.day_doy(it(iit)));
        else
            parm = dcov.f_parm_v(g.day_doy(it(iit)));
        end
        Cdist = dcov.Gneiting_dist(Ddist,parm(3),parm(6));
        Ctime = dcov.Gneiting_time(Dtime,parm(4),parm(5));
        K = parm(2) .* Cdist.*reshape(Ctime,1,1,[]);
        % K(end/2+.5,end/2+.5,end/2+.5) = K(end/2+.5,end/2+.5,end/2+.5)+parm(1);
    end

    fftnK = fftn(K).^.5;
    
    for k=1:sim_n
        tmp = real(ifftn(fftn(randn(size(K))).*fftnK));
        tmp = tmp(1:g_nlat,1:g_nlon,:);
        res(:,:,iit,k) = reshape(tmp(g_mw),[],gn*2+1);
    end
    ttoc = toc;
    disp([ num2str(iit) '/' num2str(numel(it)) ' (' num2str(round(iit/numel(it)*100)) '%) in ' num2str(ttoc/60) ' min']) 
end

%% Figure illustration 
if plotit
    figure; ha=tight_subplot(5,5); u=10;
    for iit=1:25
        axes(ha(iit)); set(gca,'ydir','normal'); hold on;
        tmp=nan(g_nlat,g_nlon);
        tmp(~g.mask_water) = res(:,u,iit,1);
        imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water)
        title(datestr(g.day(it(iit)+u)))
        borders('states','k');
        axis equal;axis([-125 -68 23 50]);
    end
    
    figure; ha=tight_subplot(5,5);
    iit=8;
    for u=1:25
        axes(ha(u)); set(gca,'ydir','normal'); hold on;
        tmp=nan(g_nlat,g_nlon);
        tmp(~g.mask_water) = res(:,u,iit,1);
        imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water)
        title(datestr(g.day(it(iit)+u)))
        borders('states','k');
        axis equal;axis([-125 -68 23 50]);
    end
end

%% Merge
Res = nan(g_nlm,gn*numel(it)+1,sim_n,'single');
t=linspace(0,pi/2,gn+2);
t=t(2:end-1);
tic
for k=1:sim_n
    for iit=1:numel(it)-1
        Res(:,it(iit),k) = res(:,end/2-.5,iit,k);
        Res(:,it(iit)+(1:gn),k) = res(:,end-gn+1:end,iit,k).*cos(t) + res(:,1:gn,iit+1,k).*sin(t);
    end
    ttoc = toc;
    disp([ num2str(k) '/' num2str(sim_n) ' (' num2str(round(k/sim_n*100)) '%) in ' num2str(ttoc/60) ' min']) 
end
Res = Res(:,1:numel(g.day),:);

%% 

if plotit
    %% Figure Illustration 
    i_v1= find(year(g.day)==2000,1);
    
    figure; hold on; set(gca,'ydir','normal'); axis equal tight
    tmp=nan(g_nlat,g_nlon);
    tmp(~g.mask_water) = Res(:,i_v1,k);
    im=imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    borders('states','k');
    axis([-125 -68 23 50]); colorbar; caxis([-3 3])
    
    for i_v=i_v1:numel(g.day)
        tmp=nan(g_nlat,g_nlon);
        tmp(~g.mask_water) = Res(:,i_v,k);
        im.CData = tmp;
        title(datestr(g.day(i_v)))
        refreshdata
        drawnow
        pause(.2);
    end
    
    figure; ha=tight_subplot(5,5);
    iit=3;
    for u=1:25
        axes(ha(u)); set(gca,'ydir','normal'); hold on;
        tmp=nan(g_nlat,g_nlon);
        tmp(~g.mask_water) = Res(:,it(iit)+u*3,k);
        imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water)
       title(datestr(g.day(it(iit)+u*3)))
        borders('states','k');
        axis equal;axis([-125 -68 23 50]);
    end
    
    end

end