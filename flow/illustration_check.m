load('data/density/inference-trans.mat')
load('data/speed/inference-trans.mat')
addpath('./functions/')

% define year
i_y=2021;



% load data
load(['data/density/est_' num2str(i_y) '.mat'],'idt_y')
g_mw = repmat(~g.mask_water,1,1,numel(idt_y));
vy = nan(numel(g.lat), numel(g.lon), numel(idt_y)); vx=vy; rho = vy;
if false
    load(['data/density/sim_' num2str(i_y) '_1.mat'])
    load(['data/speed/sim_u_' num2str(i_y) '_1.mat'])
    load(['data/speed/sim_v_' num2str(i_y) '_1.mat'])
    
    rho(g_mw) = ResCmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
    vx(g_mw) = ResCm.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
    vy(g_mw) = ResCm.v/1000*60*60; % m/s -> km/h (+) north, (-) south
    clear ResCmbT ResCm
else
    load(['data/density/est_' num2str(i_y) '.mat'])
    rho(g_mw) = EmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
    load(['data/speed/est_uv_' num2str(i_y) '.mat'])
    vx(g_mw) = Em.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
    vy(g_mw) = Em.v/1000*60*60; % m/s -> km/h (+) north, (-) south
end

gy.time = g.time(idt_y);
gy.day_id = g.day_id(idt_y)-min(g.day_id(idt_y))+1;
gy.day = g.day(unique(g.day_id(idt_y)));

% known data
vid = trans.f_inv(vidTS);

vic.u = (vici.u + vicd.u(g.day_id,:))/1000*60*60;
vic.v = (vici.v + vicd.v(g.day_id,:))/1000*60*60;

vic.u = vic.u(idt_y,:); 
vic.v = vic.v(idt_y,:);
vid = vid(idt_y,:);







%% Check histogram

figure; hold on;
histogram(log(vid(vid>0)),'Normalization','pdf','EdgeAlpha',0)
histogram(log(rho(rho>0)),'Normalization','pdf','EdgeAlpha',0)
xlim([trans.thr_T 10])
legend('Hard data','sim/est')
% set(gca,'yscale','log')

figure; hold on;
histogram(vic.u*1000/60/60,'Normalization','pdf','EdgeAlpha',0)
histogram(vx*1000/60/60,'Normalization','pdf','EdgeAlpha',0)
legend('Hard data','sim/est'); axis tight;

figure; hold on;
histogram(vic.v*1000/60/60,'Normalization','pdf','EdgeAlpha',0)
histogram(vx*1000/60/60,'Normalization','pdf','EdgeAlpha',0)
legend('Hard data','sim/est'); axis tight;

%% Check value at hard data

[~,id_cond] = min( sqrt( (radar.lon'-g.LON(~g.mask_water)).^2 + (radar.lat'-g.LAT(~g.mask_water)).^2 ) );
tmp = reshape(rho(g_mw),sum(~g.mask_water(:)),[]);
rho_hd = tmp(id_cond,:)';
tmp = reshape(vx(g_mw),sum(~g.mask_water(:)),[]);
vx_hd = tmp(id_cond,:)';
tmp = reshape(vy(g_mw),sum(~g.mask_water(:)),[]);
vy_hd = tmp(id_cond,:)';

figure;
subplot(3,1,1); imagesc((rho_hd)'); caxis([-4 4])
subplot(3,1,2); imagesc((vid)'); caxis([-4 4])
subplot(3,1,3); imagesc((vid-rho_hd)')

figure; histogram(vid-rho_hd)
figure; histogram(nanmean(vid-rho_hd))

figure; imagesc(vic.v-vy_hd)
figure; histogram(vic.v-vy_hd)
figure; plot(vic.v-vy_hd,vic.u,'.k')


% vid0 = readmatrix('data/vid.csv','NumHeaderLines',1,'TreatAsMissing' ,'NA');
[u,id]=sort(abs(vid(:)-rho_hd(:)),'descend');
id = id(~isnan(u));
figure; hold on;
for i=8:100
    [id_t,id_r]=ind2sub(size(vid),id(i));
    id_tt = (-4*24:4*24) + id_t;
    id_tt=id_tt(id_tt>0);
    clf; hold on; box on; grid on;
    ylabel(radar.name(id_r))
    plot(gy.time(id_tt),vid0(idt_y(id_tt),id_r))
    plot(gy.time(id_tt),vid(id_tt,id_r))
    plot(gy.time(id_tt),rho_hd(id_tt,id_r))
    keyboard
end

figure('position',[0 0 1600 900]);borders('states','color',[.2 .2 .2]*3);axis([-125 -68 23 50]); plot(radar.lon(id_r),radar.lat(id_r),'ok');
    for i_r=1:height(radar)
        text(radar.lon(i_r),radar.lat(i_r),[radar.name{i_r} ' (' num2str(i_r) ')'],'VerticalAlignment','top','HorizontalAlignment','center')
    end

%% Check MTR

mtr_hd = vid .* sqrt(vic.u.^2 + vic.v.^2) /4;

mtr_interp = rho .* sqrt(vx.^2 + vy.^2) /4;

midYear = datetime(['1-Jul-' num2str(i_y)]);

figure('position',[0 0 1200 600]);  
tiledlayout('flow','TileSpacing','tight','Padding','tight'); 
nexttile; hold on;
imagesc(g.lon,g.lat, sum(mtr_interp(:,:,midYear>gy.time),3,'omitnan'),'alphaData',~all(isnan(mtr_interp),3));
scatter(radar.lon,radar.lat,100,sum(mtr_hd(midYear>gy.time,:),'omitnan'),'filled','MarkerEdgeColor','k');
borders('states','k'); 
set(gca,'ydir','normal'); colorbar('south'); axis equal; axis([-125 -68 23 50]); 
nexttile; hold on;
imagesc(g.lon,g.lat, sum(mtr_interp(:,:,midYear<gy.time),3,'omitnan'),'alphaData',~all(isnan(mtr_interp),3));
scatter(radar.lon,radar.lat,100,sum(mtr_hd(midYear<gy.time,:),'omitnan'),'filled','MarkerEdgeColor','k');
borders('states','k'); 
set(gca,'ydir','normal'); colorbar('south'); axis equal; axis([-125 -68 23 50]); 



%%
% figure; plot(gy.time,nanmean(vid,2))

t_tmp= datetime('2021-9-30');
id_day = find(gy.day_id==find(gy.day==t_tmp));

% figure; plot(gy.time(id_day),vid(id_day,:))


fig=figure('position',[0 0 1200 600]);  
ha = tight_subplot(1,1,0,0,0);  hold on;
axes(ha(1));
X = (rho(:,:,id_day));
duv=10; duvs=10;
U = vx(1:duv:end,1:duv:end,id_day)/duvs;
V = vy(1:duv:end,1:duv:end,id_day)/duvs;
dr = (vid(id_day,:))';
ur = vic.u(id_day,:)'/duvs;
vr = vic.v(id_day,:)'/duvs;
h1 = imagesc(g.lon,g.lat, X(:,:,1),'alphaData',~all(isnan(X),3));
plot(radar.lon,radar.lat,'.k')
h4 = scatter(radar.lon,radar.lat,100,dr(:,1),'filled','MarkerEdgeColor','k');
h2 = quiver(g.LON(1:duv:end,1:duv:end),g.LAT(1:duv:end,1:duv:end),U(:,:,1),V(:,:,1),0,'k');
h3 = quiver(radar.lon,radar.lat,ur(:,1),vr(:,1),0,'r');
borders('states','k'); 
set(gca,'ydir','normal'); colorbar('south'); axis equal; axis([-125 -68 23 50]); 
% colormap(gca,flipud(clmap(1:end/2,:)))
caxis([10 300]);

b = uicontrol('Parent',fig,'Style','slider','units', 'normalized','Position',[0.05 0.05 .9 0.05],...
'value',1, 'min',1, 'max',size(X,3),'SliderStep', [1/size(X,3), 0.1]);

b.Callback = @(es,ed) cellfun( @(x) feval(x,es,ed), {...
    @(es,ed) set(h1,'CData', X(:,:,round(es.Value+1))), ...
    @(es,ed) set(h2,'UData', U(:,:,round(es.Value+1))), ...
    @(es,ed) set(h2,'VData', V(:,:,round(es.Value+1))), ...
    @(es,ed) set(h3,'UData', ur(:,round(es.Value+1))), ...
    @(es,ed) set(h3,'VData', vr(:,round(es.Value+1))), ...
    @(es,ed) set(h4,'CData', dr(:,round(es.Value+1))), ...
    @(es,ed) colorbar('south'),...
    @(es,ed) title(ha(1),datestr(gy.time(id_day(round(es.Value+1))))),...
   });





%% Check timeserie
i_r=40;
[~,id_lat] = min(abs(radar.lat(i_r)-g.lat));
[~,id_lon] = min(abs(radar.lon(i_r)-g.lon));

figure; hold on
plot(gy.time,vid(:,i_r))
plot(gy.time,squeeze(rho(id_lat,id_lon,:)))
legend("data", "interpolation")

%% Daily scale

rhod=nan(numel(g.lat), numel(g.lon), max(gy.day_id)); 
vxd=rhod;
vyd=rhod;
vidd=nan(max(gy.day_id),146);
vic.ud=vidd;
vic.vd=vidd;
for i = 1:max(gy.day_id)
    id = gy.day_id==i;

    rhod(:,:,i) = mean(rho(:,:,id),3,'omitnan');
%     vxd(:,:,i) = sum(vx(:,:,id) .* rho(:,:,id),3,'omitnan') ./ sum(rho(:,:,id),3,'omitnan');
%     vyd(:,:,i) = sum(vy(:,:,id) .* rho(:,:,id),3,'omitnan') ./ sum(rho(:,:,id),3,'omitnan');
    vxd(:,:,i) = mean(vx(:,:,id),3,'omitnan');
    vyd(:,:,i) = mean(vy(:,:,id),3,'omitnan');

    vidd(i,:) = mean(vid(id,:),'omitnan');
    %w = vid(id,:); w(isnan(vic.u(id)))=nan;
    vic.ud(i,:) = mean(vic.u(id,:),'omitnan');
    % w = vid(id,:); w(isnan(vic.v(id)))=nan;
    vic.vd(i,:) = mean(vic.v(id,:),'omitnan');

end


%% 

fig=figure('position',[0 0 1200 600]);  
ha = tight_subplot(1,1,0,0,0);  hold on;
axes(ha(1));
X = rhod;
duv=10;duvs=10;
U = vxd(1:duv:end,1:duv:end,:)/duvs;
V = vyd(1:duv:end,1:duv:end,:)/duvs;
dr = vidd;
ur = vic.ud'/duvs;
vr = vic.vd'/duvs;
h1 = imagesc(g.lon,g.lat, X(:,:,1),'alphaData',~isnan(X(:,:,1)));
h4 = scatter(radar.lon,radar.lat,100,dr(1,:)','filled','MarkerEdgeColor','k');
h2 = quiver(g.LON(1:duv:end,1:duv:end),g.LAT(1:duv:end,1:duv:end),U(:,:,1),V(:,:,1),0,'k');
h3 = quiver(radar.lon,radar.lat,ur(:,1),vr(:,1),0,'r');
borders('states','k'); % plot(radar.lon,radar.lat,'.k')
set(gca,'ydir','normal'); colorbar('south'); axis equal; axis([-125 -68 23 50]); 
% colormap(gca,flipud(clmap(1:end/2,:)))
caxis([0 250]);

b = uicontrol('Parent',fig,'Style','slider','units', 'normalized','Position',[0.05 0.05 .9 0.05],...
'value',1, 'min',1, 'max',size(X,3),'SliderStep', [1/size(X,3), 0.1]);

b.Callback = @(es,ed) cellfun( @(x) feval(x,es,ed), {...
    @(es,ed) set(h1,'CData', X(:,:,round(es.Value+1))), ...
    @(es,ed) set(h2,'UData', U(:,:,round(es.Value+1))), ...
    @(es,ed) set(h2,'VData', V(:,:,round(es.Value+1))), ...
    @(es,ed) set(h3,'UData', ur(:,round(es.Value+1))), ...
    @(es,ed) set(h3,'VData', vr(:,round(es.Value+1))), ...
    @(es,ed) set(h4,'CData', dr(round(es.Value+1),:)), ...
    @(es,ed) colorbar('south'),...
    @(es,ed) title(ha(1),datestr(gy.day(round(es.Value+1)))),...
   });


