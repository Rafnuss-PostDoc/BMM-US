%% Modeling Mird Migration as a Flow


%% Load data
load('data/density/inference-trans.mat','g','radar')
addpath('./functions/')
col2=brewermap([],'Paired');

%% Run the flow model
for i_y=2020
    simest='est';
    saveit=false;
    [Fd, Ts, gext, MVT_day] = sinksource(g,i_y,simest,saveit); % 5min
    
    simest='est';
    [Fd, Ts, gext, MVT_day] = sinksource(g,i_y,simest,saveit);
end

%% 

i_y=2021;

load(['data/flow/est_' num2str(i_y) '.mat'])
%load(['data/flow/sim_' num2str(i_y) '.mat'])

midYear = datetime(['1-July-' num2str(i_y)]);



%% 
% Map-daily
clmap = brewermap([],'Spectral');

k=1;

fig=figure('position',[0 0 1600 600]);  
ha = tight_subplot(1,2,0,[.05 .05],0);   
axes(ha(1));
h1=imagesc(g.lon,g.lat, Fd.takingoff(:,:,1),'alphaData',~isnan(Fd.takingoff(:,:,1)));
borders('states','k'); plot(radar.lon,radar.lat,'.k')
set(gca,'ydir','normal'); colorbar('south'); caxis([0 200000]);axis equal; axis([-125 -68 23 50]); colormap(gca,flipud(clmap(1:end/2,:)))
axes(ha(2));
h2=imagesc(g.lon,g.lat, -Fd.landing(:,:,1),'alphaData',~isnan(Fd.takingoff(:,:,1))); 
borders('states','k'); plot(radar.lon,radar.lat,'.k')
set(gca,'ydir','normal'); colorbar('south');caxis([0 200000]);axis equal; axis([-125 -68 23 50]); colormap(gca,clmap(end/2:end,:))
b = uicontrol('Parent',fig,'Style','slider','units', 'normalized','Position',[0.05 0.05 .9 0.05],...
              'value',1, 'min',1, 'max',numel(gext.day),'SliderStep', [1/numel(gext.day), 0.1]);
b.Callback = @(es,ed) cellfun( @(x) feval(x,es,ed), {...
    @(es,ed) set(h1,'CData', Fd.takingoff(:,:,round(es.Value+1))), ...
    @(es,ed) set(h2,'CData', -Fd.landing(:,:,round(es.Value))), ...
    @(es,ed) colorbar('south'),...
    @(es,ed) title(ha(1),['takingoff: ' datestr(gext.day(round(es.Value+1)))]),...
    @(es,ed) title(ha(2),['landing: ' datestr(gext.day(round(es.Value)))]),...
   });


%% Year-round accumulation of migratory birds on the ground
% Times serie daily

figure('position',[0 0 1000 600]); box on; grid on; hold on; 
% fill(datenum([gext.day(1); gext.day; gext.day(end)]),[0 cumsum((Ts.leaving_day+Ts.entering_day)/1000000,'omitnan') 0],'k','FaceAlpha',0.2);
h2=bar(datenum(gext.day),mean(Ts.landing_day,2)/1000000,'FaceColor',col2(6,:),'DisplayName','Landing');
h3=bar(datenum(gext.day),mean(Ts.takingoff_day,2)/1000000,'FaceColor',col2(2,:),'DisplayName','Takingoff');
h4=bar(datenum(gext.day),mean(Ts.leaving_day,2)/1000000,'FaceColor',col2(8,:),'DisplayName','Leaving');
h5=bar(datenum(gext.day),mean(Ts.entering_day,2)/1000000,'FaceColor',col2(10,:),'DisplayName','Entering +');
h6=bar(datenum(gext.day),mean(Ts.leaving_day+Ts.entering_day,2)/1000000,'k','DisplayName','\Delta takingoff-landing = \Delta Leaving- Entering');
for k=1:size(Ts.entering_day,2)
    h1=plot(datenum(gext.day),cumsum(mean(Ts.leaving_day(:,k)+Ts.entering_day(:,k),2)/1000000,'omitnan'),'-k','linewidth',2 ,'DisplayName','Cumulative \Delta takingoff/landing = # on the ground');
end
datetick('x'); axis tight; ylabel('Number of birdss (millions)')
legend([h2 h3 h4 h5 h6]);

figure('position',[0 0 1000 600]); box on; grid on; hold on; 
for k=1:size(Ts.entering_day,2)
    % h1=fill(datenum([gext.day ;flipud(gext.day)]),-[cumsum((Ts.leaving_day)/1000000,'omitnan') fliplr(cumsum(-(Ts.entering_day)/1000000,'omitnan'))],'k','FaceAlpha',0.2,'DisplayName','Difference');
    h2=plot(datenum(gext.day),cumsum(-mean(Ts.leaving_day(:,k),2)/1000000,'omitnan'),'Color',col2(8,:),'linewidth',2 ,'DisplayName','Cumulative Leaving');
    h3=plot(datenum(gext.day),cumsum(mean(Ts.entering_day(:,k),2)/1000000,'omitnan'),'Color',col2(10,:),'linewidth',2 ,'DisplayName','Cumulative Entering');
end
datetick('x'); xlabel('Date (2018)'); axis tight; legend([h2 h3]);
ylabel('Number of birdss (millions)')

%% 
% Maximum Number of bird in a night
for k=1:numel(Ts)
    [max_takeoff_sp(k),max_takeoff_sp_date(k)] = max(Ts.takingoff_day(gext.day<midYear));
    [max_takeoff_au(k),max_takeoff_au_date(k)] = max(Ts.takingoff_day);
end
disp(['Maximum number of bird taking-off is: ' num2str(max_takeoff_sp/1000000) ' happening on the ' datestr(gext.day(max_takeoff_sp_date)) '.'])
disp(['Maximum number of bird taking-off is: ' num2str(max_takeoff_au/1000000) ' happening on the ' datestr(gext.day(max_takeoff_au_date)) '.'])
%% 
% How many night to reach 50% of all derpature

tmp = Ts.takingoff_day(gext.day<midYear & ~isnan(Ts.takingoff_day'));
disp(['50% of all bird taking-off in Spring happened in ' num2str(sum( (cumsum(sort(tmp,'descend'),'omitnan') ./ nansum(tmp)) < .5) )])
tmp = Ts.takingoff_day(gext.day>midYear & ~isnan(Ts.takingoff_day'));
disp(['50% of all bird taking-off in Autumn happened in ' num2str(sum( (cumsum(sort(tmp,'descend'),'omitnan') ./ nansum(tmp)) < .5) )])
%% 
% Accumulation of bird in Spring

spring_entering = nansum(Ts.entering_day(gext.day<midYear)/1000000)
spring_leaving = nansum(Ts.leaving_day(gext.day<midYear)/1000000)
spring_delta = spring_entering + spring_leaving
%% 
% Accumulation (or dissipation) of bird in Autumns

autumn_entering = nansum(Ts.entering_day(gext.day>midYear)/1000000)
autumn_leaving = nansum(Ts.leaving_day(gext.day>midYear)/1000000)
autumn_delta = autumn_entering + autumn_leaving
%% 
% Recrutement = Flux out/Flux in = Bird in + reproduction / Bird in

disp(['Ratio of Autumn departure over spring arrival: ' num2str(mean(-autumn_leaving ./ spring_entering)) ])
disp(['Ratio of Autumn accumulation over spring accumulation: ' num2str(mean(-autumn_delta ./ spring_delta))])
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%% Seasonal Flow: In and Out
% *Flux per sector (15min)*

k=1;

tmp=cat(3,...
    nansum(Fd.entering(:,:,gext.day<midYear,k),3),...
    nansum(Fd.leaving(:,:,gext.day<midYear,k),3),...
    nansum(Fd.entering(:,:,gext.day>midYear,k),3),...
    nansum(Fd.leaving(:,:,gext.day>midYear,k),3)...
    );

lg={'entering spring','leaving spring','entering autumn','leaving autumn'};
figure; 
ha=tight_subplot(2,2,[.05 .03]); colormap(brewermap([],'Spectral'));
for i=1:4
    axes(ha(i)); hold on;
    imagesc(tmp(:,:,i),'AlphaData',nansum(Fd.entering(:,:,:,k),3)>0);
    set(gca,'ydir','normal'); caxis([-10^7 10^7]); set(gca,'Color','k');
    c=colorbar('south'); c.Color="w"; title(lg{i})
end

% export as png
%imwrite(mean(Fout.Fout>0,3), 'data/flow/sectors.png');
sect=imread('data/flow/sectors_classified.tif');
transect_name={'east', 'northeast', 'north','west', 'mexico','gulf'};

figure('position',[0 0 1000 400]);
imagesc(g.lon,g.lat,sect)
c=colorbar; c.Ticks=1:6; c.TickLabels=transect_name;

Fd.s=nan(gext.nat,numel(transect_name),2,size(Fd.entering,4));
for i_s=1:numel(transect_name)
    Fd.s(:,i_s,1,:) =  sum(reshape(Fd.entering(repmat(sect==i_s,1,1,gext.nat,size(Fd.entering,4))),[],gext.nat,size(Fd.entering,4)));
    Fd.s(:,i_s,2,:) =  sum(reshape(Fd.leaving(repmat(sect==i_s,1,1,gext.nat,size(Fd.leaving,4))),[],gext.nat,size(Fd.leaving,4)));
end

% Flux per sector (season)
Fd.s_season=splitapply(@(x) sum(x,'omitnan'),Fd.s,(gext.day<midYear)+1);
% season | sector | entering/leaving | simulation

% Flux per sector (daily)

figure('position',[0 0 1000 1000]); 
for i_s=1:numel(transect_name)
    subplot(numel(transect_name),1,i_s); hold on;
    bar(gext.day, mean(Fd.s(:,i_s,1),4)/1000000,'FaceColor',col2(8,:));
    bar(gext.day, -mean(Fd.s(:,i_s,2),4)/1000000,'FaceColor',col2(10,:));
    ylabel(transect_name{i_s})
    ylim([0 max(Fd.s(:))]/1000000)
    box on; grid on
end
%% 
% Figure

figure('position',[0 0 1000 400]); box on
bar(abs(Fd.s_season)'/1000000,'stacked'); axis tight
xticklabels({'Spring -', 'Spring +','Autumn -','Autumn +'})
ylabel('Number of birds crossing the transect (-) outward, (+) inward [millions]')
legend(transect_name)
tmp = [sum(Fd.s_season(:,1:2),2) sum(Fd.s_season(:,3:4),2)]/1000000;
%% 
% Spring and Autum (col) fluxes over the transects

var_name = {'Transect','Spring','Autumn','Ratio'};
table(transect_name',tmp(:,1),tmp(:,2),tmp(:,2)./tmp(:,1),'VariableNames',var_name)
%% 
% Internal change

[
    table({'Northern'},sum(tmp([1 2 6],1)),sum(tmp([1 2 6],2)),sum(tmp([1 2 6],2))./sum(tmp([1 2 6],1)),'VariableNames',var_name);
    table({'Southern'},sum(tmp(3:5,1)),sum(tmp(3:5,2)),sum(tmp(3:5,2))./sum(tmp(3:5,1)),'VariableNames',var_name);
    table({'Global'},sum(tmp(:,1)),sum(tmp(:,2)),sum(tmp(:,2))./sum(tmp(:,1)),'VariableNames',var_name);
    ]
%% 
% Sum of grouped transect

figure; subplot(1,2,1)
pie(flipud(abs(tmp(:,1))),fliplr(transect_name))
subplot(1,2,2)
pie(flipud(abs(tmp(:,2))),fliplr(transect_name))
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% 
%% *Seasonal Flow: Within the study area*
% Arrival/Departure over season

figure('position',[0 0 1000 600]); 

id = gext.day<midYear;
subplot(2,3,2); hold on
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(Fd.landing2(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON, W_tmp./g.area ); 
plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Spring landing [bird/km^2]'; caxis([-4000 4000])
subplot(2,3,1); hold on
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(Fd.takingoff(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON,W_tmp./g.area ); plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Spring takingoff [bird/km^2]'; caxis([-4000 4000])
subplot(2,3,3); hold on
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(W.day.W(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON,W_tmp./g.area ); plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Spring diff [bird/km^2]'; caxis([-2000 2000])

id = gext.day>midYear;
subplot(2,3,5); hold on
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(Fd.landing2(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON,W_tmp./g.area ); plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Autumn landing [bird/km^2]'; caxis([-4000 4000])
subplot(2,3,4); hold on; 
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow( 'landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(Fd.takingoff(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON,W_tmp./g.area ); plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Autumn takingoff [bird/km^2]'; caxis([-4000 4000])
subplot(2,3,6); hold on; 
h = worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow( 'landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
W_tmp=nansum(W.day.W(:,:,id),3);
W_tmp(W_tmp==0)=NaN;
surfm(g.LAT,g.LON,W_tmp./g.area ); plotm(data.lat,data.lon,'.k')
c=colorbar;c.Label.String = 'Autumn diff [bird/km^2]'; caxis([-2000 2000])
colormap(brewermap([],'Spectral'))
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% Appendix
%% Illutration main text
% Illustration of the fluxes (Figure 1)

files={'gt30w020n90'};%,'gt30e020n90'}; % ,'gt30e020n40','gt30w020n40'
[DEM,GeoCellRef] = geotiffread(['data/' files{1}]);
GeoCellRef_lon = linspace(GeoCellRef.LongitudeLimits(1),GeoCellRef.LongitudeLimits(2),GeoCellRef.RasterSize(2));
GeoCellRef_lat = linspace(GeoCellRef.LatitudeLimits(2),GeoCellRef.LatitudeLimits(1),GeoCellRef.RasterSize(1));
[GeoCellRef_LON,GeoCellRef_LAT] = meshgrid(GeoCellRef_lon,GeoCellRef_lat);
id_lat = GeoCellRef_lat>g.lat(1)-3 & GeoCellRef_lat<g.lat(end)+3;
id_lon = GeoCellRef_lon>g.lon(1)-3 & GeoCellRef_lon<g.lon(end)+3;
DEM(DEM>-5000&DEM<=0)=0;

figure('position',[0 0 800 600]);hold on
h = worldmap([g.lat(1)-3 g.lat(end)+3], [g.lon(1)-3 g.lon(end)+3]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
%geoshow(GeoCellRef_LAT(id_lat,id_lon), GeoCellRef_LON(id_lat,id_lon), DEM(id_lat,id_lon),'DisplayType','texturemap'); demcmap(DEM(id_lat,id_lon))
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255);
plotm([data.lat],[data.lon],'.k')
figure('position',[0 0 800 600]);hold on
h = worldmap([g.lat(1)-3 g.lat(end)+3], [g.lon(1)-3 g.lon(end)+3]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
tmp = sum(rho,3,'omitnan'); tmp(~g.mask_water)=nan;
surfm(g.LAT,g.LON, tmp)
%colormap(flipud(autumn))
figure('position',[0 0 800 600]);hold on
h = worldmap([g.lat(1)-3 g.lat(end)+3], [g.lon(1)-3 g.lon(end)+3]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
tmp1 = nansum(vy,3); tmp1(~g.mask_water)=nan;
tmp2 = nansum(vy,3); tmp2(~g.mask_water)=nan;
quiverm(g.LAT(1:3:end,1:3:end),g.LON(1:3:end,1:3:end), tmp1(1:3:end,1:3:end), tmp2(1:3:end,1:3:end),'k')
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% Supplementary Material not published
% Map of the center of gravity

tmp = -reshape(Fd.landing(repmat(g.mask_water,1,1,gext.nat)), g.nlm, gext.nat);
lat_landing = sum(tmp .* repmat(g.LAT(g.mask_water),1,gext.nat)) ./ sum(tmp);
lon_landing = sum(tmp .* repmat(g.LON(g.mask_water),1,gext.nat)) ./ sum(tmp);
tmp = -reshape(Fd.takingoff(repmat(g.mask_water,1,1,gext.nat)), g.nlm, gext.nat);
lat_takingoff = sum(tmp .* repmat(g.LAT(g.mask_water),1,gext.nat)) ./ sum(tmp);
lon_takingoff = sum(tmp .* repmat(g.LON(g.mask_water),1,gext.nat)) ./ sum(tmp);

figure('position',[0 0 1000 600]); 
h = worldmap([46 52], [1 11]);
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
for i_t=1:gext.nat-1
    h1=plotm([lat_takingoff(i_t)', lat_landing(i_t)], [lon_takingoff(i_t) lon_landing(i_t)'], '-k');
end
h2=scatterm(lat_landing', lon_landing',[], datenum(gext.day), 'v','filled');
h3=scatterm(lat_takingoff', lon_takingoff',[], datenum(gext.day),'^','filled');
legend([h1 h2 h3],{'displacement of mass','depature','landing'});
c=colorbar; c.TickLabels=datestr(c.Ticks);

r = lldistkm([lat_landing', lon_landing'], [lat_takingoff', lon_takingoff']);
ThetaInDegrees = atan2d( (lat_landing-lat_takingoff) , (lon_landing-lon_takingoff) );
figure;
for i_b=1:numel(data.block.date)
    polarscatter(deg2rad(ThetaInDegrees(gext.day_b==i_b)),r(gext.day_b==i_b),[],data.block.col(i_b,:),'filled'); hold on;
end
axis tight; legend(data.block.date_str,'Location','north')
%% 
% load sunrise data


%% 
% Passage vs stop



id = ~isnan(W.day.W) & ~isnan(MVT_day) & MVT_day~=0; 
figure; hold on;
scatter(abs(W.day.W(id)),MVT_day(id),'.k');
% ksdensity([abs(W.day.W(id)),MVT_day(id)])
p=polyfit(abs(W.day.W(id)),MVT_day(id),1);
h2=plot(1:100:max(abs(W.day.W(id))), polyval(p,1:100:max(abs(W.day.W(id)))),'LineWidth',2);
xlabel('Daily difference of landing and takingoff (absolute value). [bird/day]')
ylabel('Average MTR of the night [bird/km/hr]')
legend(h2,['correlation of ' num2str(corr(abs(W.day.W(id)),MVT_day(id)))])
nancov(Ts.entering_day',Ts.leaving_day')
nancov(Ts.entering_day(2:end)',Ts.leaving_day(1:end-1)')
nancov(Ts.entering_day(1:end-1)',Ts.leaving_day(2:end)')
%% 
% Timeseries-total space high resolution

figure('position',[0 0 1000 600]);
hold on;legend
ylabel('Flux from/to the ground: landing(-) and takingoff (+)')
plot(gext.time(2:end),Ts.W,'--k','DisplayName','\Delta takingoff/landing');
plot(gext.time(2:end),Ts.landing,'Color',col2(1,:),'DisplayName','landing');
plot(gext.time(2:end),Ts.takingoff,'Color',col2(end,:),'DisplayName','takingoff');
yyaxis('right')
plot(gext.time(2:end),cumsum(Ts.W),'-k','DisplayName','Cumulative Flux delta landing/takingoff');
axis tight;
%% 
%% 
% Total Number of birds in the air along the year

tmp = hour(gext.time)*60+minute(gext.time);
tmp(tmp>720) = tmp(tmp>720)-24*60;
id_tmp = tmp/15+12*60/15;
nb=nan(24*60/15,gext.nat);
for i_t=1:gext.nat
    idt=gext.day_id(1:end-1)==i_t;
    nb(id_tmp(idt),i_t)=reshape(sum(sum(rho0(:,:,idt),1),2),[],1);
end

% knownday = ismember(datenum(g.day), unique(data.day(data.day_id)));
% nb(:,~knownday)=nan;

figure('position',[0 0 1000 600]); hold on;
imagesc(datenum(g.day),(-12:12)*60,nb/1000000);%,'AlphaData',~isnan(nb)); 
set(gca,'ydir','normal')
datetick('x'); axis tight; ylim([-480 480]);
yticks(-480:120:480); yticklabels(-8:2:8);
xlabel('Date'); ylabel('Hours before and after midnight');
c=colorbar;c.Label.String='Number of birds in the air [millions]';

plot(datenum(g.day),hour(mean(gg_dusk,1))*60+minute(mean(gg_dusk,1))-24*60,'-k')
plot(datenum(g.day),hour(mean(gg_dawn,1))*60+minute(mean(gg_dawn,1)),'-k')
plot(datenum(g.day),mean(hour(gg_sunset)*60+minute(gg_sunset))-24*60,'--k')
plot(datenum(g.day),mean(hour(gg_sunrise)*60+minute(gg_sunrise)),'--k')
%% 
% takingoff/lading during the nights
% 
% Over the next 15minutes. 

gg_dawn = datenum(g.dawn(:,Locb(gext.day_id)));
gg_dusk = datenum(g_dusk(:,Locb(gext.day_id)-1));

figure('position',[0 0 1000 600]);
tmp = reshape(W.takingoff(repmat(g.mask_water,1,1,gext.nt-1)),g.nlm,gext.nt-1);
tmp2 = (repmat(datenum(gext.time(1:end-1)),1,g.nlm)'-datenum(gg_dusk(:,1:end-1)) ) * 24;
% figure; scatter(tmp2(:),tmp(:),'.k')
[G,ID] = findgroups(round(tmp2(:),2));
tmp_grp = splitapply(@nansum,tmp(:),G);
cumtmp_grp=cumsum(tmp_grp);
subplot(2,1,1); hold on
fill([0 15/60 15/60 0],[0 0 6 6],'r')
plot(ID+15/60,cumtmp_grp/1000000000,'k','LineWidth',2); 
xlim([0 10]); xlabel('Hours of dusk'); ylabel('Number of birds takeoff in the next 15min [billion of bird]')

tmp = reshape(-W.landing(repmat(g.mask_water,1,1,gext.nt-1)),g.nlm,gext.nt-1);
tmp2 = (repmat(datenum(gext.time(1:end-1)),1,g.nlm)'-datenum(gg_dawn(:,1:end-1)) ) * 24;
% figure; scatter(tmp2(:),tmp(:),'.k')
[G,ID] = findgroups(round(tmp2(:),2));
tmp_grp = splitapply(@nansum,tmp(:),G);
cumtmp_grp=cumsum(tmp_grp);
cumtmp_grp(ID+15/60==0)=cumtmp_grp(end);
subplot(2,1,2); hold on;
fill([-15/60 0 0 -15/60],[0 0 6 6],'r')
plot(ID+15/60,cumtmp_grp/1000000000,'k','LineWidth',2); xlim([-10 0]); ylim([0 6])
xlabel('Hours before dawn'); ylabel(['Number of birds landing in the next 15min [billions of bird]'])

