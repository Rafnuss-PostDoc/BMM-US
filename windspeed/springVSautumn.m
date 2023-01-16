cd('/Users/raphael/Library/CloudStorage/Box-Box/BMM-US/')
set(0, 'DefaultAxesBox', 'on');
addpath(genpath('functions'))

load('data/density/inference-trans.mat')
vid = trans.f_inv(vidTS);

viv = readmatrix('data/viv.csv','NumHeaderLines',1,'TreatAsMissing' ,'NA');
viu = readmatrix('data/viu.csv','NumHeaderLines',1,'TreatAsMissing' ,'NA');
viwv = readmatrix('data/viwv.csv','NumHeaderLines',1,'TreatAsMissing' ,'NA');
viwu = readmatrix('data/viwu.csv','NumHeaderLines',1,'TreatAsMissing' ,'NA');

% transform into complex variable
gs = viu + 1i*viv;
ws = viwu + 1i*viwv;
as = gs-ws;

clear viv viu viwu viwv

%%
% compute daily average
vidd = splitapply(@nansum,vid,g.day_id);
gsd=splitapply(@nansum,gs.*vid,g.day_id) ./ vidd;
wsd=splitapply(@nansum,ws.*vid,g.day_id) ./ vidd;
wsdw = splitapply(@nanmean,ws,g.day_id);
asd=splitapply(@nansum,as.*vid,g.day_id) ./ vidd;

% compute seasonal average of every year
[season,~,season_id] = unique([g.day_doy>180 year(g.day)],'rows');
vids = splitapply(@(x)sum(x, 1,'omitnan'),vidd,season_id);
gss=splitapply(@(x)sum(x, 1,'omitnan'),gsd.*vidd,season_id) ./ vids;
wss=splitapply(@(x)sum(x, 1,'omitnan'),wsd.*vidd,season_id) ./ vids;
wssw = splitapply(@(x)mean(x, 1,'omitnan'),wsdw,season_id);
ass=splitapply(@(x)sum(x, 1,'omitnan'),asd.*vidd,season_id) ./ vids;

% compute day of year average
viddoy = splitapply(@(x)sum(x, 1,'omitnan'),vidd,g.day_doy);
gsdoy=splitapply(@(x)sum(x, 1,'omitnan'),gsd.*vidd,g.day_doy) ./ viddoy;
wsdoy=splitapply(@(x)sum(x, 1,'omitnan'),wsd.*vidd,g.day_doy) ./ viddoy;
wsdoyw = splitapply(@(x)mean(x, 1,'omitnan'),wsdw,g.day_doy);
asdoy=splitapply(@(x)sum(x, 1,'omitnan'),asd.*vidd,g.day_doy) ./ viddoy;

%% 
parm={'gss','wss','ass','wssw'};
parm_name={'Groundspeed','Windspeed','Airspeed','Available Windspeed'};
figure;  tiledlayout(numel(parm),1,'TileSpacing','tight','Padding','tight')
 
col_s = [51 137 86;229 103 7]/255;
x_axis = (0:.01:20)'; pdf_width=1.2;
for i_p=1:numel(parm)
    nexttile; hold on;
    p = eval(parm{i_p});
    for i_s=0:1
        id_s=find(season(:,1)==i_s);
    
        pdf_p=nan(numel(x_axis),numel(id_s));
        for i_y=1:numel(id_s)
            if ~all(isnan(abs(p(id_s(i_y),:))))
                pdf_p(:,i_y) = ksdensity(abs(p(id_s(i_y),:)),x_axis);
            end
        end
        PDF_p = [pdf_p; -flipud(pdf_p)];
        X_axis = [x_axis;flipud(x_axis)];
        
        for i_y=1:numel(id_s)
            fill(season(id_s(i_y),2)+.25+.5*i_s+PDF_p(:,i_y)*pdf_width,X_axis,col_s(i_s+1,:),'FaceAlpha',0.5)
        end
        % plot(season(id_s,2)+.25+.5*i_s,abs(p(id_s,:)),'.','Color',col_s(i_s+1,:))
        plot(season(id_s,2)+.25+.5*i_s, nanmean(abs(p(id_s,:)),2),'.','markersize',15,'Color',col_s(i_s+1,:))
        yline(nanmean(abs(p(id_s,:)),'all'),'Color',col_s(i_s+1,:))
    end
    xlim([1995 2021]); box on;
    ylabel([parm_name{i_p} " [m/s]"])
end

%% 
% parm={'gss','wss','ass'}; parm_name={'Groundspeed','Windspeed','Airspeed'};
figure('position',[0 0 1200 800]); 
tiledlayout(numel(parm),2,'TileSpacing','tight','Padding','tight')

for i_p=1:numel(parm)
    p = eval(parm{i_p});
    for i_s=0:1
        id_s=find(season(:,1)==i_s);
        nexttile; hold on; box on
        
        borders('states','edgecolor','k','facecolor',[.8 .8 .8])
        scatter(radar.lon,radar.lat,sum(vids(id_s,:),1)/1e4,mean(abs(p(id_s,:)),1,'omitnan'),'filled','MarkerEdgeColor','k');
        quiver(radar.lon,radar.lat,mean(real(p(id_s,:)),'omitnan')',mean(imag(p(id_s,:)),'omitnan')','k','LineWidth',1)
        axis equal tight;
        axis([min(radar.lon)-1 max(radar.lon)+1 min(radar.lat)-1 max(radar.lat)+1 ]);
        title(parm_name{i_p})
        caxis([0 15]); 
        %c=colorbar;% c.Label.String="Spring-Autumn [m/s]";
        xticklabels('');yticklabels('');
        colormap(crameri('batlow'))
    end
end

%% Animated Map

p = wsdoy;
tt = string(g.day_doy);
figure('position',[0 0 1200 800]); 
tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
borders('states','edgecolor','k','facecolor',[.8 .8 .8])
axis equal tight;
axis([min(radar.lon)-1 max(radar.lon)+1 min(radar.lat)-1 max(radar.lat)+1 ]);
caxis([0 15]); 
%c=colorbar;% c.Label.String="Spring-Autumn [m/s]";
xticklabels('');yticklabels('');
colormap(crameri('batlow'))

id_s=1;
s=scatter(radar.lon,radar.lat,sum(vids(id_s,:),1)/1e4,mean(abs(p(id_s,:)),1,'omitnan'),'filled','MarkerEdgeColor','k');
q=quiver(radar.lon,radar.lat,mean(real(p(id_s,:)),'omitnan')',mean(imag(p(id_s,:)),'omitnan')','k','LineWidth',1);
t=title(parm_name{i_p} + string(tt(id_s)));

for id_s=1:size(p,1)
    
    title(parm_name{i_p})
    scatter(radar.lon,radar.lat,sum(vids(id_s,:),1)/1e4,mean(abs(p(id_s,:)),1,'omitnan'),'filled','MarkerEdgeColor','k');
    quiver(radar.lon,radar.lat,mean(real(p(id_s,:)),'omitnan')',mean(imag(p(id_s,:)),'omitnan')','k','LineWidth',1)
    title(parm_name{i_p} + string(tt(id_s)))

end
end
%%
col_gwa = [  82, 43, 41;204, 164, 59;17, 138, 178;67, 100, 54; 228, 87, 46;204, 164, 59]/255;

figure('position',[0 0 1000 600]);hold on; box on; grid on
parm={'gsdoy','wsdoy','asdoy','wsdoyw'};
for i_p=1:numel(parm)
    p = eval(parm{i_p});
    plot(datetime(0,1,1:366),smooth(mean(abs(p(:,radar.lat<40)),2,'omitnan')),'Color',col_gwa(i_p,:),'LineWidth',2)
end
legend(parm_name)

%%

% group radar
lat_gr=[0 30 35 40 45 Inf];
radar.group(:)=nan;
for i_lat=1:(numel(lat_gr)-1)
    radar.group(radar.lat>lat_gr(i_lat)&radar.lat<lat_gr(i_lat+1))=i_lat;
end
col_lat = parula(numel(lat_gr)-1);

parm={'gsdoy','wsdoy','asdoy'};


figure('position',[0 0 800 600]); tiledlayout(4,1,'TileSpacing','tight','Padding','tight')
for i_p=1:numel(parm)
    nexttile; 
    hold on; box on; grid on
    p = eval(parm{i_p});
    for i_lat=1:(numel(lat_gr)-1)
        plot(datetime(0,1,1:366),smooth(mean(abs(p(:,radar.group==i_lat)),2,'omitnan')),'color',col_lat(i_lat,:),'linewidth',2)
    end
    ylabel(parm_name{i_p})
end

nexttile;
hold on; box on   
borders('states','edgecolor','k','facecolor',[.8 .8 .8])
scatter(radar.lon,radar.lat,100,radar.group,'filled','MarkerEdgeColor','k');
axis equal tight;
axis([min(radar.lon)-1 max(radar.lon)+1 min(radar.lat)-1 max(radar.lat)+1 ]);
xticklabels('');yticklabels('');


%%

