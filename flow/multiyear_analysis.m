
load('data/density/inference-trans.mat','g','radar')
addpath(genpath('functions/'))
col2=brewermap([],'Paired');
clmap = brewermap([],'Spectral');
clmapd = crameri('berlin'); %hex2rgb(['#9e0142'; '#a00342'; '#a20543'; '#a40843'; '#a60a44'; '#a90d44'; '#ab0f45'; '#ad1245'; '#af1446'; '#b21746'; '#b41947'; '#b61c47'; '#b81e48'; '#bb2148'; '#bd2349'; '#bf2649'; '#c1284a'; '#c42b4b'; '#c62d4b'; '#c8304c'; '#ca324c'; '#cd354d'; '#cf374d'; '#d13a4e'; '#d33c4e'; '#d03c4d'; '#c73a4a'; '#bf3746'; '#b63543'; '#ad3240'; '#a5303d'; '#9c2d3a'; '#932b36'; '#8b2833'; '#822530'; '#79232d'; '#712029'; '#681e26'; '#5f1b23'; '#561920'; '#4e161d'; '#451419'; '#3c1116'; '#340f13'; '#2b0c10'; '#220a0c'; '#1a0709'; '#110506'; '#080203'; '#000000'; '#020507'; '#040b0f'; '#061017'; '#08161e'; '#0a1b26'; '#0c212e'; '#0e2635'; '#102c3d'; '#123145'; '#14374d'; '#163d54'; '#18425c'; '#1a4864'; '#1c4d6c'; '#1e5373'; '#20587b'; '#225e83'; '#24638a'; '#266992'; '#286f9a'; '#2a74a2'; '#2c7aa9'; '#2e7fb1'; '#3085b9'; '#3286bc'; '#3484bb'; '#3682ba'; '#387fb9'; '#3a7db8'; '#3b7bb6'; '#3d78b5'; '#3f76b4'; '#4174b3'; '#4371b2'; '#446fb1'; '#466db0'; '#486aaf'; '#4a68ae'; '#4c66ad'; '#4d63ab'; '#4f61aa'; '#515fa9'; '#535ca8'; '#555aa7'; '#5658a6'; '#5855a5'; '#5a53a4'; '#5c51a3'; '#5e4fa2']);

bd = load('borderdata.mat');
Anan = cellfun(@(x) [x(:);NaN],bd.lat(247:302),'un',0);
latb = cell2mat(Anan(:));
Anan = cellfun(@(x) [x(:);NaN],bd.lon(247:302),'un',0);
lonb = cell2mat(Anan(:));

%% Export daily map

for i_y=2010:2021
    load(['data/flow/est_' num2str(i_y) '.mat'],'Fd','gext');  
    tmp = -(Fd.takingoff+Fd.landing);

    for i_m=1:12
        id_m = find(month(gext.day)==i_m);
        figure('position',[0 0 2200 1200]); tiledlayout('flow','TileSpacing','none','Padding','none')
        colormap(clmap)
        for i_d=1:31
            nexttile; hold on; xticks([]); yticks([]); box on;set(gca,'Visible','off'); set(gca, 'color', 'none');
            if i_d<=numel(id_m)
                axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
                mlabel off; plabel off; gridm off; framem; tightmap; box off;
                 setm(gca,'frame','off'); 
                im2=geoshow(g.LAT,g.LON,smooth2a(Fd.takingoff(:,:,id_m(i_d)) + Fd.landing(:,:,id_m(i_d)),1),'DisplayType','surface');
                plot3m(latb,lonb,1000000000,"k",LineWidth=.1);
                textm(25,-120,datestr(gext.day(id_m(i_d))),'FontSize',16,'verticalalignment',"bottom");
                clim([-1 1]*2e5)
            end
        end
        exportgraphics(gcf, "figures/flow/daily_maps/"+num2str(i_y)+"_"+num2str(i_m)+".png")
    end
    close all
end

%% Accumulation
Y=11; K=1; D=366;
Ts_leaving_day=nan(D,K,Y);
Ts_entering_day=nan(D,K,Y);
t_axis = datenum(datetime('2010-1-1') + (0:(D-1)))';

for i_y=2010:2021
    load(['data/flow/est_' num2str(i_y) '.mat'],'Ts','gext')
    Ts_leaving_day(1:gext.nat,:,i_y-1999) = Ts.leaving_day;
    Ts_entering_day(1:gext.nat,:,i_y-1999) = Ts.entering_day;
end

accum = cumsum(Ts_leaving_day+Ts_entering_day)/10^9;

f_col = @(y) interp1(linspace(2000,2021,height(colormap)),colormap,y);
figure('position',[0 0 1000 600]); box on; grid on; hold on;
for i_y=2010:2021
    tmp_q10 = quantile(accum(:,:,i_y-1999),.1,2);
    tmp_q90 = quantile(accum(:,:,i_y-1999),.9,2);
    
    plot(t_axis,mean(accum(:,:,i_y-1999),2),'color',f_col(i_y),'LineWidth',2);
    fill([t_axis ; flipud(t_axis)] , [tmp_q10 ;flipud(tmp_q90)],f_col(i_y),'FaceAlpha',0.1,'EdgeAlpha',0);
end
datetick('x'); axis tight; ylabel('Number of birdss (billions)')




t_axis_sa = t_axis<datenum(datetime('2010-7-1'));

figure('position',[0 0 1000 600]); 
ha=tight_subplot(1,2);
axes(ha(1)); box on; grid on; hold on; 
for i_y=2000:2021
    accum_s = accum(t_axis_sa,:,i_y-1999)./ accum(find(t_axis_sa,1,'last'),:,i_y-1999);
    plot(t_axis(t_axis_sa),mean(accum_s,2),'color',f_col(i_y),'LineWidth',2);
    fill([t_axis(t_axis_sa) ; flipud(t_axis(t_axis_sa))] , [quantile(accum_s,.1,2) ;flipud(quantile(accum_s,.9,2))],f_col(i_y),'FaceAlpha',0.1,'EdgeAlpha',0);
end
datetick('x'); axis tight; ylabel('Number of birdss')

axes(ha(2)); box on; grid on; hold on; 
for i_y=2000:2021
    accum_s = accum(~t_axis_sa,:,i_y-1999)-accum(find(~t_axis_sa,1),:,i_y-1999);
    accum_s = -accum_s ./ accum_s(end-1,:,:);
    accum_s(isnan(accum_s))=-1;
    plot(t_axis(~t_axis_sa),mean(accum_s,2),'color',f_col(i_y),'LineWidth',2);
    fill([t_axis(~t_axis_sa) ; flipud(t_axis(~t_axis_sa))] , [quantile(accum_s,.1,2) ;flipud(quantile(accum_s,.9,2))],f_col(i_y),'FaceAlpha',0.1,'EdgeAlpha',0);
end
datetick('x'); axis tight; ylabel('Number of birdss')

%% transect/area

gext.mask_out=movmean(movmean(~gext.mask_water,[1 1],2),[1 1],1)>0;
gext.mask_out(~gext.mask_water)=false;

% IN/OUT: North South
mask_NS = double(gext.mask_out);
mask_NS(gext.mask_out&gext.lat<41) = 2;
mask_NS(gext.mask_out&gext.lat<45&gext.lon'>-75) = 2;

% imagesc(mask_NS) % 1: north, 2: south
    
% IN/OUT
sect=imread('data/flow/sectors_classified.tif');
transect_name={'east', 'northeast', 'north','west', 'mexico','gulf'};


% Area States
bd = load('borderdata.mat'); 
st.lat = bd.lat(247:302); 
st.lon = bd.lon(247:302); 
st.name = bd.places(247:302);
mask_states = nan(g.nlat,g.nlon);
for i_st=1:numel(st.lat)
    tmp = reshape(inpolygon(g.LON(:),g.LAT(:),st.lon{i_st},st.lat{i_st}),size(g.LAT));
    mask_states(tmp) = i_st;
end
states_id = unique(mask_states(~isnan(mask_states)));
% imagesc(mask_states)

% Area BCR
BCR = shaperead('data/flow/bcr_terrestrial_shape/BCR_Terrestrial_master.shp');
mask_bcr = nan(g.nlat,g.nlon);
for i_st=1:numel(BCR)
    tmp = reshape(inpolygon(g.LON(:),g.LAT(:),BCR(i_st).X,BCR(i_st).Y),size(g.LAT));
    mask_bcr(tmp) = BCR(i_st).BCR;
end
% imagesc(mask_bcr); mapshow(BCR)
bcr_id = unique(mask_bcr(~isnan(mask_bcr)));

%% Seasonal flux transect

FsInoutNS=nan(2,2,K,21);
FsInoutSect=nan(numel(transect_name),2,K,21);
FsAreaSt = nan(numel(states_id),2,2,K,21);
FsAreaBcr = nan(numel(bcr_id),2,2,K,21);

for i_y=2010:2021
    i_y
    load(['data/flow/est_' num2str(i_y) '.mat']) 
    midYear = datetime(['1-July-' num2str(i_y)]);

    Fdinout=Fd.entering+Fd.leaving;
    for i_s=1:2
        tmp = sum(reshape(Fdinout(repmat(mask_NS==i_s,1,1,gext.nat,K)),[],gext.nat,K));
        FsInoutNS(i_s,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsInoutNS(i_s,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
    end
    for i_s=1:numel(transect_name)
        tmp =  sum(reshape(Fdinout(repmat(sect==i_s,1,1,gext.nat,K)),[],gext.nat,K));
        FsInoutSect(i_s,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsInoutSect(i_s,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
    end

    for i_s=1:numel(states_id)
        tmp =  sum(reshape(Fd.landing(repmat(mask_states==states_id(i_s),1,1,gext.nat,K)),[],gext.nat,K),'omitnan');
        FsAreaSt(i_s,1,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsAreaSt(i_s,1,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
        tmp =  sum(reshape(Fd.takingoff(repmat(mask_states==states_id(i_s),1,1,gext.nat,K)),[],gext.nat,K),'omitnan');
        FsAreaSt(i_s,2,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsAreaSt(i_s,2,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
    end

    for i_s=1:numel(bcr_id)
        tmp =  sum(reshape(Fd.landing(repmat(mask_bcr==bcr_id(i_s),1,1,gext.nat,K)),[],gext.nat,K),'omitnan');
        FsAreaBcr(i_s,1,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsAreaBcr(i_s,1,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
        tmp =  sum(reshape(Fd.takingoff(repmat(mask_bcr==bcr_id(i_s),1,1,gext.nat,K)),[],gext.nat,K),'omitnan');
        FsAreaBcr(i_s,2,1,:,i_y-1999) = sum(tmp(:,midYear>gext.day,:),2,'omitnan');
        FsAreaBcr(i_s,2,2,:,i_y-1999) = sum(tmp(:,midYear<gext.day,:),2,'omitnan');
    end

end

% save('data/flow/Fs','FsInoutNS','FsInoutSect','FsAreaSt','FsAreaBcr')
load('data/flow/Fs')

%% Visualization

% FsInoutNS: North/South, Spring/Autum, sim, years
% FsInoutSect: Transect, Spring/Autum, sim, years

yy=2000:2021;

% figure('position',[0 0 1000 600]);  hold on; grid on
% tmp = permute(reshape(FsInoutNS,4,K,21),[3 1 2]);
% bar(yy,mean(tmp,3))

figure('position',[0 0 1000 600]);  hold on; grid on
s_label={'North','South'}; t_label={'Spring','Autumn'};u=1; clear l
for i_s=2:-1:1
    for i_t=1:2
        tmp = permute(FsInoutNS(i_s,i_t,:,:),[4 3 1 2]);
        l(u)=bar(yy,mean(tmp,2),'displayname',[s_label{i_s} ' - ' t_label{i_t}]);
        % errorbar(yy,mean(tmp,2),quantile(tmp,.1,2),quantile(tmp,.9,2),'.k','LineWidth',2);
        errorbar(yy,mean(tmp,2),std(tmp,[],2),'.k','LineWidth',2);
        u=u+1;
    end
end
legend(l,'location','north','Orientation','horizontal');
ylim([-4.5 4.5]*10^9)

% compute the ratio like Adriaan
figure('position',[0 0 1000 600]);  hold on; grid on; clear h
for i_s=2:-1:1
    tmp = permute(FsInoutNS(i_s,2,:,:) ./ FsInoutNS(i_s,1,:,:),[4 3 1 2]);
    h(i_s)=errorbar(yy,mean(abs(tmp),2),std(tmp,[],2),'o','LineWidth',2); 
    yline(nanmean(mean(abs(tmp),2)),'--','color',h(i_s).Color); 
end
legend(h,{'South','North'},'location','north','Orientation','horizontal');
title('Ratio Autunn/Spring - recrutement')


figure('position',[0 0 1000 600]);  hold on; grid on; clear h
for i_s=2:-1:1
    tmp = permute(FsInoutNS(i_s,1,:,2:end) ./ FsInoutNS(i_s,2,:,1:end-1),[4 3 1 2]);
    h(i_s)=errorbar(yy(1:end-1),mean(abs(tmp),2),std(tmp,[],2),'o','LineWidth',2,'displayname',[s_label{i_s} ' - s/a']); 
    yline(nanmean(mean(abs(tmp),2)),'--','color',h(i_s).Color)
end
legend(h,{'South','North'},'location','north','Orientation','horizontal');
title('Ratio Spring/Autunn - wintering')

% Area

% FsAreaSt: States, Landing/Takingoff, Spring/Autum, sim, years
% FsAreaBcr: BCR, Landing/Takingoff, Spring/Autum, sim, years


figure('position',[0 0 1000 600]);  hold on; grid on
for i_s=1:numel(states_id)
    tmp = permute(nansum(FsAreaSt(i_s,:,:,:,:),2),[5 4 3 1 2]);
    tmp = tmp(:,:,2)./ tmp(:,:,1);
    %h=errorbar(yy,mean(abs(tmp),2),std(tmp,[],2),'o','LineWidth',2,'displayname',[st.name{states_id(i_s)} ' - a/s']); 
    plot(yy,mean(abs(tmp),2),'LineWidth',2,'displayname',[st.name{states_id(i_s)} ' - a/s']); 
end

%% Map of MVT
MVT_season=nan(g.nlat,g.nlon,21,2);
for i_y=2010:2021
    i_y
    load(['data/flow/est_' num2str(i_y) '.mat'],'MVT_day','gext') 
    midYear = datetime(['1-July-' num2str(i_y)]);
    MVT_season(:,:,i_y-1999,1) = sum(MVT_day(:,:,midYear>gext.day),3,'omitnan');
    MVT_season(:,:,i_y-1999,2) = sum(MVT_day(:,:,midYear<gext.day),3,'omitnan');
end

figure('position',[0 0 1000 600]); 
tiledlayout('flow','TileSpacing','tight','Padding','tight')
nexttile; 
imagesc(g.lon,g.lat,(mean(MVT_season(:,:,:,1),3,'omitnan')),'alphadata',~g.mask_water); 
axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);colorbar;
title('MVT Spring')

nexttile; 
imagesc(g.lon,g.lat,(mean(MVT_season(:,:,:,2),3,'omitnan')),'alphadata',~g.mask_water); 
axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]); colorbar;
title('MVT Autumn')


%% Map Seasonal
Fs_takingoff=nan(g.nlat,g.nlon,21,2); Fs_landing=Fs_takingoff;
for i_y=2010:2021
    load(['data/flow/est_' num2str(i_y) '.mat'],'Fd','gext')
    i_y
    midYear = datetime(['1-July-' num2str(i_y)]);

    Fs_takingoff(:,:,i_y-1999,1) = mean(sum(Fd.takingoff(:,:,gext.day<midYear,:),3,'omitnan'),4,'omitnan');
    Fs_landing(:,:,i_y-1999,1) = mean(sum(Fd.landing(:,:,gext.day<midYear,:),3,'omitnan'),4,'omitnan');

    Fs_takingoff(:,:,i_y-1999,2) = mean(sum(Fd.takingoff(:,:,gext.day>midYear,:),3,'omitnan'),4,'omitnan');
    Fs_landing(:,:,i_y-1999,2) = mean(sum(Fd.landing(:,:,gext.day>midYear,:),3,'omitnan'),4,'omitnan');
end

% save('data/flow/Fs_map','Fs_takingoff','Fs_landing')
load('data/flow/Fs_map')
clmap = brewermap([],'Spectral');


% tmp = Fs_takingoff(:,:,:,1)./g.area;
%tmp = Fs_landing(:,:,:,2)./g.area; 
tmp = (Fs_takingoff(:,:,:,1)+Fs_landing(:,:,:,1));%./g.area;
%tmp = (Fs_takingoff(:,:,:,2)+Fs_landing(:,:,:,2))./(Fs_takingoff(:,:,:,1)+Fs_landing(:,:,:,1));

figure('position',[0 0 1000 600]); 
tiledlayout('flow','TileSpacing','none','Padding','none')
for i_y=2010:2021
    nexttile;hold on; title(num2str(i_y)); box on; xticks([]); yticks([])
    imagesc(g.lon,g.lat,smooth2a(tmp(:,:,i_y-1999),2),'alphadata',~g.mask_water);
    axis tight equal;
    %colorbar;
    % caxis([-1 1]*max(abs(tmps(:)))); %colormap(gca,clmap);
    % caxis([0 1]*max(abs(tmp(:))));
    caxis([-1 1]*3e6)
end

%% Average accross all year
figure('position',[0 0 800 800]); 
tiledlayout('flow','TileSpacing','tight','Padding','tight')
nexttile;  box on; xticks([]); yticks([])
imagesc(g.lon,g.lat,mean(Fs_takingoff(:,:,:,1),3,'omitnan')./g.area,'alphadata',~g.mask_water); 
axis tight equal;caxis([0 35000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Taking-off Spring'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
imagesc(g.lon,g.lat,mean(Fs_takingoff(:,:,:,2),3,'omitnan')./g.area,'alphadata',~g.mask_water); 
axis tight equal;caxis([0 35000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Taking-off Autumn'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
imagesc(g.lon,g.lat,mean(-Fs_landing(:,:,:,1),3,'omitnan')./g.area,'alphadata',~g.mask_water); 
axis tight equal;caxis([0 35000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Landing Spring'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
imagesc(g.lon,g.lat,mean(-Fs_landing(:,:,:,2),3,'omitnan')./g.area,'alphadata',~g.mask_water); 
axis tight equal;caxis([0 35000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Landing Autumn'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
tmp = mean(-Fs_takingoff(:,:,:,1)-Fs_landing(:,:,:,1),3,'omitnan')./g.area;
tmps = imgaussfilt(tmp,2,'FilterSize',31);
imagesc(g.lon,g.lat,tmps,'alphadata',~g.mask_water); axis tight equal; caxis([-4000 4000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Diff Spring'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
tmp = mean(-Fs_takingoff(:,:,:,2)-Fs_landing(:,:,:,2),3,'omitnan')./g.area;
tmps = imgaussfilt(tmp,2,'FilterSize',15);
imagesc(g.lon,g.lat,tmps,'alphadata',~g.mask_water); 
axis tight equal;caxis([-4000 4000]); set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Diff Autumn'); xticks([]); yticks([])

nexttile;  box on; xticks([]); yticks([])
tmp = mean(-Fs_takingoff(:,:,:,1)-Fs_landing(:,:,:,1),3,'omitnan')./mean(-Fs_landing(:,:,:,1),3,'omitnan');
tmp(isnan(tmp))=0;
tmps = imgaussfilt(tmp,2,'FilterSize',15);
imagesc(g.lon,g.lat,(tmps),'alphadata',~g.mask_water); 
axis tight equal;  set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Diff/landing Spring'); xticks([]); yticks([])

nexttile;  box on
tmp = mean(-Fs_takingoff(:,:,:,2)-Fs_landing(:,:,:,2),3,'omitnan')./mean(-Fs_landing(:,:,:,2),3,'omitnan');
tmp(isnan(tmp))=0;
tmps = imgaussfilt(tmp,2,'FilterSize',15);
imagesc(g.lon,g.lat,(tmps),'alphadata',~g.mask_water); 
axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);
title('Diff/landing Autumn'); xticks([]); yticks([])


%% Accumulation Map
Fs_accum=nan(g.nlat,g.nlon,366,21);
gext_dayY=NaT(366,21);
for i_y=2010:2021
    load(['data/flow/est_' num2str(i_y) '.mat'],'Fd','gext')
    Fs_accum(:,:,1:size(Fd.takingoff,3),i_y-1999) = mean(Fd.takingoff,4,'omitnan')+mean(Fd.landing,4,'omitnan');
    gext_dayY(1:size(Fd.takingoff,3),i_y-1999) = gext.day;
end


fig=figure('position',[0 0 1200 800]); tiledlayout('flow','TileSpacing','none','Padding','none'); nexttile;
im=imagesc(g.lon,g.lat,Fs_accum(:,:,1,i_y-1999)./g.area,'alphadata',~g.mask_water); 
axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]); caxis([-300 300]); 
colormap(clmap); xticks([]); xticks([]); set(gca, 'color', 'none'); set(gca,'Visible','off')
t = text(-120,25,datestr(gext_dayY(i_t,i_y-1999)),'FontSize',40);

v = VideoWriter('figures/flow/accumulation/accumulation_daily.mp4','MPEG-4');
v.Quality = 50;
v.FrameRate = 5;
open(v)

for i_y=2010:2021
    for i_t=1:366
        try
            im.CData = smooth2a(Fs_accum(:,:,i_t,i_y-1999)./g.area,2);
            t.String=datestr(gext_dayY(i_t,i_y-1999));
            drawnow
            % pause(.1)
            writeVideo(v,getframe(fig))
        end
    end
end
close(v)

%% Daily change per area 

FdInoutNS=nan(2,366,21);
FdInoutSect=nan(numel(transect_name),366,21);
FdAreaSt = nan(numel(states_id),2,366,21);
FdAreaBcr = nan(numel(bcr_id),2,366,21);

for i_y=2010:2021
    i_y
    load(['data/flow/est_' num2str(i_y) '.mat']) 

    Fdinout=mean(Fd.entering+Fd.leaving,4,'omitnan');
    for i_s=1:2
        FdInoutNS(i_s,1:size(Fdinout,3),i_y-1999) = sum(reshape(Fdinout(repmat(mask_NS==i_s,1,1,gext.nat)),[],gext.nat));
    end
    for i_s=1:numel(transect_name)
        FdInoutSect(i_s,1:size(Fdinout,3),i_y-1999) =  sum(reshape(Fdinout(repmat(sect==i_s,1,1,gext.nat)),[],gext.nat));
    end
    tmp1 = mean(Fd.landing,4,'omitnan');
    tmp2 = mean(Fd.takingoff,4,'omitnan');
    for i_s=1:numel(states_id)
        FdAreaSt(i_s,1,1:size(tmp1,3),i_y-1999) = sum(reshape(tmp1(repmat(mask_states==states_id(i_s),1,1,gext.nat)),[],gext.nat),'omitnan');
        FdAreaSt(i_s,2,1:size(tmp2,3),i_y-1999) = sum(reshape(tmp2(repmat(mask_states==states_id(i_s),1,1,gext.nat)),[],gext.nat),'omitnan');
    end

    for i_s=1:numel(bcr_id)
        FdAreaBcr(i_s,1,1:size(tmp1,3),i_y-1999) = sum(reshape(tmp1(repmat(mask_bcr==bcr_id(i_s),1,1,gext.nat)),[],gext.nat),'omitnan');
        FdAreaBcr(i_s,2,1:size(tmp2,3),i_y-1999) =  sum(reshape(tmp2(repmat(mask_bcr==bcr_id(i_s),1,1,gext.nat)),[],gext.nat),'omitnan');
    end

end

% save('data/flow/Fd','FdInoutNS','FdInoutSect','FdAreaSt','FdAreaBcr')
load('data/flow/Fd')

%%
%i_s=find(states_id==find(st.name=="Hawai"));

for i_s=1:numel(states_id)
    figure('position',[0 0 1200 600]); tiledlayout('flow','TileSpacing','none','Padding','none'); nexttile; hold on;grid on; box on;
    FdAreaStrs = -squeeze(sum(FdAreaSt(i_s,:,:,:),2));
    x = datetime(2000,1,1:366);
    % errorbar(1:366,mean(FdAreaStrs,2),std(FdAreaStrs,[],2));
    patch([x fliplr(x)],[smooth(quantile(FdAreaStrs,.4,2),10); flipud(smooth(quantile(FdAreaStrs,.6,2),10))]','k','facealpha',.2,'edgealpha',.1)
    patch([x fliplr(x)],[smooth(quantile(FdAreaStrs,.25,2),10); flipud(smooth(quantile(FdAreaStrs,.75,2),10))]','k','facealpha',.2,'edgealpha',.1)
    axis tight; y_lim=ylim;
    patch([x fliplr(x)],[smooth(quantile(FdAreaStrs,.10,2),10); flipud(smooth(quantile(FdAreaStrs,.90,2),10))]','k','facealpha',.2,'edgealpha',.1)
    text(x(1)+2,max(abs(y_lim)),st.name{states_id(i_s)},'fontsize',36,'verticalalignment','top')
    scatter(repmat((x)',22,1),FdAreaStrs(:),20,'k','filled','MarkerFaceAlpha',.2)
    bar(x,smooth(nanmedian(FdAreaStrs,2),10),'facecolor','k')
    datetick('x', 'mmm'); axis tight; ylim([-1 1]*max(abs(y_lim)));
    ylabel('Change on the ground (Nb. birds)',FontSize=18);
    exportgraphics(gcf, "figures/flow/states/daily_accum_"+st.name{states_id(i_s)}+".png")
end

%% Along North South

Fd_takingoff=nan(g.nlat,366,21);
Fd_landing=nan(g.nlat,366,21);
Fd_accum=nan(g.nlat,366,21);
gext_dayY=NaT(366,21);


for i_y=2010:2021
    i_y
    load(['data/flow/est_' num2str(i_y) '.mat'],'Fd','gext')
    tmp1 = nan(size(Fd.takingoff)); tmp2=tmp1;
    for i_t=1:size(Fd.takingoff,3)
        tmp1(:,:,i_t) = smooth2a(double(Fd.takingoff(:,:,i_t)./g.area),2);
        tmp2(:,:,i_t) = smooth2a(double(Fd.landing(:,:,i_t)./g.area),2);
    end
    id_lon = g.lon>-105;
    %id_lon = g.lon>-98 & g.lon<-92;%105
    Fd_takingoff(:,1:size(Fd.takingoff,3),i_y-1999) = squeeze(mean(tmp1(:,id_lon,:),2,'omitnan'));
    Fd_landing(:,1:size(Fd.takingoff,3),i_y-1999) = squeeze(mean(tmp2(:,id_lon,:),2,'omitnan'));
    Fd_accum(:,1:size(Fd.takingoff,3),i_y-1999) = squeeze(mean(tmp1(:,id_lon,:)+tmp2(:,id_lon,:),2,'omitnan'));
    gext_dayY(1:size(Fd.takingoff,3),i_y-1999) = gext.day;
end


col = viridis(180);
figure; tiledlayout(2,6,'TileSpacing','none','Padding','tight');
for i=2010:2021
nexttile; hold on; box on; grid on;
for i_t=1:180
    plot(sum(-Fd_accum(:,185+(1:i_t),i-1999),2,'omitnan'),g.lat,'color',col(i_t,:));
end
if i>2010
    yticks([])
end
title(i)
 ylim([30 48]); xlim(1500*[-1 1])
end

figure; 
nexttile; hold on;
for i_t=1:180
    plot(g.lat,sum(-nanmean(Fd_accum(:,(1:i_t),:),3),2,'omitnan'),'color',col(i_t,:));
end
xlim([30 48]); box on; grid on; view([90 -90]); title('Spring')
nexttile; hold on;
for i_t=1:180
    plot(g.lat,sum(-nanmean(Fd_accum(:,185+(1:i_t),:),3),2,'omitnan'),'color',col(181-i_t,:));
end
xlim([30 48]); box on; grid on; view([90 -90]); title('Autumn')

nexttile; hold on;
for i_t=1:180
    plot(g.lat,sum(-nanmean(Fd_takingoff(:,(1:i_t),:),3),2,'omitnan'),'color',col(i_t,:));
    plot(g.lat,sum(-nanmean(Fd_landing(:,(1:i_t),:),3),2,'omitnan'),'color',col(i_t,:));
end
xlim([30 48]); box on; grid on; view([90 -90]);
nexttile; hold on;
for i_t=1:180
    plot(g.lat,sum(-nanmean(Fd_takingoff(:,185+(1:i_t),:),3),2,'omitnan'),'color',col(181-i_t,:));
    plot(g.lat,sum(-nanmean(Fd_landing(:,185+(1:i_t),:),3),2,'omitnan'),'color',col(181-i_t,:));
end
xlim([30 48]); box on; grid on; view([90 -90]); 

