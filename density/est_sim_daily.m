% Script to produce estimation and simulation


%% ----- PART 1: Daily scale

% Load data
load('data/density/inference-trans.mat')
load('data/density/inference-daily.mat')
addpath('functions/'); 


%% ---- ON SERVER


%% Unconditional simulation
gn=15;
sim_n=0;
plotit=false;
Resd=simulationDailyUnconditional('d', g, gn, sim_n, dcov, plotit); %~15min
% save('data/density/simuncond-daily','Resd','dcov','gn','-v7.3')


%% Simulation and Estimation

Resd=[];
% load('data/density/simuncond-daily')
gn=15; % number of date to simulation together
gnm=4; % margin around those date
condVal = vidTSddt; % conditional value
[Ed,S,ResdC] = estimationSimulationDailyConditional('d', g, gn, gnm, radar, dcov, condVal, Resd);


%% Transform inverse
%Compute maps

daytrend_sf_std_LL = nan(g.nlm,366);
daytrend_sf_mean_LL = nan(g.nlm,366);
for i_m=1:numel(daytrend.sf_m)
    daytrend_sf_std_LL(:,daytrend.sf_m(i_m)) = daytrend.sf_std{i_m}([g.LON(~g.mask_water(:)) g.LAT(~g.mask_water(:))]);
    daytrend_sf_mean_LL(:,daytrend.sf_m(i_m)) = daytrend.sf_mean{i_m}([g.LON(~g.mask_water(:)) g.LAT(~g.mask_water(:))]);
end
daytrend_sf_std_LL=fillmissing(daytrend_sf_std_LL,'linear',2);
daytrend_sf_mean_LL=fillmissing(daytrend_sf_mean_LL,'linear',2);


% figure;
% for i_m=1:366
%     imagesc(daytrend_sf_std_LL(:,:,i_m),'alphadata',~g.mask_water)
%     title(i_m); pause(0.1)
% end

f_trans_inver = @(X) ( X ... 
    .* daytrend_sf_std_LL(:,g.day_doy) ...
    .* daytrend.std_doy(g.day_doy)' )...
    + daytrend_sf_mean_LL(:,g.day_doy)...
    + daytrend.f_dem(daytrend.g_dem(~g.mask_water))...
    + daytrend.m_doy(g.day_doy)';

Edm = f_trans_inver(Ed);
% Sdm = Sd .* abs(daytrend.sf_std([g.LON(~g.mask_water),g.LAT(~g.mask_water)])) .* abs(daytrend.f_std_doy(g.day_doy))';

ResdCm = f_trans_inver(ResdC);

save('data/density/estsim-daily','Edm','ResdCm','-v7.3') % 5Go





%% ----
load('data/density/estsim-daily')


%% Figure
load('data/density/inference-daily-intra.mat','vidTSd')


% sfv=ResdCm; % simulation
% sfv=Edm; % estimation
% ptv=vidTSd;

 % as stransformed back 
sfv=trans.f_inv(Edm);
ptv=trans.f_inv(vidTSd);

% Check histogram
% Not good, but always during the winter, not good filtering... 
figure; hold on;
histogram(sfv(:))
histogram(ptv(:))



k=1;
i_v1= find(year(g.day)==2010,1)+0;

figure('position',[0 0 1600 900]);
tiledlayout('flow','TileSpacing','tight','Padding','tight')
for i=1:1
    ha(i)=nexttile; hold on; set(gca,'ydir','normal'); axis equal tight
    tmp=nan(size(g.mask_water));
    tmp(~g.mask_water) = sfv(:,i_v1,k);
    im(i)=imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    borders('states','k');
    s(i)=scatter(radar.lon,radar.lat,100,ptv(i_v1,:),'filled','MarkerEdgeColor','k');
    axis([-125 -68 23 50]); colorbar; 
    %caxis([-3 3])
    caxis([10 500])
end
t= text(-120,25,datestr(g.day(i_v1)),'color','k','FontSize',40);

%v = VideoWriter('test');
%v.FrameRate=4;open(v);
for i_v=i_v1+(1:300)
    for i=1:1
        axes(ha(i));
        if i==1
            tmp(~g.mask_water) = sfv(:,i_v);
            title(datestr(g.day(i_v)),'Color','white')
        else
            tmp(~g.mask_water) = sfv(:,i_v,i);
        end
        im(i).CData = tmp;
        delete(s(i))
        tmp2 = ~isnan(ptv(i_v,:));
        s(i)=scatter(radar.lon(tmp2),radar.lat(tmp2),100,ptv(i_v,tmp2),'filled','MarkerEdgeColor','k');
    end
    t.String=datestr(g.day(i_v));
    refreshdata
    drawnow
    pause(0.2);
    % keyboard
    % writeVideo(v,getframe)
end
%close(v);

%% Export video for BMM-US
% Not working because need to use surfm to get transparency, but then,
% scatter plot below image... stupide... 

i_v1= find(year(g.day)==2010,1)+90;

figure('position',[0 0 1600 900]); tight_subplot(1,1,0,0,0);
%worldmap([23 50],[-125 -68]); 
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;

mlabel off; plabel off; gridm off; framem; tightmap; box off; 
set(gcf, 'color', 'none'); set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
tmp=nan(size(g.mask_water,1),size(g.mask_water,2),3);
clm = colormap; 
caxis([0 600]); c_axis = get(gca,'CLim');
Gs = @(x) round(fillmissing(interp1(linspace(c_axis(1),c_axis(end),size(clm,1)),1:size(clm,1),round(x)),'nearest'));
tmp(repmat(~g.mask_water,1,1,3)) = clm(Gs(sfv(:,i_v1,1)),:);
%im=surfm(g.lat,g.lon,tmp,'alphadata',~g.mask_water); 
im=geoshow(g.LAT,g.LON,tmp,'DisplayType','image');
% bordersm('states','k');
tmp2 = ~isnan(ptv(i_v1,:));
a_xis=axis;
s=scatterm(radar.lat(tmp2),radar.lon(tmp2),100,ptv(i_v1,tmp2),'filled','MarkerEdgeColor','w');

t= textm(25,-120,datestr(g.day(i_v1)),'color','w','FontSize',40);

for i_y=2019:2020

    filename="figures/estimation_daily/"+num2str(i_y)+"/"; %  estimation_daily
    mkdir(filename)
    % v = VideoWriter('data/estimation_daily.mp4','MPEG-4');
    % v.Quality = 95;
    % v.FrameRate = 2;
    % open(v)

    i_v1= find(year(g.day)==i_y,1);

    for i_v=i_v1+(1:366)%numel(g.day)
        tmp(repmat(~g.mask_water,1,1,3)) = clm(Gs(sfv(:,i_v,1)),:);
        t.String=datestr(g.day(i_v));
        im.CData = tmp;
        delete(s(i))
        tmp2 = ~isnan(ptv(i_v,:));
        tmp2(20)=false;
        s=scatterm(radar.lat(tmp2),radar.lon(tmp2),100,ptv(i_v,tmp2),'filled','MarkerEdgeColor','w');
        refreshdata
        drawnow
        axis(a_xis);
        %pause(0.2);

        % writeVideo(v,getframe(gcf))

       %[imind,cm] = frame2im(getframe(gcf)); 
    %   % Write to the GIF File 
    %   if i_v == i_v1 
    %       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    %   else 
    %       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    %   end 
        % export_fig(gca,[filename num2str(i_v) '.png'], '-transparent')
       exportgraphics(gca,filename+num2str(i_v)+ ".png")
    end
end

close(v)

%

% Emo = 10.^interp1(trans.denst_axis,trans.dens_axis,Em,'pchip');
% vidtdo = 10.^interp1(trans.denst_axis,trans.dens_axis,vidtd,'pchip');
% 
% 
% ii=[10 75 105];
% tmp=nan(size(g.LAT));
% 
% figure; ha=tight_subplot(3,numel(ii));
% for i=1:numel(ii)
% 
%     axes(ha(i)); hold on; set(gca,'ydir','normal'); axis equal tight
%     tmp(~g.mask_water) = E(:,ii(i));
%     im=imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
%     s=scatter(radar.lon,radar.lat,100,vidtddtv(ii(i),:),'filled','MarkerEdgeColor','k');
%     borders('states','k'); axis([-125 -68 23 50]); colorbar; title(datestr(g.day(ii(i))))
%     ylabel('Normalized density')
% 
%     axes(ha(i+numel(ii))); hold on; set(gca,'ydir','normal'); axis equal tight
%     im=imagesc(g.lon,g.lat,Em(:,:,ii(i)),'alphadata',~g.mask_water); 
%     s=scatter(radar.lon,radar.lat,100,vidtd(ii(i),:),'filled','MarkerEdgeColor','k');
%     borders('states','k'); axis([-125 -68 23 50]); colorbar; title(datestr(g.day(ii(i))))
%     ylabel('Transformed density')
% 
%     axes(ha(i+numel(ii)*2)); hold on; set(gca,'ydir','normal'); axis equal tight
%     im=imagesc(g.lon,g.lat,Emo(:,:,ii(i)),'alphadata',~g.mask_water); 
%     s=scatter(radar.lon,radar.lat,100,vidtdo(ii(i),:),'filled','MarkerEdgeColor','k');
%     borders('states','k'); axis([-125 -68 23 50]); colorbar; title(datestr(g.day(ii(i))))
%     ylabel('Density')
% end


%%

g2=g;
g2.day = g2.day(1:400);
g2.day_doy = g2.day_doy(1:400);

load('data/density/inference-daily.mat')
Resd_normal=simulationDailyUnconditional('d', g2, gn, sim_n, dcov, plotit); 
[~,~,ResC_normal] = estimationSimulationDailyConditional('d', g, gn, gnm, radar, dcov, condVal, Resd_normal);

load('data/test-nugget-daily-sim/inference-daily_max_nugget.mat')
Resd_max    =simulationDailyUnconditional('d', g2, gn, sim_n, dcov, plotit); 
[~,~,ResC_max] = estimationSimulationDailyConditional('d', g, gn, gnm, radar, dcov, condVal, Resd_max);

load('data/test-nugget-daily-sim/inference-daily_no_nugget')
Resd_no=simulationDailyUnconditional('d', g2, gn, sim_n, dcov, plotit); 
[~,~,ResC_no] = estimationSimulationDailyConditional('d', g, gn, gnm, radar, dcov, condVal, Resd_no);


seas_pattern = daytrend.sf_season([g.LON(~g.mask_water),g.LAT(~g.mask_water)]);
f_trans_inver2 = @(X) ( X ... 
    .* daytrend.sf_std([g.LON(~g.mask_water),g.LAT(~g.mask_water)])...
    .* daytrend.f_std_doy(g2.day_doy)' )...
    + daytrend.f_season(g2.day_doy)'.*seas_pattern...
    + daytrend.sf_mean([g.LON(~g.mask_water),g.LAT(~g.mask_water)])...
    + daytrend.f_dem(daytrend.g_dem(~g.mask_water))...
    + daytrend.f_doy(g2.day_doy)';

ResdCm_normal = f_trans_inver2(ResC_normal);
ResdCm_max = f_trans_inver2(ResC_max);
ResdCm_no = f_trans_inver2(ResC_no);


save('data/test-nugget-daily-sim/data.mat','ResdCm_normal','ResdCm_max','ResdCm_no','g2')
load('data/test-nugget-daily-sim/data.mat')
%%
k=1;
i_v=[290 258 281]; % 258 281 124

figure('position',[0 0 1600 900]);
ha=tight_subplot(numel(i_v),3);
for i=1:numel(i_v)
    for j=1:3
        axes(ha((i-1)*3+j));     
        hold on; set(gca,'ydir','normal'); axis equal tight
        tmp=nan(size(g.mask_water));
        if j==1
            tmp2 = ResdCm_normal;
            title(['Normal - ' datestr(g2.day(i_v(i)))])
        elseif j==2
            tmp2=ResdCm_max;
            title(['Max Nugget - ' datestr(g2.day(i_v(i)))])
        else
            tmp2=ResdCm_no;
            title(['No nugget - ' datestr(g2.day(i_v(i)))])
        end
        tmp(~g.mask_water) = tmp2(:,i_v(i),k);
        im(i)=imagesc(g.lon,g.lat,10.^tmp,'alphadata',~g.mask_water); 
        borders('states','k');
        id = ~isnan(vidtd(i_v(i),:));
        s(i)=scatter(radar.lon(id),radar.lat(id),100,10.^vidtd(i_v(i),id),'filled','MarkerEdgeColor','k');
        % text(radar.lon,radar.lat,radar.name)
        axis([-125 -68 23 50]); 
        %  axis([-90 -68 40 50]); 
        colorbar; %caxis([0 500])
    end
end


figure; hold on;
plot(g2.day,permute(sum(ResdCm_normal,1),[2 3 1]),'b','displayname','normal')
plot(g2.day,permute(sum(ResdCm_max,1),[2 3 1]),'r','displayname','max nugget')
plot(g2.day,permute(sum(ResdCm_no,1),[2 3 1]),'g','displayname','no nugget')
legend

figure; hold on;
plot(g2.day,permute(sum(10.^ResdCm_normal,1),[2 3 1]),'b','displayname','normal')
plot(g2.day,permute(sum(10.^ResdCm_max,1),[2 3 1]),'r','displayname','max nugget')
plot(g2.day,permute(sum(10.^ResdCm_no,1),[2 3 1]),'g','displayname','no nugget')
legend

