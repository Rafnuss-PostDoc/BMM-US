
%% ----- INTRA-NIGHT

% Load data
load('data/density/inference-trans.mat','g','radar','thr_nnt','trans')
load('data/density/inference-daily-intra.mat','vidTSi','itrend')
load('data/density/inference-intra','icov')
load('data/density/estsim-daily')
addpath('./functions/'); 



%% ------- SERVER
% Number of simulation
sim_n=0;

% conditional value
condVal=vidTSi;

% plotit
plotit=false;

for i_y=2021:-1:1995
    % find the index of the datapoint to estimate in that year
    % Note that I'm using g.day(g.day_id) to get all the datapoint point whose mid-night belong to that year (thus including the datapoint from the
    % night of new year at the begining
    idt_y_tmp = find(year(g.day(g.day_id))==i_y);

    % Compute the gridded NNT
    NNTg_tmp = twilightNNT(g.time(idt_y_tmp), double(g.LON(~g.mask_water)), double(g.LAT(~g.mask_water)))';
    NNTg_tmp(NNTg_tmp<thr_nnt(1)|NNTg_tmp>thr_nnt(2))=nan;
    % limit the point to simulated with those during the night
    id_night = ~all(isnan(NNTg_tmp));
    idt_y = idt_y_tmp(id_night);
    NNTg = NNTg_tmp(:,id_night);

    % Unconditional simulation
    % Resi=simulationIntraUnconditional('d', g, radar, idt_y, sim_n, icov, plotit); % 20min for 12 sim , 10min for 2sim
    % save(['data/density/simuncond-intra_' num2str(i_y)],'Resi','icov','i_y','i_yy','-v7.3') % 13Go
    Resi = [];

    % Estimation and conditioning
    [Ei,Si,ResiC]=estimationSimulationIntraConditional('d', g, radar, idt_y, icov, condVal, Resi); % ~170 min for 12 sim,  for 2sim
    % clear Resi

    % 
    id_nnt = round(NNTg(:)/itrend.dnnt)+1/itrend.dnnt+1;
    id_nnt_doy = sub2ind(size(itrend.nnt_m), repelem(g.time_doy(idt_y),sum(~g.mask_water(:)),1), id_nnt);
    nnt_std_template = nan(size(NNTg));
    nnt_std_template(~isnan(id_nnt_doy)) = itrend.nnt_std(id_nnt_doy(~isnan(id_nnt_doy)));
    nnt_m_template = nan(size(NNTg));
    nnt_m_template(~isnan(id_nnt_doy)) = itrend.nnt_m(id_nnt_doy(~isnan(id_nnt_doy)));
   
    % Combine with daily scale
    Em = (Edm(:,g.day_id(idt_y)) + Ei.*nnt_std_template)+nnt_m_template;
    % Sm = Sdm(g:,g.day_id(i_yy)) + Sim;
    % back transform
    EmbT = trans.f_inv(Em);
    % Save
    save(['data/density/est_' num2str(i_y)],'EmbT','idt_y','-v7.3') % 800Mb
    % clear EmbT Em Eim Sim

    for k=1:sim_n
        ResCm = (ResdCm(:,g.day_id(idt_y),k) + ResiC(:,:,k).*nnt_std_template)+nnt_m_template;
        ResCmbT = trans.f_inv(ResCm);
        % save(['data/density/sim_' num2str(i_y) '_' num2str(k)],'ResCmbT','idt_y','-v7.3') % 800Mb
    end
end









%% Illustration
load('data/density/inference-trans.mat','vidTS')


i_y = 2021;
load(['data/density/est_' num2str(i_y)])

% vidT = vidti + vidtd(g.day_id,:);
% vidL = interp1(trans.denst_axis,trans.dens_axis,vidT,'pchip');
% vidLd = splitapply(@(x) mean(x,'omitnan'),vidL,g.day_id);

sfv = EmbT;
ptv = trans.f_inv(vidTS);

%find(any(~isnan(Em)))
i_s=1584;%10050;%numel(i_yy)-15;
figure('position',[0 0 1600 800]); 
tiledlayout('flow','TileSpacing','tight','Padding','tight')
for u=0:15
    nexttile; hold on; set(gca,'ydir','normal');
    tmp=nan(size(g.mask_water));
    tmp(~g.mask_water) = sfv(:,i_s+u); % ResC(:,i_s+u,k);% 
    imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    borders('states','k');
    id=~isnan(ptv(idt_y(i_s)+u,:)); % vidti  vidtditv
    scatter(radar.lon(id),radar.lat(id),100,ptv(idt_y(i_s)+u,id)','filled','MarkerEdgeColor','k');
    title(datestr(g.time(idt_y(i_s)+u)))
    axis([-125 -68 23 50]); colorbar;%caxis([-3 3])
end

%% 
figure; hold on; set(gca,'ydir','normal');
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = S(:,i_s+u); % std(ResC(:,i_s+u,:),[],3);%  
imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); colorbar;


