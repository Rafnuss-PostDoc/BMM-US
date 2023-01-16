
%% ----- INTRA-NIGHT

% Load data
load('data/density/inference-trans.mat','g','radar')
load('data/speed/inference-trans.mat')
load('data/speed/inference-intra.mat')
addpath('./functions/'); 





%% ----- SERVER

% Number of simulation
sim_n=0;

% plotit
plotit=false;

uv={'u','v'};

for i_y=2020%2000:5:2020

    % Get index of point to estimated from density
    load(['data/density/est_' num2str(i_y)],'idt_y');

    for i_uv = 1:2

        % conditional value
        condVal=victidt.(uv{i_uv});

        load(['data/speed/estsim-daily_' uv{i_uv}])

        % unconditional realization
        % Resi.(uv{i_uv}) = simulationIntraUnconditional(uv{i_uv}, g, radar, idt_y, sim_n, icov, plotit); % 20min
        Resi.(uv{i_uv}) = [];
        
        % estimation and conditional realization
        [Ei.(uv{i_uv}),Si.(uv{i_uv}),ResiC.(uv{i_uv})]=estimationSimulationIntraConditional(uv{i_uv}, g, radar, idt_y, icov, condVal, Resi.(uv{i_uv})); % ~85 min
        % clear Resi

        f_trans_inver = @(X) interp1(transi.([uv{i_uv} '_axis_t']),transi.([uv{i_uv} '_axis']),( X ... 
        .* victi_trend.(['t_std_' uv{i_uv}])(g.time_doy(idt_y))' ...
        .* victi_trend.(['sf_std_' uv{i_uv}])(double([g.LON(~g.mask_water) g.LAT(~g.mask_water)])))...
        ,'pchip');
    
        % Estimation & std
        % transformation
        Eim.(uv{i_uv}) = f_trans_inver(Ei.(uv{i_uv}));
        Sim.(uv{i_uv}) = [];%f_trans_inver(SiU);
        % add daily
        Em.(uv{i_uv}) = Eim.(uv{i_uv}) + Edm.(uv{i_uv})(:,g.day_id(idt_y));% + ;
        
        % Sm = Sdm(g.day_id,:) + Sim;
        % clear Em Sm Eim Sim

        % Simulation
        for k=1:sim_n
            ResiCm = f_trans_inver(ResiC.(uv{i_uv})(:,:,k));
            ResCm.(uv{i_uv}) = ResiC.(uv{i_uv})(:,:,k) + ResdCm.(uv{i_uv})(:,g.day_id(idt_y),k);
        end
    end

    % Save
    save(['data/speed/est_uv_' num2str(i_y)],'Em','idt_y','-v7.3') 

    for k=1:sim_n
            save(['data/speed/sim_uv_' num2str(i_y) '_' num2str(k)],'ResCm','idt_y','-v7.3')
        end
        
end








%% 
i_y=2021;


%% Illustration

vi.u =  vici.u + vicd.u(g.day_id,:);% +
vi.v = vici.v + vicd.v(g.day_id,:);%+ 

i_s=6500;
i_step = 4;
i_uv=1;
figure('position',[0 0 1600 900]); ha=tight_subplot(3,3);
for u=0:8
    axes(ha(u+1))
    hold on; set(gca,'ydir','normal');
    tmp=nan(size(g.mask_water));
    tmp(~g.mask_water) = ResCm.(uv{i_uv})(:,i_s+u*i_step,1); % EmbT(:,i_s+u,k);% 
    imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    borders('states','k');
    id=~isnan(vi.(uv{i_uv})(idt_y(i_s)+u*i_step,:)); % vidti  vidtditv
    scatter(radar.lon(id),radar.lat(id),100,vi.(uv{i_uv})(idt_y(i_s)+u*i_step,id)','filled','MarkerEdgeColor','k');
    title(datestr(g.time(idt_y(i_s)+u*i_step)))
    axis([-125 -68 23 50]); colorbar;  % caxis([-5 11])
end


%%

[~,id_cond] = min( sqrt( (radar.lon'-g.LON(~g.mask_water)).^2 + (radar.lat'-g.LAT(~g.mask_water)).^2 ) );

tmp = ResCm.(uv{i_uv})(:,:,1);

tmp_hd = tmp(id_cond,:)';
tmp_hd(isnan(victidt.(uv{i_uv})(idt_y,:)))=nan;

figure; histogram(tmp_hd,'Normalization','pdf')
hold on; histogram(vi.(uv{i_uv})(idt_y,:),'Normalization','pdf')

figure; histogram(tmp_hd-victidt.(uv{i_uv})(idt_y,:))




%%

sfv = Ei.v;
ptv = victidt.v;

%find(any(~isnan(Em)))
i_s=6550;%10050;%numel(i_yy)-15;
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
figure; hold on; set(gca,'ydir','normal'); axis equal tight
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = Ei(:,i_ey(1));
im=imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
borders('states','k');
s=scatter(radar.lon,radar.lat,100,tmpE(:,i_s),'filled','MarkerEdgeColor','k');
axis([-125 -68 23 50]); colorbar; caxis([-3 3])

for i_d=1:numel(unique(day_id))
    i_en = find(g.day_id==day_id_unique(i_d));
    i_en(all(isnan(NNT(i_en,:)),2))=[];
    i_ey = find(ismember(idt_y,i_en));
    
    tmpE = vidtditv(i_en,:)';
    tmpEisnan = find(~isnan(tmpE(:)));
        
    for i_s=1:numel(i_ey)
        tmp(~g.mask_water) = Ei(:,i_ey(i_s));
        title(datestr(g.time(idt_y(i_ey(i_s)))))
        im.CData = tmp;
        delete(s)
        tmp2 = ~isnan(tmpE(:,i_s));
        s=scatter(radar.lon(tmp2),radar.lat(tmp2),100,tmpE(tmp2,i_s),'filled','MarkerEdgeColor','k');
        refreshdata
        drawnow
        % pause(0.2);
           [imind,cm] = frame2im(getframe(gcf)); 
  % Write to the GIF File 
  if i_d == 1 
      imwrite(imind,cm,'data/estimation_intra','gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,'data/estimation_intra','gif','WriteMode','append'); 
  end
    end
end





