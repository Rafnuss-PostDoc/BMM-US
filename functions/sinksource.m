function [Fd, Ts, gext, MVT_day]=sinksource(g,i_y,simest,saveit)

plotit=false;
%%
% *Build extended grid*
%
% Build a grid |gext| with a cell more on each side of the domain (ext=extanded)
% in order to compute the fluxes in/out of the existing grid

gext.lat = fillmissing([nan(1,1) ;g.lat' ;nan(1,1)],'linear');
gext.lon = fillmissing([nan(1,1); g.lon'; nan(1,1)],'linear');
gext.nlat = numel(gext.lat);
gext.nlon = numel(gext.lon);
[gext.LAT, gext.LON] = ndgrid(gext.lat,gext.lon);
gext.mask_water = padarray(g.mask_water,[1 1],true);
gext.mask_out=movmean(movmean(~gext.mask_water,[1 1],2),[1 1],1)>0;
gext.mask_out(~gext.mask_water)=false;

%%
% *Build extended time grid*
%
% First, we select only the data during the year i_y and the night data.
% This should already be computed and loaded with the estimation.
load(['data/density/est_' num2str(i_y) '.mat'],'idt_y')

% The extended grid also include the time between two time step of the previous
% grid. Yet, because there were no time step between the last time step of the
% end of night (sunrise) and the first of the next day (sunset), we add one at
% noon every day.
% For the US data, we set noon at 1800 (anytime between the gap on the plot
% below:
% plot(unique(gext.time-dateshift(g.time,'start','day')))
gext.time_s = g.time(idt_y);
g_dmd = g.day(unique(g.day_id(idt_y)))+18/24;
gext.time = sort( [gext.time_s; g_dmd]);
gext.time_s_id = ismember(gext.time,gext.time_s);

gext.nt = numel(gext.time);
gext.day_id=nan(gext.nt,1);
gext.day_id(gext.time_s_id) = g.day_id(idt_y);
gext.day_id = fillmissing(gext.day_id,'next');
gext.day_id(end)=gext.day_id(end-1);
gext_day_id_unique = unique(gext.day_id);
gext.day = g.day(gext_day_id_unique);
gext.nat = numel(gext.day);

% The time resolution is 15 min -> 15/60 hr
dt=hours(unique(diff(g.time)));

% Prealocate variable
vy = nan(g.nlat, g.nlon, gext.nt); vx=vy;
rho = nan(g.nlat, g.nlon, gext.nt);
% mask of water
g_mw = repmat(~g.mask_water,1,1,gext.nt); g_mw(:,:,~gext.time_s_id)=false;



if simest=="est"
    % *Compute flight speed and density*
    load(['data/density/est_' num2str(i_y) '.mat'])
    rho(g_mw) = EmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
    clear EmbT
    load(['data/speed/est_uv_' num2str(i_y) '.mat'])
    vx(g_mw) = Em.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
    vy(g_mw) = Em.v/1000*60*60; % m/s -> km/h (+) north, (-) south
    clear Em
    K=1;
else
    dirsim = dir(['data/density/sim_' num2str(i_y) '_*']);
    K=numel(dirsim);
end


Fd.takingoff=nan(g.nlat,g.nlon,gext.nat,K);
Fd.landing=nan(g.nlat,g.nlon,gext.nat,K);
% Fd.landing2=nan(g.nlat,g.nlon,gext.nat);
Fd.entering=nan(gext.nlat,gext.nlon,gext.nat,K);
Fd.leaving=nan(gext.nlat,gext.nlon,gext.nat,K);
% F_day=nan(gext.nlat,gext.nlon,gext.nat);

MVT_day = nan(g.nlat,g.nlon,gext.nat,K);



tic
for k=1:K

    if simest=="sim"
        load(['data/density/sim_' num2str(i_y) '_' num2str(k) '.mat'])
        rho(g_mw) = ResCmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
        clear ResCmbT
        load(['data/speed/sim_u_' num2str(i_y) '_' num2str(k) '.mat'])
        vx(g_mw) = ResCm.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
        load(['data/speed/sim_v_' num2str(i_y) '_' num2str(k) '.mat'])
        vy(g_mw) = ResCm.v/1000*60*60; % m/s -> km/h (+) north, (-) south
        clear ResCm
    end

    % *Compute the flux* $\Phi =\mathit{\mathbf{v}}\rho$ *at +/- 1/2 grid cell.*
    % First, add a nan layer in lat or lon direction
    Phiy_pad = padarray(rho .* vy .* repmat(g.dx,1,g.nlon,gext.nt) ,[1 0 0],nan); % bird/km^2 * km/h *km ->  bird/h
    Phix_pad = padarray(rho .* vx .* g.dy,[0 1 0],nan);

    % and then, compute the flux at +/- 1/2 even if either previous or next cells
    % is nan.
    Phiy_h = movmean(Phiy_pad,[0 1],1,'omitnan','Endpoints','discard');
    Phix_h = movmean(Phix_pad,[0 1],2,'omitnan','Endpoints','discard');
    % clear Phix_pad Phiy_pad

    % *Compute the delta flux / delta distance* $\frac{\Delta \Phi }{\Delta \left(\textrm{lat},\textrm{lon}\right)}$
    % First, add 0 padding for the outer zone (allows to compute boundary cell)
    Phiy_h_0=padarray(Phiy_h,[1 1 0],0);
    Phix_h_0=padarray(Phix_h,[1 1 0],0);
    % clear Phiy_h Phix_h

    % Then, replace the nan by zero to compute the change of fluxes even no fluxes
    % on the other side (i.e. boundary cell)
    Phiy_h_0(isnan(Phiy_h_0))=0;
    Phix_h_0(isnan(Phix_h_0))=0;

    % Finally, compute the delta flux over delta distance.
    dPhiy = diff(Phiy_h_0,1,1);
    dPhix = diff(Phix_h_0,1,2);
    % clear Phix_h_0 Phiy_h_0

    % *Compute the variation of bird in/out each cell*
    F = (dPhiy + dPhix ).*dt; % bird/h * hr -> bird
    % clear dPhiy dPhix

    % Crank-nicolson Move F at t+1/2:
    % F(F==0)=nan;
    % F = movmean( F ,[0 1],3,'omitnan');
    % F(isnan(F))=0;

    % F comprise both the intral change of bird and the bird going out of the area
    % Inner flux: remove the boundary cells
    Fin = F;
    %Fin(~gext.mask_water)=0;
    Fin(repmat(gext.mask_water,1,1,gext.nt)) = 0;

    % Outer Flux. Note that Philat_h_0 and Philat_h_0 are 0 for the cell outside
    % the domain, such at dPhilatdlat is equal to the flux at the outer cell (with
    % the sign corresponding to the direction).
    Fout = F;
    %Fout.Fout(gext.mask_water)=0;
    Fout(repmat(~gext.mask_water,1,1,gext.nt)) = 0;
    % clear F

    % *Compute the landing/derpature* $W$
    % |rho0| is used (NaN replaced by 0) to account for landing as sunset and sunrise.
    % Unit of W is [bird =bird (*15min/15min)], when (+) Birds entering the airspace
    % (takingoff) and (-) Birds leaving the raispace (landing)
    rho0 = rho; rho0(isnan(rho)) = 0; % bird/km^2
    W = diff(rho0 .* repmat(g.area,1,1,gext.nt),1,3) + Fin(2:end-1,2:end-1,1:end-1); % bird/km^2 * km^2  + bird -> bird
    % clear rho0 Fin

    % *Check methodology*
    %
    % Check for a step
    if plotit
    i=20; i=70; i=56;% i=1095;
    i=2605;i=10000;
    gext.time(i);
    figure('position',[0 0 1200 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc(rho(:,:,i),'AlphaData',~isnan(rho(:,:,i))); set(gca,'ydir','normal'); colorbar; title('\rho')
    nexttile; imagesc(vy(:,:,i),'AlphaData',~isnan(vy(:,:,i))); set(gca,'ydir','normal'); colorbar; title('vy')
    nexttile; imagesc(vx(:,:,i),'AlphaData',~isnan(vx(:,:,i))); set(gca,'ydir','normal'); colorbar; title('vx')
    figure('position',[0 0 1200 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc(Phix_pad(:,:,i),'AlphaData',~isnan(Phix_pad(:,:,i))); set(gca,'ydir','normal'); colorbar; title('\Phi_{x}')
    nexttile; imagesc(Phiy_pad(:,:,i),'AlphaData',~isnan(Phiy_pad(:,:,i))); set(gca,'ydir','normal'); colorbar; title('\Phi_{y}')
    nexttile; imagesc(dPhix(:,:,i),'AlphaData',~(0==dPhix(:,:,i))); set(gca,'ydir','normal'); colorbar; title('\Delta\Phi_{x}/\Delta x')
    nexttile; imagesc(dPhiy(:,:,i),'AlphaData',~(0==dPhiy(:,:,i))); set(gca,'ydir','normal'); colorbar; title('\Delta\Phi_{y}/\Delta y')
    nexttile; imagesc(F(:,:,i),'AlphaData',~(0==F(:,:,i))); set(gca,'ydir','normal'); colorbar; title('F=(\Delta\Phi_{y}+\Delta\Phi_{x})\Delta t')

    % Check the spatial aggregated mass balence.
%     a = reshape(sum(sum(Fout,1),2),1,[]);
%     b = reshape(sum(sum(rho0 .* repmat(g.area,1,1,gext.nt),1),2),1,[]);
%     c = reshape(sum(sum(W,1),2),1,[]);
%     figure('position',[0 0 1000 200]); hold on; xlabel('time');ylabel('Error')
%     plot(b(2:end) - b(1:end-1) - c - a(1:end-1))

    % See a single timestep.
    figure('position',[0 0 1300 600]);tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc( rho0(:,:,i),'AlphaData',~(0==rho0(:,:,i))); title('\rho(t)');set(gca,'ydir','normal'); colorbar;
    nexttile; imagesc( rho0(:,:,i+1),'AlphaData',~(0==rho0(:,:,i+1))); title('\rho(t+1)'); set(gca,'ydir','normal'); colorbar;
    nexttile; imagesc( rho0(:,:,i+1)-rho0(:,:,i),'AlphaData',~(0==rho0(:,:,i))); title('\rho(t+1) - \rho(t)'); set(gca,'ydir','normal'); colorbar;
    nexttile; imagesc(Fin(:,:,i),'AlphaData',~(0==Fin(:,:,i)));  title('F_{in}'); set(gca,'ydir','normal'); colorbar;
    nexttile; imagesc(W(:,:,i),'AlphaData',~(0==W(:,:,i)));  title('W^{t\rightarrow t+1}=\rho(t+1) - \rho(t)+F_{in}'); set(gca,'ydir','normal'); colorbar;
    nexttile; imagesc(Fout(:,:,i),'AlphaData',~(0==Fout(:,:,i)));  title('F_{out}'); set(gca,'ydir','normal'); colorbar;
    end
    %% *Migratory Processes*
    % *Seperature the landing/landing (-) and takingoff/take-off (+)*
    landing=W; 
    landing(landing>=0)=0;
    % landing2 
%     landing2 = landing;
%     landing2(gext.mask_water(2:end-1,2:end-1,1:end-1))=0;
    takingoff=W;
    takingoff(takingoff<=0)=0;
    leaving=Fout; 
    leaving(leaving>=0)=0;
    entering=Fout;
    entering(entering<=0)=0;

    for i_t=1:gext.nat-1
        idt=gext.day_id(1:end-1)==gext_day_id_unique(i_t);
        Fd.takingoff(:,:,i_t,k)=sum(takingoff(:,:,idt),3,'omitnan');
        Fd.landing(:,:,i_t,k)=sum(landing(:,:,idt),3,'omitnan');
        % Fd.landing2(:,:,i_t)=sum(landing2(:,:,idt),3);
        Fd.entering(:,:,i_t,k)=sum(entering(:,:,idt),3,'omitnan');
        Fd.leaving(:,:,i_t,k)=sum(leaving(:,:,idt),3,'omitnan');
        % F_day(:,:,i_t)=sum(F(:,:,idt),3,'omitnan');
    end
    Fd.takingoff(Fd.takingoff==0)=nan;
    Fd.landing(Fd.landing==0)=nan;
    % Fd.landing2(Fd.landing2==0)=nan;
    % Fd.W=Fd.takingoff+Fd.landing;

    % *Compute timeseries*
    Ts.W(:,k) = reshape(sum(sum(W,1,'omitnan'),2,'omitnan'),1,[]);
    Ts.landing(:,k) = reshape(sum(sum(landing,1,'omitnan'),2,'omitnan'),1,[]);
    Ts.takingoff(:,k) = reshape(sum(sum(takingoff,1,'omitnan'),2,'omitnan'),1,[]);
    Ts.entering(:,k) = reshape(sum(sum(entering,1,'omitnan'),2,'omitnan'),1,[]);
    Ts.leaving(:,k) = reshape(sum(sum(leaving,1,'omitnan'),2,'omitnan'),1,[]);
    Ts.landing_day(:,k) = reshape(sum(sum(Fd.landing(:,:,:,k),1,'omitnan'),2,'omitnan'),1,[]);
    Ts.takingoff_day(:,k) = reshape(sum(sum(Fd.takingoff(:,:,:,k),1,'omitnan'),2,'omitnan'),1,[]);
    Ts.entering_day(:,k) = reshape(sum(sum(Fd.entering(:,:,:,k),1,'omitnan'),2,'omitnan'),1,[]);
    Ts.leaving_day(:,k) = reshape(sum(sum(Fd.leaving(:,:,:,k),1,'omitnan'),2,'omitnan'),1,[]);

    % MTR
    MVT = rho .* sqrt(vy.^2 + vx.^2) .* dt; % bird/km^2 .* km/hr .* hr -> bird/km
    for i_t=1:gext.nat-1
        idt=gext.day_id(1:end-1)==gext_day_id_unique(i_t);
        MVT_day(:,:,i_t,k) = sum(MVT(:,:,idt),3,'omitnan');
    end
    
    ttoc = toc;
    disp(['sink/source - ' num2str(i_y) ' - ' simest ' - ' num2str(k) '-' num2str(ttoc/60) ' min'])  
end

%  *Save*
if saveit
    if simest=="est"
        save(['data/flow/est_' num2str(i_y) '.mat'],'Fd','Ts','gext','MVT_day','-v7.3');%,'vy','vx','rho');
        % save(['data/flow/est_' num2str(i_y) '_full.mat'],'rho','W','gext','MVT_day','-v7.3')
    else
        save(['data/flow/sim_' num2str(i_y) '.mat'],'Fd','gext','MVT_day','-v7.3');%,'vy','vx','rho');
    end
end


end