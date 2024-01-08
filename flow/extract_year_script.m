
load('../data/density/inference-trans.mat')
%load('../data/speed/inference-trans.mat')
addpath('../functions/')
%vid = trans.f_inv(vidTS);

for i_y=1995:2021

    % Load flow data
    load(['../data/flow/est_' num2str(i_y) '.mat'],'gext')
    
    % load intra-night data
    load(['../data/density/est_' num2str(i_y) '.mat'],'idt_y')
    g_mw = repmat(~g.mask_water,1,1,numel(idt_y));
    vy = nan(numel(g.lat), numel(g.lon), numel(idt_y)); vx=vy; rho = vy;
    
    load(['../data/density/est_' num2str(i_y) '.mat'])
    rho(g_mw) = EmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
    load(['../data/speed/est_uv_' num2str(i_y) '.mat'])
    vx(g_mw) = Em.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
    vy(g_mw) = Em.v/1000*60*60; % m/s -> km/h (+) north, (-) south
    
    gy.time = g.time(idt_y);
    gy.day_id = g.day_id(idt_y)-min(g.day_id(idt_y))+1;
    gy.day = g.day(unique(g.day_id(idt_y)));

    rhod=nan(numel(g.lat), numel(g.lon), numel(gy.day));
    vxd=rhod;
    vyd=rhod;
    for i = 1:numel(gy.day)
        id = gy.day_id==i;
        rhod(:,:,i) = mean(rho(:,:,id),3,'omitnan');
        vxd(:,:,i) = mean(vx(:,:,id),3,'omitnan');
        vyd(:,:,i) = mean(vy(:,:,id),3,'omitnan');
    end

    save("../data/flow/rho-velocity-"+i_y+".mat","rhod","vyd","vxd","-v7.3")

end
