function [E,S,ResC]=estimationSimulationIntraConditional(duv, g, radar, idt_y, icov, condVal, Res)

% Estimation is perform per block of gn consecutive days adding the +/- gm days in the conditioning values. 
% Define groupping duration
% gn % number of date to simulation togeter
% gnm % margin around those date


%% Setup
% Number of grid element
g_nlm = sum(~g.mask_water(:));
% nr = height(radar);

% Pre-build distance matrix
Ddist_rr=squareform(pdist([radar.lat radar.lon],@lldistkm));
Ddist_rg = pdist2([radar.lat radar.lon], double([g.LAT(~g.mask_water) g.LON(~g.mask_water)]), @lldistkm);

if ~isempty(Res)
    % get conditioning point
    [~,id_cond] = min( sqrt( (radar.lon'-g.LON(~g.mask_water)).^2 + (radar.lat'-g.LAT(~g.mask_water)).^2 ) );
    Res_cond = Res(id_cond,:,:);
    ResC = nan(size(Res),'single');
    assert(size(Res,2)==numel(idt_y),'issue with the time axes between the unconditional simulation Res and the current estmation grid')
else
    ResC=[];
end

% Estimated and varaince value
E = nan(g_nlm,numel(idt_y));
S = nan(g_nlm,numel(idt_y));

% find the day of year
day_id = g.day_id(idt_y);
day_id_unique = unique(day_id);
    
tic 
disp("warning")
for i_d=100:100%numel(day_id_unique)
    
    % Find the index of the element to sim/est for the year (i.e., in idt_y)
    i_ey = find(day_id==day_id_unique(i_d));

    % Find the index of the element to estimate on the full time (absolute id)
    i_en = idt_y(i_ey);
    
    % compute their delta time
    idt = datenum(g.time(i_en)-g.day(g.day_id(i_en)))*24*4;
    
    % find conditioning values
    tmpE = condVal(i_en,:)';
    tmpEisnotnan = find(~isnan(tmpE(:)));
    [r_s, r_t]=ind2sub(size(tmpE), tmpEisnotnan);

    % Define the parameter of the covariance
    if duv=='d'
        parm = icov.f_parm(g.day_doy(day_id(day_id_unique(i_d))));
    elseif duv=='u'
        parm = icov.f_parm_u(g.day_doy(day_id(day_id_unique(i_d))));
    else
        parm = icov.f_parm_v(g.day_doy(day_id(day_id_unique(i_d))));
    end

    Cdist_rr = icov.Gneiting_dist(Ddist_rr(r_s,r_s),parm(3),parm(6));
    Ctime_rr = icov.Gneiting_time(abs(r_t-r_t')/4/24,parm(4),parm(5));  
    Crr = parm(2) * Cdist_rr .* Ctime_rr;
    Crr = Crr + parm(1) * eye(height(Crr));

    Cdist_rg = icov.Gneiting_dist(Ddist_rg,parm(3),parm(6));

    % Compute the Kriging Weight
    %tic; Lambda = inv(Crr)*Crg(tmpEisnan,:); toc

    for i_t=1:numel(i_en)
        
        Dtime_rg=abs(idt(r_t)-idt(i_t))/4/24;
        
        % limit to delta time below the range time paramter
        id_dr = Dtime_rg < parm(4);
        % further limit to 4*145 pts for computational reason
        while sum(id_dr)>145*4
            dt_unique = sort(unique(Dtime_rg(id_dr)));
            id_dr = Dtime_rg <= dt_unique(end-1);
        end
        Ctime_rg = icov.Gneiting_time(Dtime_rg(id_dr),parm(4),parm(5));
        Crg = parm(2) * Cdist_rg(r_s(id_dr),:) .* repelem(Ctime_rg,1,g_nlm);
        
        % Solve Kriging system
        Lambda = Crr(id_dr,id_dr) \ Crg;
        
        % Kriging estimation
        E(:,i_ey(i_t)) = Lambda' * tmpE(tmpEisnotnan(id_dr));
        
        % Kriging Variance
        S(:,i_ey(i_t)) = sqrt(parm(2) - sum(Lambda.*Crg));
 
        % Conditional simulation
        if ~isempty(Res)
            for k=1:size(Res_cond,3)
                tmpRc = Res_cond(:,i_ey,k);
                ResC(:,i_ey(i_t),:) = Res(:,i_ey(i_t),:) + Lambda' * (tmpE(tmpEisnotnan(id_dr))-tmpRc(tmpEisnotnan(id_dr))); 
            end
        end
    end
    
    ttoc = toc;
    disp(['est/sim - ' datestr(g.day(day_id_unique(i_d))) ' - ' num2str(ttoc/60,3) ' min'])   
end

end