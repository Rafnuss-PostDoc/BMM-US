function [E,S,ResC]=estimationSimulationDailyConditional(duv, g, gn, gnm, radar, dcov, condVal, Res)

% Estimation is perform per block of gn consecutive days adding the +/- gm days in the conditioning values. 
% Define groupping duration
% gn % number of date to simulation togeter
% gnm % margin around those date


%% Setup Estiamtion grid
% Number of grid element
g_nlm = sum(~g.mask_water(:));
nr = height(radar);

% Pre-build distance matrix
Ddist_rr=squareform(pdist([radar.lat radar.lon],@lldistkm));
Ddist_rg = pdist2([radar.lat radar.lon], double([g.LAT(~g.mask_water) g.LON(~g.mask_water)]), @lldistkm);

%% Prepare values

% Estimated and varaince value
E = nan(g_nlm,numel(g.day),'single');
S = nan(g_nlm,numel(g.day),'single');

%% Simulated value
if ~isempty(Res)
    % get conditioning point
    [~,id_cond] = min( sqrt( (radar.lon'-g.LON(~g.mask_water)).^2 + (radar.lat'-g.LAT(~g.mask_water)).^2 ) );
    Res_cond = Res(id_cond,:,:);
    ResC = nan(size(Res),'single');
else
    ResC=[];
end


% Covariance matrix is the same from one year to the next but change over the course of a year. So, to save computational coast, we estimate for the
% same period of all the years. 
% First we loop over the day of year, and then through all years.
tic 
for i_eg=1:gn:366
    
    % estimated index
    i_g=i_eg+(1:gn)'-1;
    i_g=i_g(i_g<366);
    
    % conditioning index
    i_n = ((i_g(1)-gnm):(i_g(end)+gnm))';
    % number of conditioning
    nn = numel(i_n);

    if duv=='d'
        % Define the parameter of the covariance
        parm = mean(dcov.f_parm(g.day_doy(i_g)'),2);

        % Compute the cross-covriance radar-radar
        Cdist_rr = parm(2)*(Ddist_rr==0)+dcov.Gneiting_dist(Ddist_rr,parm(4),parm(7));
        Ctime_rr = dcov.Gneiting_time(squareform(pdist(i_n)),parm(5),parm(6));
        Crr = parm(3) * repmat(Cdist_rr,nn,nn) .* repelem(Ctime_rr,nr,nr);
        Crr = Crr + parm(1) * eye(height(Crr)) + parm(2)*repmat(Cdist_rr==0,nn,nn);
        
        % Cross-covariance of the radar data
        Cdist_rg = parm(2)*(Ddist_rg==0)+dcov.Gneiting_dist(Ddist_rg,parm(4),parm(7));
        Ctime_rg = dcov.Gneiting_time(pdist2(i_n,i_g),parm(5),parm(6));
        % Build the full covariance matrix
        Crg = parm(3) * repmat(Cdist_rg,nn,numel(i_g)) .* repelem(Ctime_rg,nr,g_nlm);
    else
        if duv=='u'
            parm = mean(dcov.f_parm_u(g.day_doy(i_g)'),2);
        else
            parm = mean(dcov.f_parm_v(g.day_doy(i_g)'),2);
        end

        % Compute the cross-covriance radar-radar
        Cdist_rr = dcov.Gneiting_dist(Ddist_rr,parm(3),parm(6));
        Ctime_rr = dcov.Gneiting_time(squareform(pdist(i_n)),parm(4),parm(5));
        Crr = parm(2) * repmat(Cdist_rr,nn,nn) .* repelem(Ctime_rr,nr,nr);
        Crr = Crr + parm(1) * eye(height(Crr));
        
        % Cross-covariance of the radar data
        Cdist_rg = dcov.Gneiting_dist(Ddist_rg,parm(3),parm(6));
        Ctime_rg = dcov.Gneiting_time(pdist2(i_n,i_g),parm(4),parm(5));
        % Build the full covariance matrix
        Crg = parm(2) * repmat(Cdist_rg,nn,numel(i_g)) .* repelem(Ctime_rg,nr,g_nlm);
    end

    disp('warning')
   for i_y=2010:2010
        % find the index of the first day of the year to estimate
        i_yy = find(year(g.day)==i_y,1);
        
        % Find the index of the position to estimate
        i_s = i_yy+i_g-1;
        
        % find possible neighborhood (index of conditioning data)
        i_e=i_yy+i_n-1;
        i_eid = i_e>0 & i_e<numel(g.day);
        i_e = i_e(i_eid);
        
        % find conditioning values
        tmpE = nan(nr,numel(i_n));
        tmpE(:,i_eid) = condVal(i_e,:)';
        tmpEisnotnan = ~isnan(tmpE(:));

        % Compute the Kriging Weight
        %tic; Lambda = inv(Crr(tmpEisnan,tmpEisnan))*Crg(tmpEisnan,:); toc
        Lambda = Crr(tmpEisnotnan,tmpEisnotnan)\Crg(tmpEisnotnan,:);
        
        % Kriging estimation
        E(:,i_s) = reshape(Lambda' * tmpE(tmpEisnotnan),[],numel(i_g));
        % Kriging Variance
        if duv=='d'
            %S(:,i_s) = reshape( sqrt(parm(3) - sum(Lambda.*Crg(tmpEisnotnan,:))),[],numel(i_g));
        else
            %S(:,i_s) = reshape( sqrt(parm(2) - sum(Lambda.*Crg(tmpEisnotnan,:))),[],numel(i_g));
        end
        
        % Conditional simulation
        if ~isempty(Res)
            for k=1:size(Res_cond,3)
                tmpRc = nan(nr,numel(i_n));
                tmpRc(:,i_eid) = Res_cond(:,i_e,k);
                ResC(:,i_s,k) = Res(:,i_s,k) + reshape(Lambda' * (tmpE(tmpEisnotnan)-tmpRc(tmpEisnotnan)),[],numel(i_g));  
            end
        end
        
        tmp = toc;
        disp([num2str(i_eg)  '-' num2str(i_y) '-' num2str(tmp/60)])
    end
end

end