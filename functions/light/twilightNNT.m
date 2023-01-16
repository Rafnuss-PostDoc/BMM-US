function [NNT,twlSet,twlRise,date] = twilightNNT(tm, lon, lat, varargin) 
%TWILIGHT estimate time of sunrise or sunset for a given day and location
%   Twilight uses an iterative algorithm to estimate times of sunrise and
%   sunset. Note that these functions return the twilight that occurs on
%   the same date GMT as tm, and so sunset may occur before sunrise,
%   depending upon latitude. Solar declination and equation of time vary
%   slowly over the day, and so the values of the Solar declination and
%   equation of time at sunrise/sunset are well approximated by their
%   values at 6AM/6PM local time. The sun's hour angle and hence
%   sunrise/sunset for the required zenith can then be caclulates from
%   these approximations. The calculation is then repeated using the
%   approximate sunrise/sunset times to derive more accurate values of the
%   Solar declination and equation of time and hence better approximations
%   of sunrise/sunset. The process is repreated and is accurate to less
%   than 2 seconds within 2 or 3 iterations.

if nargin==3
    removeDay=true;
else
    removeDay =varargin{1};
end

% The angular diameter of the sun is about 0.536 degrees, therefore the moment of sunrise/sunset corresponds to half that elevation at -0.268 degrees.
zenith = 90; %90+0.268; % 96;
tol=3; % change of max 3 min
cond = true;

tm=tm(:);
lat=lat(:)';
lon=lon(:)';
lon = mod(lon + 180,360) - 180;

date = ((dateshift(min(tm),"start",'day')-2):(dateshift(max(tm),"start",'day')+2))';
% date = dateshift
rise = [true(size(date));false(size(date))];
date2 = [date;date];

twl = date2 + 240 .* (90+180*~rise - lon)/60/60/24;

while cond
    s = solar(twl);
    s.solarTime = mod(s.solarTime,360);
    solarTime = 4 * twilightSolartime(s, lon, lat, rise, zenith) - s.eqnTime;
    tmp = 60 * solarTime / 60/60/24;
    tmp(tmp<.5&~rise)=tmp(tmp<.5&~rise)+1;
    twlN = date2 + tmp;
    if max(min(abs([minutes(twl(:)-twlN(:)) minutes(twl(:)-twlN(:)+1) minutes(twl(:)-twlN(:)-1)]),[],2))<tol
        cond=false;
    end
    twl = twlN;
end

twlRise = twl(1:end/2,:);
twlSet = twl(end/2+1:end,:);

NNT=nan(numel(tm),numel(lat));

for i_r=1:numel(lat)
    SetPrev = interp1(twlSet(:,i_r),1:numel(date),tm,'previous');
    SetNext = interp1(twlSet(:,i_r),1:numel(date),tm,'next');
    % RisePrev = interp1(twlRise(:,i_r),1:numel(date),tm,'previous');
    RiseNext = interp1(twlRise(:,i_r),1:numel(date),tm,'next');
    
    tmp = (tm - mean([twlRise(RiseNext,i_r)  twlSet(SetPrev,i_r)],2)) ./ (twlRise(RiseNext,i_r) - twlSet(SetPrev,i_r)) *2;
    if removeDay
        isDay = twlSet(SetNext,i_r)<twlRise(RiseNext,i_r);
        tmp(isDay) = nan;
    end
        NNT(:,i_r) = tmp;
%     figure; hold on;
%     %imagesc(datenum(tm),[0,1],isDay')
%     plot(datenum(tm(1:1000)),NNT(1:1000,i_r),'-k')
%     xline(datenum(twlSet(1:1000,i_r)),'g')
%     xline(datenum(twlRise(1:1000,i_r)),'r')

end


end