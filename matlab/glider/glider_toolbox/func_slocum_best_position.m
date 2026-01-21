function [data] = func_slocum_best_position(data,op);
% function [data] = func_slocum_best_position(data,op);
% 
% GEOMAR SVN $Id: func_slocum_best_position.m 950 2022-10-14 15:42:13Z gkrahmann@geomar.de $
%
% best estimate of glider position from all available data
%
% input  :  data                    - data structure with the following fields
%           .lat                      glider's m_lat (dead-reckoning & GPS)
%           .lon                      glider's m_lon (dead-reckoning & GPS)
%           .gps_lat                  glider's m_gps_lat (GPS positions)
%           .gps_lon                  glider's m_gps_lon (GPS positions)
%           .inv_lat                  glider's m_gps_invalid_lat
%           .inv_lon                  glider's m_gps_invalid_lon
%           .far_lat                  glider's m_gps_toofar_lat
%           .far_lon                  glider's m_gps_toofar_lon
%           .tim                      glider's m_present_time
%           .dep                      glider's m_depth
%           .ign_lat                  glider's m_gps_ignored_lat
%           .ign_lon                  glider's m_gps_ignored_lon
%           .dos_id                   glider's dos formatted bd-file id
%           use_median       [1]   - 1:use median for gps surface, else extrapol.
%           min_dive_length  [250] - minimum dive length (in whatever units, 250=10 min)
%           display_spurious [0]   - 1:display all problem cases  0:do not display just fix them
%
% output :  nlat        - best estimate of latitude
%           nlon        - best estimate of longitude
%
% version 12.1.0  last change 25.01.2023

% G.Krahmann, GEOMAR, Sep 2006

% modified from older routine                                          GK, 23.08.2012  0.5-->1
% check indices                                                        GK, 12.03.2013  1-->2
% add check for too short gooddr and more out of bounds                GK, 09.10.2013  2-->3
% no-plot option                                                       GK, 27.05.2014  3-->4
% cleanup of disp values                                               GK, 24.02.2017  4-->5
% more text for problem                                                GK, 06.09.2017  5-->6
% fixed the case of thrown-out good fixes within the dive when the glider does
% very shallow upper returns.                                          GK, 18.09.2017  6-->7
% allow for lat_range and lon_range options                            GK, 31.01.2018  7-->8
% catch cases when at the same time stamp an invalid GPS and a good GPS position
% is stored. In these cases we throw both away.                        GK, 19.02.2018  8-->9
% better discriminate short/shallow dead reckoning intervals           GK, 21.02.2018  9-->10
% save final plot, improve short dive plot                             GK, 15.01.2019  10-->11
% improvements                                                         GK, 18.12.2020  11-->12.0.0
% catch gps_overrun setting                                            GK, 04.01.2021  12.0.0-->12.0.1
% catch ignored when last value                                        GK, 22.09.2022  12.0.1-->12.0.2
% changed function call name                                           GK, 25.01.2023  12.0.2-->12.1.0


%
% catch all Webb dummy values in lat/lon 
%
fn = fieldnames(data);
for n=1:length(fn)
  if ~isempty(findstr(fn{n},'_lat'))
    dummy = getfield(data,fn{n});
    bad = find(dummy==69696969);
    if ~isempty(bad)
      dummy(bad) = nan;
      data = setfield(data,fn{n},dummy);
    end
  end
  if ~isempty(findstr(fn{n},'_lon'))
    dummy = getfield(data,fn{n});
    bad = find(dummy==69696969);
    if ~isempty(bad)
      dummy(bad) = nan;
      data = setfield(data,fn{n},dummy);
    end
  end
end


%
% set latitude and longitude ranges to check
%
if ~isfield(op,'lat_range')
  lat_range = [-90,90];
else
  lat_range = op.lat_range;
end
if ~isfield(op,'lon_range')
  lon_range = [-361,361];
else
  lon_range = op.lon_range;
end


%
% fill arrays, as they might be stored as sparse vectors
% thereafter convert to decimal degrees from Webb's odd convention
%
tim = data.m_present_time;
dep = data.m_depth;
[lat,lon] = parse_structures(data.m_lat,data.m_lon);
orig_lat = lat;
orig_lon = lon;
[gps_lat,gps_lon] = parse_structures(data.m_gps_lat,data.m_gps_lon);
[inv_lat,inv_lon] = parse_structures(data.m_gps_invalid_lat,data.m_gps_invalid_lon);
[far_lat,far_lon] = parse_structures(data.m_gps_toofar_lat,data.m_gps_toofar_lon);
[ign_lat,ign_lon] = parse_structures(data.m_gps_ignored_lat,data.m_gps_ignored_lon);

%figure(6)
%clf
%plot(tim,lat,'xb')
%hold on
%plot(tim,gps_lat,'r+')
%ylim([14.5,15])
%pause

%

%
% combine gps utc information to get the gps timestamp
% find all GPS fixes which are older than 1 minute and discard them
%
if ~isfield(op,'gps_overrun')
  gps_dn = datenum(2000+data.m_gps_utc_year,data.m_gps_utc_month,data.m_gps_utc_day,...
    data.m_gps_utc_hour,data.m_gps_utc_minute,data.m_gps_utc_second);
  present_time = func_webbtime2mattime(data.m_present_time);
  ind = find(present_time-gps_dn>1/1440);
  if ~isempty(ind)
    disp(['Found ',int2str(length(ind)),' gps fixes older than one minute. Discarding all of them.'])
    pause(2)
    gps_lon(ind) = nan;
    gps_lat(ind) = nan;
    inv_lon(ind) = nan;
    inv_lat(ind) = nan;
    ign_lon(ind) = nan;
    ign_lat(ind) = nan;
    far_lon(ind) = nan;
    far_lat(ind) = nan;
    bad_list.old_fix = length(ind);
  else
    bad_list.old_fix = 0;
  end
else
  disp('op.gps_overrun was set. Not checking for correctness of GPS timestamps.')
end


%
% find bad values and replace them by NaN
%
bad = find(gps_lon<lon_range(1) | gps_lon>lon_range(2) | (gps_lat==0 & gps_lon==0) |...
    gps_lat>lat_range(2) | gps_lat<lat_range(1));
if ~isempty(bad)
  gps_lon(bad) = nan;
  gps_lat(bad) = nan;
  bad_list.gps_range = length(bad);
else
  bad_list.gps_range = 0;
end
bad = find(far_lon<lon_range(1) | far_lon>lon_range(2) | (far_lat==0 & far_lon==0) |...
    far_lat>lat_range(2) | far_lat<lon_range(1));
if ~isempty(bad)
  far_lon(bad) = nan;
  far_lat(bad) = nan;
  bad_list.gps_far_range = length(bad);
else
  bad_list.gps_far_range = 0;
end
bad = find(ign_lon<lon_range(1) | ign_lon>lon_range(2) | (ign_lat==0 & ign_lon==0) |...
    ign_lat>lat_range(2) | ign_lat<lat_range(1));
if ~isempty(bad)
  ign_lon(bad) = nan;
  ign_lat(bad) = nan;
  bad_list.gps_ign_range = length(bad);
else
  bad_list.gps_ign_range = 0;
end
bad = find(inv_lon<lon_range(1) | inv_lon>lon_range(2) | (inv_lat==0 & inv_lon==0) |...
    inv_lat>lat_range(2) | inv_lat<lat_range(1));
if ~isempty(bad)
  inv_lon(bad) = nan;
  inv_lat(bad) = nan;
  bad_list.gps_inv_range = length(bad);
else
  bad_list.gps_inv_range = 0;
end
bad = find(lon<lon_range(1) | lon>lon_range(2) | (lat==0 & lon==0) | lat>lat_range(2) | lat<lat_range(1));
if ~isempty(bad)
  lon(bad) = nan;
  lat(bad) = nan;
  bad_list.dr_range = length(bad);
else
  bad_list.dr_range = 0;
end


% 
% catch cases when there is an invalid GPS and a good GPS stored at the same time
% Then we throw both away.
%
ind = find( ~isnan(inv_lon) & ~isnan(gps_lon) );
if ~isempty(ind)
  disp(['Found ',int2str(length(ind)),' times with invalid and valid positions. Discarding all of them.'])
  gps_lon(ind) = nan;
  gps_lat(ind) = nan;
  inv_lon(ind) = nan;
  inv_lat(ind) = nan;
  ign_lon(ind) = nan;
  ign_lat(ind) = nan;
  far_lon(ind) = nan;
  far_lat(ind) = nan;
end
ind = find( ~isnan(ign_lon) & ~isnan(gps_lon) );
if ~isempty(ind)
  disp(['Found ',int2str(length(ind)),' times with ignored and valid positions. Discarding all of them.'])
  gps_lon(ind) = nan;
  gps_lat(ind) = nan;
  inv_lon(ind) = nan;
  inv_lat(ind) = nan;
  ign_lon(ind) = nan;
  ign_lat(ind) = nan;
  far_lon(ind) = nan;
  far_lat(ind) = nan;
end
ind = find( ~isnan(far_lon) & ~isnan(gps_lon) );
if ~isempty(ind)
  disp(['Found ',int2str(length(ind)),' times with toofar and valid positions. Discarding all of them.'])
  gps_lon(ind) = nan;
  gps_lat(ind) = nan;
  inv_lon(ind) = nan;
  inv_lat(ind) = nan;
  ign_lon(ind) = nan;
  ign_lat(ind) = nan;
  far_lon(ind) = nan;
  far_lat(ind) = nan;
end


%
% remove all position info from very short bd files as they
% often are bad and stem only from Iridium calls, when there is
% no regular GPS reception
%
[ids] = unique(data.dos_id);
for n=1:length(ids)
  ind = find(ids(n)==data.dos_id);
  if length(ind)<20
    disp(' ')
    disp('Found short BD file. Setting these positions to NaN.')
    disp(['Yo number : ',int2str(n),'   removed values : ',int2str(length(ind))])
    inv_lat(ind) = nan;
    inv_lon(ind) = nan;
    gps_lat(ind) = nan;
    gps_lon(ind) = nan;
    ign_lat(ind) = nan;
    ign_lon(ind) = nan;
    far_lat(ind) = nan;
    far_lon(ind) = nan;
    lat(ind) = nan;
    lon(ind) = nan;
  end
end


%
% Calculate the apparent speed of the glider between GPS fixes the glider trusts.
% Occasionally the glider will trust a wrong position. They usually stem from a 
% previous surfacing and somehow find their way into the trusted GPS values.
% Could be that they are connected to aborts.
% Here we find them by looking for apparent speeds larger than 5m/s.
% Usually the first position of any bad position pair (needed to get the velocity)
% is the bad one. Here we set these to NaN.
%
% If more than one value is given for a maximum allowed velocity, we repeat the
% checking process. This should take care of cases when up to length(max_vel)
% bad positions occur in sequence.
%
max_vel = [5,5];
allbad = [];
gps_lon_orig = gps_lon;
for m=1:length(max_vel)
  good = find(~isnan(gps_lon) & gps_lon~=-999);
  gps_u = cosd(gps_lat(good(1:end-1))).*diff(gps_lon(good))./diff(tim(good));
  gps_v = diff(gps_lat(good))/diff(tim(good));
  gps_vel = sqrt(gps_u.^2+gps_v.^2);
  gps_vel = gps_vel*60*1852;
  bad = find(gps_vel>max_vel(m));
  allbad = [allbad,bad];
  if ~isempty(bad)
    disp(' ')
    disp(['Found ',length(bad),' apparent velocities larger than ',int2str(max_vel(m)),' m/s.'])
    disp('Usually these are cases of ''leftover'' positions from the previous surfacing.')
    for n=1:length(bad)
      disp(['case #',int2str(n)])
      ind = bad(n)+[-50:50];
      ind = ind(find(ind>0 & ind<=length(good)));
      if op.no_plot~=1
        figure(2)
        clf
        subplot(3,1,1)
        title(['Case #',int2str(n)])
        plot(tim(good(ind)),gps_lat(good(ind)),'g.')
        xlim(tim(good(bad(n)))+[-60,60])
        xlabel('time')
        ylabel('latitude')
        hold on
        plot(tim(good(bad(n))),gps_lat(good(bad(n))),'rx')
        subplot(3,1,2)
        plot(tim(good(ind)),gps_lon(good(ind)),'g.')
        xlim(tim(good(bad(n)))+[-60,60])
        xlabel('time')
        ylabel('longitude')
        hold on
        plot(tim(good(bad(n))),gps_lon(good(bad(n))),'rx')
        subplot(3,1,3)
        plot((tim(good(ind(1:end-1)))+tim(good(ind(2:end))))/2,gps_vel(ind(1:end-1)),'x')
        xlim(tim(good(bad(n)))+[-60,60])
        disp('paused for 1 seconds')
        pause(1)
      end
    end
  end
  if op.no_plot~=1
    figure(1)
    clf
    subplot(2,1,1)
    title('Too large apparent speed check')
    plot(tim(good),gps_lat(good),'g.')
    xlabel('time')
    ylabel('latitude')
    hold on
    badbad = find(gps_lon==-999);
    plot(tim(badbad),gps_lat(badbad),'rx')
    subplot(2,1,2)
    plot(tim(good),gps_lon(good),'g.')
    xlabel('time')
    ylabel('longitude')
    hold on
    plot(tim(badbad),gps_lon_orig(badbad),'rx')
    pause(1)
    gps_lon(good(bad)) = -999;
  end
end    
bad = find(gps_lon==-999);
gps_lon(bad) = nan;
gps_lat(bad) = nan;

%
% Assume all glider-trusted GPS positions now to be ok and remove all positions in the other 
% position vectors (ignored, invalid, toofar, deadreckoning) that
% deviate more than 20 miles from all known trusted ones.
% This check will fail only if a dive was indeed longer than 20 miles.
%
% We are also applying this check to the 'toofar' positions. These positions should inherently 
% all fail this check. In rare circumstances the glider gets confused with his positions and
% might inadvertantly believe the wrong positions more than the good ones. There might thus
% be hidden good positions in the 'toofar' vector.
%
max_dist = 20;
good = find(~isnan(gps_lat));
if length(good)>1
  filled_gps_lat = interp1(good,gps_lat(good),[1:length(gps_lat)],'nearest','extrap');
  filled_gps_lon = interp1(good,gps_lon(good),[1:length(gps_lon)],'nearest','extrap');
elseif length(good)==1
  filled_gps_lat = gps_lat(good);
  filled_gps_lon = gps_lon(good);
else
  error('found no valid gps position')
end

% get the indices of the bad ones
bad_inv = find(sqrt( cosd(filled_gps_lat).^2.*(filled_gps_lat-inv_lat).^2+...
	(filled_gps_lon-inv_lon).^2)>max_dist/60);
bad_ign = find(sqrt( cosd(filled_gps_lat).^2.*(filled_gps_lat-ign_lat).^2+...
	(filled_gps_lon-ign_lon).^2)>max_dist/60);
bad_far = find(sqrt( cosd(filled_gps_lat).^2.*(filled_gps_lat-far_lat).^2+...
	(filled_gps_lon-far_lon).^2)>max_dist/60);
bad_dead = find(sqrt( cosd(filled_gps_lat).^2.*(filled_gps_lat-lat).^2+...
	(filled_gps_lon-lon).^2)>max_dist/60);

if op.no_plot~=1
  figure(2)
  clf
  subplot(2,1,1)
  plot(tim,lat,'g.')
  title('Positions too far away from trusted GPS positions')
  xlabel('time')
  ylabel('latitude')
  hold on
  plot(tim,inv_lat,'g.')
  plot(tim,ign_lat,'g.')
  plot(tim,far_lat,'g.')
  if ~isempty(bad_inv)
    disp(' ')
    disp(['Found ',int2str(length(bad_inv)),...
      '  m_gps_invalid_lat/lon  positions more than 20 miles away from trusted GPS.'])
    plot(tim(bad_inv),inv_lat(bad_inv),'rx')
  end
  if ~isempty(bad_ign)
    disp(' ')
    disp(['Found ',int2str(length(bad_ign)),...
      '  m_gps_ignored_lat/lon  positions more than 20 miles away from trusted GPS.'])
    plot(tim(bad_ign),ign_lat(bad_ign),'rx')
  end
  if ~isempty(bad_far)
    disp(' ')
    disp(['Found ',int2str(length(bad_far)),...
      '  m_gps_toofar_lat/lon  positions more than 20 miles away from trusted GPS.'])
    plot(tim(bad_far),far_lat(bad_far),'rx')
  end
  if ~isempty(bad_dead)
    disp(' ')
    disp(['Found ',int2str(length(bad_inv)),...
      '  m_lat/lon  positions more than 20 miles away from trusted GPS.'])
    plot(tim(bad_dead),lat(bad_dead),'rx')
  end
  subplot(2,1,2)
  plot(tim,lon,'g.')
  xlabel('time')
  ylabel('longitude')
  hold on
  plot(tim,ign_lon,'g.')
  plot(tim,inv_lon,'g.')
  plot(tim,far_lon,'g.')
  if ~isempty(bad_inv)
    plot(tim(bad_inv),inv_lon(bad_inv),'rx')
  end
  if ~isempty(bad_ign)
    plot(tim(bad_ign),ign_lon(bad_ign),'rx')
  end
  if ~isempty(bad_far)
    plot(tim(bad_far),far_lon(bad_far),'rx')
  end
  if ~isempty(bad_dead)
    plot(tim(bad_dead),lon(bad_dead),'rx')
  end
end

% finally set the bad ones to NaN
if ~isempty(bad_inv)
  inv_lat(bad_inv) = nan;
  inv_lon(bad_inv) = nan;
end
if ~isempty(bad_ign)
  ign_lat(bad_ign) = nan;
  ign_lon(bad_ign) = nan;
end
if ~isempty(bad_far)
  far_lat(bad_far) = nan;
  far_lon(bad_far) = nan;
end
if ~isempty(bad_dead)
  lat(bad_dead) = nan;
  lon(bad_dead) = nan;
end


%
% One sign of a really bad position is that it did not change at
% all from the previous position.
% Here we check for invalid positions that occur more than once
% and remove all deadreckoning positions that deviate by less than 20 cm.
% This >0 limit allows for rounding differences.
% For computation speed we do not correct for the latitudinal effect on
% the distance in x-direction. This should however not be a problem as
% it only makes the difference between positions smaller.
%
ii = inv_lat+sqrt(-1)*inv_lon;
good = find(~isnan(ii));
ig = lat+sqrt(-1)*lon;
[dummy,ind] = unique(ii(good));
for n=1:length(ind)
  ind2 = find( ii(good) == ii(good(ind(n))) );
  if length(ind2)>1
    ind3 = find( abs(ig-ii(good(ind(n))))<0.0001/60 );
    if ~isempty(ind3)
      disp(' ')
      disp(['Found ',int2str(length(ind3)),' positions when deadreckoned and invalid agree '])
      disp(['index ',int2str(ind3)])
      lat(ind3) = nan;
      lon(ind3) = nan;
    end
  end
end

figure(4)
clf
plot(tim,lat,'xb')
hold on
plot(tim,gps_lat,'r+')

%
% For unknown reasons deepy was sometimes unable to establish
% a proper GPS position. The result was a waiting time for a
% trusted fix with only invalid fixes that had the last trusted
% GPS position from the previous surfacing. After a timeout period
% it appears that the glider took one of the invalid positions
% as starting point by making it trusted with the new time stamp. 
% After that the glider
% continued to get GPS fixes and marked ALL of them toofar
% as these new ones had the correct position but only a very short time
% difference to the wrongly trusted position.
%
% Here we try to fix that by finding all cases where there is
% only a single m_gps_lat position within 1150 scans.
% After that single wrong valid fix all following good
% fixes are declared toofar by the glider.
% Here this one m_gps_lat fix is declared invalid and
% the toofar fixes are moved over to trusted m_gps_lat.
% The following deadrecoking fixes are corrected by the 
% same shift.
%
badvalid = find(~isnan(gps_lat));
for m=1:length(badvalid)
  surround = badvalid(m)+[-150:1000];
  surround = surround(find(surround>=1 & surround<=length(dep)));
  goodvalid = find(~isnan(gps_lat(surround)));
  if length(goodvalid)>1
    badvalid(m) = nan;
  end
  if ~isnan(badvalid(m))
    if badvalid(m)<length(ign_lat)
      ign = ign_lat(badvalid(m)+[-20:1]);
    else
      ign = ign_lat(badvalid(m)+[-20:0]);
    end
    if any( abs( gps_lat(badvalid(m)) - ign ) > 1e-6 )
      badvalid(m) = nan;
    end
  end
end
badvalid = badvalid(find(~isnan(badvalid)));
for m=1:length(badvalid)
  disp(' ')
  disp('Found bad valid position')
  if op.no_plot~=1
    figure(3)
    clf
    surround = badvalid(m)+[-150:10000];
    surround = surround(find(surround>=1 & surround<=length(dep)));
    subplot(2,1,1)
    plot(tim(surround),gps_lat(surround),'r.','markersize',20)
    hold on
    plot(tim(surround),lat(surround),'.','markersize',20)
    plot(tim(surround),inv_lat(surround),'xg')
    plot(tim(surround),far_lat(surround),'xk')
    plot(tim(badvalid(m)),gps_lat(badvalid(m)),'r.','markersize',20)
    xlabel('r. valid   b. deadreckoned   xg invalid   xk toofar')  
  end
  maxind = min([10000,length(lat)-badvalid(m)-1]);
  ind_dr = find(~isnan(lat(badvalid(m)+[1:maxind])));
  count = 1;
  while count>=1
    if count<length(ind_dr)
      check_dr = find(~isnan(lat(badvalid(m)+ind_dr(count)+[0:20])));
      num_dr = length(check_dr);
      if num_dr<15
%      disp(num_dr)
%      disp(lat(badvalid(m)+ind_dr(count)))
%disp(badvalid(m)+ind_dr(count))
        lat(badvalid(m)+ind_dr(count)+[0:20]) = nan; 
        lon(badvalid(m)+ind_dr(count)+[0:20]) = nan; 
        count = count+1;
      else
        count = 0;
      end
    else
      count = 0;
    end
  end
  next_dr = badvalid(m)+min(find(~isnan(lat(badvalid(m)+1:end))));
  if ~isempty(next_dr)
    corr_lat = nmean(far_lat(next_dr+[-2:2]))-lat(next_dr);
    corr_lon = nmean(far_lon(next_dr+[-2:2]))-lon(next_dr);
    ind = next_dr+[0:10000];
    ind = ind(find(ind<=length(lat)));
    dummy = find(~isnan(lat(ind)));
    cont_dr = dummy(find(diff(dummy)>10));
    if ~isempty(cont_dr)
      cont_dr = next_dr+[0:cont_dr(1)+1];
      if op.no_plot~=1
        plot(tim(cont_dr),lat(cont_dr),'+c')
        plot(tim(next_dr),lat(next_dr),'xr')
      end
    end
  else
    corr_lat = 0;
    corr_lon = 0;
  end

  if op.no_plot~=1
    subplot(2,1,2)
    plot(tim(badvalid(m)),gps_lat(badvalid(m)),'+r','markersize',10)
    hold on
  end

  if ~isempty(next_dr)
    lat(cont_dr) = lat(cont_dr)+corr_lat;
    lon(cont_dr) = lon(cont_dr)+corr_lon;
  end
  maxind = min([1000,length(far_lat)-badvalid(m)-1]);
  toofar = find(~isnan(far_lat(badvalid(m)+[1:maxind])));
  gps_lat(badvalid(m)+toofar) = far_lat(badvalid(m)+toofar);
  gps_lon(badvalid(m)+toofar) = far_lon(badvalid(m)+toofar);
  gps_lat(badvalid(m)) = nan;
  gps_lon(badvalid(m)) = nan;
  lat(badvalid(m)) = nan;
  lon(badvalid(m)) = nan;

  if op.no_plot~=1
    plot(tim(surround),gps_lat(surround),'r.','markersize',20)
    if ~isempty(next_dr)
      plot(tim(cont_dr),lat(cont_dr),'+c')
      plot(tim(next_dr),lat(next_dr),'xr')
      plot(tim(cont_dr),lat(cont_dr),'+m')
      disp([gps_lat(badvalid(m)),lat(next_dr),corr_lat,corr_lon])
    end
    plot(tim(surround),lat(surround),'.','markersize',20)
    plot(tim(surround),inv_lat(surround),'xg')
    plot(tim(surround),far_lat(surround),'xk')
    plot(tim(surround),gps_lat(surround),'r.','markersize',20)
  end
end

%figure(5)
%clf
%plot(tim,lat,'xb')
%hold on
%plot(tim,gps_lat,'r+')
%pause


%
% prepare result vectors
%
nlat = gps_lat;
nlon = gps_lon;


%
% Look for all gaps in the GPS positions that are longer
% than 10 minutes, assume that these are dives. Here we correct the respective
% dead-reckoned positions linearly from 0 at the beginning to max at the
% end so that at the end the dea-reckoned positions agree with the new
% GPS fixes.
%
disp(' ')
disp('Modifying dead-reckoned positions.')
good = find(~isnan(gps_lon));
dgood = diff(good);
for n=1:length(dgood)

  %
  % if there is a large gap (>10 mins) in GPS data we assume a dive
  %
  if dgood(n)>op.min_dive_length		% a value of 250 is about 10 minutes

    % get the indices for the dive (last and first time there is a GPS
    % position
    st = good(n);
    en = good(n+1);
    ind = [st:en];

 
    % get the indices of all the GPS fixes within 2min after the 
    % first GPS fix after the dive
    extra_ind = en+[0:30];
    extra_ind = extra_ind( find(extra_ind<=length(gps_lon)) );
    extra_ind = extra_ind( find( tim(extra_ind)-tim(extra_ind(1)) < 120 ) );
    goodgps = en + find(~isnan(gps_lon(extra_ind))) - 1;


    % get the indices when the glider did dead-reckoning
    % during the dive
    gooddr = find(~isnan(lat(ind)));


    % this is the last time during the dive when the glider did
    % dead-reckoning. This is NOT the surface, but about 7m depth.
    % At this moment the battery moves and the movement of the 
    % glider becomes unclear. Surfacing occurs 1-2minutes later
    if ~isempty(gooddr)
      lastdr = ind(gooddr(end));
      if lastdr==en & length(gooddr)>1
        lastdr = ind(gooddr(end-1));
      end

      % thus the best guess for the time of the surfacing is
      indsurf = lastdr + 20;

      % if there are multiple GPS fixes we calculate the surface drift
      % and extrapolate back to the time of the surfacing
      % if there are too few fixes, we just take these fixes
      % if there are no fixes, we do not correct
      if op.use_median==1
        extra_ind = en+[0:500];
        extra_ind = extra_ind( find(extra_ind<=length(gps_lon)) );
        goodgps = en + find(~isnan(gps_lon(extra_ind))) - 1;
        fit_gps_lat = nmedian(gps_lat(goodgps));
        fit_gps_lon = nmedian(gps_lon(goodgps));
      elseif length(goodgps)>5
        warning off  % with time values so high and close to each other Matlab thinks they are all the same...
        uu = (gps_lon(goodgps(end))-gps_lon(goodgps(1)))/(tim(goodgps(end))-tim(goodgps(1)))*1852*60;
        vv = (gps_lat(goodgps(end))-gps_lat(goodgps(1)))/(tim(goodgps(end))-tim(goodgps(1)))*1852*60;
        if abs(uu)>5 | abs(vv)>5
          disp('Found velocity > 5m/s')
          keyboard
        end
        plon = polyfit(tim(goodgps),gps_lon(goodgps),1);
        fit_gps_lon = polyval(plon,tim([indsurf,extra_ind]));
        plat = polyfit(tim(goodgps),gps_lat(goodgps),1);
        fit_gps_lat = polyval(plat,tim([indsurf,extra_ind]));

        warning on
      elseif isempty(goodgps)
        fit_gps_lat = nan;
        fit_gps_lon = nan;
      else
        fit_gps_lat = gps_lat(goodgps(1));
        fit_gps_lon = gps_lon(goodgps(1));
      end
    
      % now correct the dead-reckoned positions in a way that
      % the start point agrees with GPS (that's always the case,
      % since the last GPS fix is used as start of the dead-reckoning)
      % and the end point of the linearly corrected dead-reckoing
      % agrees with the estimated surfacing position
      if (length(gooddr)>30 & ~isnan(fit_gps_lon)) | all(dep(gooddr)<3)
        lastp = lastdr;
        firstp = ind( gooddr(1) );
        indp = [firstp:lastp];
        da = lat(lastp) - lat(firstp);
        do = lon(lastp) - lon(firstp);
        dt = tim(lastp) - tim(firstp);
        dga = fit_gps_lat(1) - gps_lat(st);
        dgo = fit_gps_lon(1) - gps_lon(st);
        dgt = tim(lastp) - tim(st);
        dat = da/dt;
        dot = do/dt;
        dgat = dga/dgt;
        dgot = dgo/dgt;
        t = tim(indp)-tim(indp(1));
        corr_lat = t*(dgat-dat);
        corr_lon = t*(dgot-dot);
        bad = find( isnan( gps_lat(indp) ) );
        nlat(indp(bad)) = lat(indp(bad)) + corr_lat(bad);
        nlon(indp(bad)) = lon(indp(bad)) + corr_lon(bad);
      else
        disp('Only short or shallow dead-reckoning interval. Likely at the surface. Not correcting.')
        if all(dep(ind)<3)
          disp('all depths shallower than 3m')
        end
        if length(gooddr)
          disp('found fewer than 30 dead reckoned positions')
        end
        disp('dos_ids')
        disp([data.dos_id(ind(1)),data.dos_id(ind(end))])
        figure(1)
        clf
        plot(tim(ind),lon(ind),'b+')
        hold on
        plot(tim(ind),gps_lon(ind),'r+','markersize',20)
        plot(tim(extra_ind),gps_lon(extra_ind),'rx')
        plot(tim(st-[1:30]),gps_lon(st-[1:30]),'rx')
        plot(tim([indsurf,extra_ind]),fit_gps_lon,'k-')
        xlabel('glider time')
        ylabel('longitude [degrees]')
        title('b+: dead reckoned  r+: used gps  rx: other gps  k-:fit')
        ax = axis;
        figure(3)
        clf
        plot(tim(ind),dep(ind),'b+')
        xlabel('glider time')
        ylabel('depth [m]')
        ax2 = axis;
        axis([ax(1:2),ax2(3:4)])
        pause(1)
        nlat(gooddr) = lat(gooddr);
        nlon(gooddr) = lon(gooddr);
      end 
    else
      disp('Problem with nav. No deadreckoning position.') 
    end

    % debugging plot
    if 0
      figure(1)
      clf
      plot(tim(ind),lon(ind),'b+')
      hold on
      plot(tim(ind),gps_lon(ind),'r+')
      plot(tim(extra_ind),gps_lon(extra_ind),'rx')
      plot(tim(st-[1:30]),gps_lon(st-[1:30]),'rx')
      plot(tim([indsurf,extra_ind]),fit_gps_lon,'k-')
      plot(tim(lastp),lon(lastp),'gx')
      plot(tim(indp),nlon(indp),'go')
      figure(2)
      clf
      plot(tim(ind),lat(ind),'b+')
      hold on
      plot(tim(ind),gps_lat(ind),'r+')
      plot(tim(extra_ind),gps_lat(extra_ind),'rx')
      plot(tim(st-[1:30]),gps_lat(st-[1:30]),'rx')
      plot(tim([indsurf,extra_ind]),fit_gps_lat,'k-')
      plot(tim(lastp),lat(lastp),'gx')
      plot(tim(indp),nlat(indp),'go')
      ax = axis;
      figure(3)
      clf
      plot(tim(ind),dep(ind),'b+')
      ax2 = axis;
      axis([ax(1:2),ax2(3:4)])
      pause
    end


  end
end  

data.n_lat = nlat;
data.n_lon = nlon;


%
% clean up data
%
%clear lon dos_id dep gps_lon inv_lon far_lon ign_lon gps_u gps_v gps_vel
%clear dgood ii ig allgpslat allgpslon


if op.no_plot~=1
  figure(1)
  clf
  plot(tim,orig_lat,'m+')
  hold on
  plot(tim,lat,'.','markersize',1)
  plot(tim,gps_lat,'r+')
  plot(tim,inv_lat,'bx')
  plot(tim,far_lat,'c+')
  plot(tim,ign_lat,'kx')
  plot(tim,nlat,'g.')
  xlabel('g.: new   m+: old   r+: measured  bx: invalid  c+: toofar  kx: ignored')
  ylabel('latitude [degrees]')
  print('-djpeg',['plots',filesep,'step_05_best_position_data_lat.jpg'])
  figure(2)
  clf
  plot(tim,orig_lon,'m+')
  hold on
  plot(tim,lon,'.','markersize',1)
  plot(tim,gps_lon,'r+')
  plot(tim,gps_lon,'rx')
  plot(tim,inv_lon,'bx')
  plot(tim,far_lon,'c+')
  plot(tim,ign_lon,'kx')
  plot(tim,nlon,'g.')
  xlabel('g.: new   m+: old   r+: measured  bx: invalid  c+: toofar  kx: ignored')
  ylabel('longitude [degrees]')
  print('-djpeg',['plots',filesep,'step_05_best_position_data_lon.jpg'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%]
% input variables are coming in a structure. The fields of the structure
% can be cell structures. Here some processing is applied to each element
% of the cell structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%]

function [latout,lonout] = parse_structures(latin,lonin)

if ~iscell(latin)
  if issparse(latin)
    latin = func_slocum_sparse(latin);
    lonin = func_slocum_sparse(lonin);
  end
  latout = func_slocum_degmin(latin);
  lonout = func_slocum_degmin(lonin);
else
  for n=1:length(latin)
    [latout{n},lonout{n}] = parse_structures(latin{n},lonin{n})
  end
end
