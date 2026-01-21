function [] = func_analyze_tso_data(glider,ref,comp_type,op,reference_file)
% analyze the already set up T,S,O differences
% 
% GEOMAR SVN $Id: func_analyze_tso_data.m 890 2021-12-15 15:38:48Z gkrahmann@geomar.de $
%
% function [] = func_analyze_tso_data(glider,ref,comp_type,op,reference_file)
%
% version 4.1.0  last change 17.11.2021

% G.Krahmann, GEOMAR,  Feb 2014

% rewrite for standard data sets                GK, 16.10.2014  1-->2
% handle PC/UNIX and folder plots               GK, 22.03.2017  2-->3
% allow for different maximum ctd-glider distance. Previously this was fixed to
% 0.2 degrees.                                  GK, 16.01.2018  3-->4
% changed header                                GK, 17.11.2021  4-->4.1.0

%
% check whether we have two different oxygen data sets (that is reference data from another glider)
% or only a single one (that is reference data from a mooring or CTDs). In case it is only
% a single data set, we simply double it.
%
% if isfield(ref,'no')
%   ref.noa = ref.no;
%   ref.nog = ref.no;
%   ref.noad = ref.nod;
%   ref.nogd = ref.nod;
% end


%
% handle on press or dens
%
if ~isempty(findstr(comp_type,'dens'))
  checks = [[1022:0.05:1025],[1025.01:0.01:1033]];
else
  checks = [10:10:1000];
end


%
% reduce reference data to overlapping data
%
gl_tmin = nmin(glider.ntim);
gl_tmax = nmax(glider.ntim);
gl_latmin = nmin(glider.nlat);
gl_latmax = nmax(glider.nlat);
gl_lonmin = nmin(glider.nlon);
gl_lonmax = nmax(glider.nlon);
lat_factor = cosd(nmean(glider.nlat));
if ~isfield(op,'max_ctd_comp_dist')
  max_dist = 0.2;
else
  max_dist = op.max_ctd_comp_dist;
end
ind = find( ref.ntim+2>=gl_tmin & ref.ntim-2<=gl_tmax &...
  ref.nlat+max_dist/lat_factor>gl_latmin & ref.nlat-max_dist/lat_factor<gl_latmax &...
  ref.nlon+max_dist/lat_factor>gl_lonmin & ref.nlon-max_dist/lat_factor<gl_lonmax);
if isempty(ind)
  disp('found no suitable reference data')
end
ref.ntim = ref.ntim(ind);
ref.nlat = ref.nlat(ind);
ref.nlon = ref.nlon(ind);
ref.np = ref.np(:,ind);
ref.nt = ref.nt(:,ind);
ref.ns = ref.ns(:,ind);
% ref.noa = ref.noa(:,ind);
% ref.nog = ref.nog(:,ind);
% ref.ntd = ref.ntd(:,ind);
% ref.nsd = ref.nsd(:,ind);
% ref.noad = ref.noad(:,ind);
% ref.nogd = ref.nogd(:,ind);



%
% create time and distance matrix
%
ds = repmat(nan,length(ref.ntim),length(glider.ntim));
dt = repmat(nan,length(ref.ntim),length(glider.ntim));
[i_glider,i_ref] = meshgrid([1:length(glider.ntim)],[1:length(ref.ntim)]);
for n=1:length(ref.ntim)
  ds(n,:) = sqrt( lat_factor^2*(ref.nlon(n)-glider.nlon).^2 + (ref.nlat(n)-glider.nlat).^2 );
  dt(n,:) = abs( ref.ntim(n) - glider.ntim );
end


%
% loop over max time deviations
%
disp('------------------------------------------------------------------------------------ ')
count = 1;
for check_time=[2,1,0.5,0.2]

  good = find( dt < check_time & ds < max_dist );

  % limit this to something manageable
  if length(good)>1e6
    good = good(1:10:end);
  elseif length(good)>1e5
    good = good(1:3:end);
  end

  if ~isempty(good)
    comp_p = repmat(nan,length(checks),length(good));
    comp_t = repmat(nan,length(checks),length(good));
    comp_td = repmat(nan,length(checks),length(good));
    comp_s = repmat(nan,length(checks),length(good));
    comp_sd = repmat(nan,length(checks),length(good));
%     comp_og = repmat(nan,length(checks),length(good));
%     comp_ogd = repmat(nan,length(checks),length(good));
%     comp_oa = repmat(nan,length(checks),length(good));
%     comp_oad = repmat(nan,length(checks),length(good));
    comp_ds = nan*ones(1,length(good));
    comp_dt = nan*ones(1,length(good));
  else
    comp_p = repmat(nan,length(checks),1);
    comp_t = repmat(nan,length(checks),1);
    comp_td = repmat(nan,length(checks),1);
    comp_s = repmat(nan,length(checks),1);
    comp_sd = repmat(nan,length(checks),1);
%     comp_og = repmat(nan,length(checks),1);
%     comp_ogd = repmat(nan,length(checks),1);
%     comp_oa = repmat(nan,length(checks),1);
%     comp_oad = repmat(nan,length(checks),1);
    comp_ds = nan;
    comp_dt = nan;
  end
% v1 = comp_og;
% v2 = comp_og;

  for n=1:length(good)

    comp_p(:,n) = glider.np(:,i_glider(good(n)))-ref.np(:,i_ref(good(n)));
    comp_t(:,n) = glider.nt(:,i_glider(good(n)))-ref.nt(:,i_ref(good(n)));
%     comp_td(:,n) = glider.ntd(:,i_glider(good(n)))-ref.ntd(:,i_ref(good(n)));
    comp_s(:,n) = glider.ns(:,i_glider(good(n)))-ref.ns(:,i_ref(good(n)));
%     comp_sd(:,n) = glider.nsd(:,i_glider(good(n)))-ref.nsd(:,i_ref(good(n)));
%     comp_og(:,n) = glider.nog(:,i_glider(good(n)))-ref.nog(:,i_ref(good(n)));
%     comp_ogd(:,n) = glider.nogd(:,i_glider(good(n)))-ref.nogd(:,i_ref(good(n)));
%     comp_oa(:,n) = glider.noa(:,i_glider(good(n)))-ref.noa(:,i_ref(good(n)));
%     comp_oad(:,n) = glider.noad(:,i_glider(good(n)))-ref.noad(:,i_ref(good(n)));
    comp_ds(n) = ds(good(n));
    comp_dt(n) = dt(good(n));

%     v1(:,n) = glider.nogd(:,i_glider(good(n)));
%     v2(:,n) = ref.nogd(:,i_ref(good(n)));

  end

  disp(['calculating for time difference shorter than : ',num2str(check_time),' days'])
  dummy = ones(length(checks),1)*comp_ds;

  %
  % calculate mean TS data by averaging all nearby profiles
  %
  mean_dis = nmean(comp_ds);
  ind_gli = unique(i_glider(good));
  glider_nearby_mean_t = nmean(glider.nt(:,ind_gli)')';
  glider_nearby_mean_s = nmean(glider.ns(:,ind_gli)')';
  ind_ref = unique(i_ref(good));
  ref_nearby_mean_t = nmean(ref.nt(:,ind_ref)')';
  ref_nearby_mean_s = nmean(ref.ns(:,ind_ref)')';
  if strcmp('ctd_on_dens',comp_type)
    figure(9)
  else
    figure(10)
  end
  if count==1
    clf
  end
  subplot(2,2,count)
  plot(glider_nearby_mean_s,glider_nearby_mean_t,'k')
  hold on
  plot(ref_nearby_mean_s,ref_nearby_mean_t,'r')
  grid on
  orient landscape
  ax = axis;
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.9,['# gli : ',int2str(length(ind_gli))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.8,['# ref : ',int2str(length(ind_ref))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.7,['mean dis : ',int2str(mean_dis*60),'nm'])



  dtmax = 0.5;
  good = find(abs(comp_t)<dtmax);
  bad = find(abs(comp_t)>dtmax);
  good_deep = find(abs(comp_td)<dtmax);
  bad_deep = find(abs(comp_td)>dtmax);
  if ~isempty(good)
    res(count,1) = nmean(comp_t(good));
    res(count,2) = nmedian(comp_t(good));
    res(count,3) = meanweight(comp_t(good),1./(dummy(good)+0.05),1);
    dis(count,1) = nmean(dummy(good));
    dis(count,2) = nmedian(dummy(good));
    dis(count,3) = meanweight(dummy(good),1./(dummy(good)+0.05),1);
  else
    res(count,1:3) = [nan,nan,nan];
    dis(count,1:3) = [nan,nan,nan];
  end
  if ~isempty(good_deep)
    res_deep(count,1) = nmean(comp_td(good_deep));
    res_deep(count,2) = nmedian(comp_td(good_deep));
    res_deep(count,3) = meanweight(comp_td(good_deep),1./(dummy(good_deep)+0.05),1);
    dis_deep(count,1) = nmean(dummy(good_deep));
    dis_deep(count,2) = nmedian(dummy(good_deep));
    dis_deep(count,3) = meanweight(dummy(good_deep),1./(dummy(good_deep)+0.05),1);
  else
    res_deep(count,1:3) = [nan,nan,nan];
    dis_deep(count,1:3) = [nan,nan,nan];
  end

if 0
  figure(1)
  clf
  plot(dummy,comp_t,'x')
end

  figure(2)
  if count==1
    clf
  end
  subplot(2,2,count)
  [aa,bb] = hist(comp_t(good),[-dtmax:0.01:dtmax]);
  [dummy2,indmax] = nmax(aa);
  res2(count,1) = bb(indmax);
  hist(comp_t(good),[-dtmax:0.01:dtmax])
  h = findobj(gca,'Type','patch');
  hold on
  hist(comp_td(good_deep),[-dtmax:0.01:dtmax])
  set(h,'FaceColor','r','EdgeColor','w')
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none')
  ax = axis;
  axis([dtmax*[-1.2,1.2],ax(3:4)])
  plot(res(count,3)*[1,1],ax(3:4),'g','linewidth',1)
  plot(res_deep(count,3)*[1,1],ax(3:4),'k','linewidth',1)
  ax = axis;
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.9,['av-all : ',num2str(res(count,3))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.8,['av-deep: ',num2str(res_deep(count,3))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.7,['di-all : ',int2str(dis(count,3)*60),'nm'])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.6,['di-deep: ',int2str(dis_deep(count,3)*60),'nm'])
  ax = axis;
  title(['\Delta time < ',num2str(check_time),' days    ',...
    int2str(length(bad)/(length(good)+length(bad))*100),'% of values out of range'])
  grid on


  dsmax = 0.2;
  good = find(abs(comp_s)<dsmax);
  bad = find(abs(comp_s)>dsmax);
  good_deep = find(abs(comp_sd)<dsmax);
  bad_deep = find(abs(comp_sd)>dsmax);
  if ~isempty(good)
    res(count,4) = nmean(comp_s(good));
    res(count,5) = nmedian(comp_s(good));
    res(count,6) = meanweight(comp_s(good),1./(dummy(good)+0.05),1);
    dis(count,4) = nmean(dummy(good));
    dis(count,5) = nmedian(dummy(good));
    dis(count,6) = meanweight(dummy(good),1./(dummy(good)+0.05),1);
  else
    res(count,4:6) = [nan,nan,nan];
    dis(count,4:6) = [nan,nan,nan];
  end
  if ~isempty(good_deep)
    res_deep(count,4) = nmean(comp_sd(good_deep));
    res_deep(count,5) = nmedian(comp_sd(good_deep));
    res_deep(count,6) = meanweight(comp_sd(good_deep),1./(dummy(good_deep)+0.05),1);
    dis_deep(count,4) = nmean(dummy(good_deep));
    dis_deep(count,5) = nmedian(dummy(good_deep));
    dis_deep(count,6) = meanweight(dummy(good_deep),1./(dummy(good_deep)+0.05),1);
  else
    res_deep(count,4:6) = [nan,nan,nan];
    dis_deep(count,4:6) = [nan,nan,nan];
  end

if 0 
  figure(3)
  clf
  plot(dummy,comp_s,'x')
end

  figure(4)
  if count==1
    clf
  end
  subplot(2,2,count)
  [aa,bb] = hist(comp_s(good),[-dsmax:0.005:dsmax]);
  [dummy2,indmax] = nmax(aa);
  res2(count,2) = bb(indmax);
  hist(comp_s(good),[-dsmax:0.005:dsmax])
  h = findobj(gca,'Type','patch');
  hold on
  hist(comp_sd(good_deep),[-dsmax:0.005:dsmax])
  set(h,'FaceColor','r','EdgeColor','w')
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none')
  ax = axis;
  axis([dsmax*[-1.2,1.2],ax(3:4)])
  plot(res(count,6)*[1,1],ax(3:4),'g','linewidth',1)
  plot(res_deep(count,6)*[1,1],ax(3:4),'k','linewidth',1)
  ax = axis;
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.9,['av-all : ',num2str(res(count,6))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.8,['av-deep: ',num2str(res_deep(count,6))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.7,['di-all : ',int2str(dis(count,6)*60),'nm'])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.6,['di-deep: ',int2str(dis_deep(count,6)*60),'nm'])
  title(['\Delta time < ',num2str(check_time),' days    ',...
    int2str(length(bad)/(length(good)+length(bad))*100),'% of values out of range'])
  grid on

  domax = 50;
  good = find(abs(comp_og)<domax);
  bad = find(abs(comp_og)>domax);
  good_deep = find(abs(comp_ogd)<domax);
  bad_deep = find(abs(comp_ogd)>domax);
  if ~isempty(good)
    res(count,7) = nmean(comp_og(good));
    res(count,8) = nmedian(comp_og(good));
    res(count,9) = meanweight(comp_og(good),1./(dummy(good)+0.05),1);
    dis(count,7) = nmean(dummy(good));
    dis(count,8) = nmedian(dummy(good));
    dis(count,9) = meanweight(dummy(good),1./(dummy(good)+0.05),1);
  else
    res(count,7:9) = [nan,nan,nan];
    dis(count,7:9) = [nan,nan,nan];
  end
  if ~isempty(good_deep)
    res_deep(count,7) = nmean(comp_ogd(good_deep));
    res_deep(count,8) = nmedian(comp_ogd(good_deep));
    res_deep(count,9) = meanweight(comp_ogd(good_deep),1./(dummy(good_deep)+0.05),1);
    dis_deep(count,7) = nmean(dummy(good_deep));
    dis_deep(count,8) = nmedian(dummy(good_deep));
    dis_deep(count,9) = meanweight(dummy(good_deep),1./(dummy(good_deep)+0.05),1);
  else
    res_deep(count,7:9) = [nan,nan,nan];
    dis_deep(count,7:9) = [nan,nan,nan];
  end

if 0
  figure(5)
  clf
  plot(dummy,comp_og,'x')
end

  figure(6)
  if count==1
    clf
  end
  subplot(2,2,count)
  [aa,bb] = hist(comp_og(good),[-domax:0.2:domax]);
  [dummy2,indmax] = nmax(aa);
  res2(count,3) = bb(indmax);
  hist(comp_og(good),[-domax:0.2:domax])
  h = findobj(gca,'Type','patch');
  hold on
  hist(comp_ogd(good_deep),[-domax:0.2:domax])
  set(h,'FaceColor','r','EdgeColor','w')
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none')
  ax = axis;
  axis([domax*[-1.2,1.2],ax(3:4)])
  plot(res(count,9)*[1,1],ax(3:4),'g','linewidth',1)
  plot(res_deep(count,9)*[1,1],ax(3:4),'k','linewidth',1)
  ax = axis;
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.9,['av-all : ',num2str(res(count,9))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.8,['av-deep: ',num2str(res_deep(count,9))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.7,['di-all : ',int2str(dis(count,9)*60),'nm'])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.6,['di-deep: ',int2str(dis_deep(count,9)*60),'nm'])
  title(['\Delta time < ',num2str(check_time),' days    ',...
    int2str(length(bad)/(length(good)+length(bad))*100),'% of values out of range'])
  grid on



  good = find(abs(comp_oa)<domax);
  bad = find(abs(comp_oa)>domax);
  good_deep = find(abs(comp_oad)<domax);
  bad_deep = find(abs(comp_oad)>domax);
  if ~isempty(good)
    res(count,10) = nmean(comp_oa(good));
    res(count,11) = nmedian(comp_oa(good));
    res(count,12) = meanweight(comp_oa(good),1./(dummy(good)+0.05),1);
    dis(count,10) = nmean(dummy(good));
    dis(count,11) = nmedian(dummy(good));
    dis(count,12) = meanweight(dummy(good),1./(dummy(good)+0.05),1);
  else
    res(count,10:12) = [nan,nan,nan];
    dis(count,10:12) = [nan,nan,nan];
  end
  if ~isempty(good_deep)
    res_deep(count,10) = nmean(comp_oad(good_deep));
    res_deep(count,11) = nmedian(comp_oad(good_deep));
    res_deep(count,12) = meanweight(comp_oad(good_deep),1./(dummy(good_deep)+0.05),1);
    dis_deep(count,10) = nmean(dummy(good_deep));
    dis_deep(count,11) = nmedian(dummy(good_deep));
    dis_deep(count,12) = meanweight(dummy(good_deep),1./(dummy(good_deep)+0.05),1);
  else
    res_deep(count,10:12) = [nan,nan,nan];
    dis_deep(count,10:12) = [nan,nan,nan];
  end

if 0
  figure(7)
  clf
  plot(dummy,comp_oa,'x')
end

  figure(8)
  if count==1
    clf
  end
  subplot(2,2,count)
  [aa,bb] = hist(comp_oa(good),[-domax:0.2:domax]);
  [dummy2,indmax] = nmax(aa);
  res2(count,4) = bb(indmax);
  hist(comp_oa(good),[-domax:0.2:domax])
  h = findobj(gca,'Type','patch');
  hold on
  hist(comp_oad(good_deep),[-domax:0.2:domax])
  set(h,'FaceColor','r','EdgeColor','w')
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none')
  ax = axis;
  axis([domax*[-1.2,1.2],ax(3:4)])
  plot(res(count,12)*[1,1],ax(3:4),'g','linewidth',1)
  plot(res_deep(count,12)*[1,1],ax(3:4),'k','linewidth',1)
  ax = axis;
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.9,['av-all : ',num2str(res(count,12))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.8,['av-deep: ',num2str(res_deep(count,12))])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.7,['di-all : ',int2str(dis(count,12)*60),'nm'])
  text(ax(1)+(ax(2)-ax(1))*0.05,ax(3)+(ax(4)-ax(3))*0.6,['di-deep: ',int2str(dis_deep(count,12)*60),'nm'])
  title(['\Delta time < ',num2str(check_time),' days    ',...
    int2str(length(bad)/(length(good)+length(bad))*100),'% of values out of range'])
  grid on

%keyboard

  count = count+1;
end

disp(' ')
disp('Differences glider data minus reference data')
disp([op.deplname,'  minus  ',reference_file])
disp(comp_type)
disp(' ')
ind = find(~isnan(res(:,3)));
fprintf(1,['Temperature   :  <2.0d     <1.0d     <0.5d     <0.2d\n']);
if ~isempty(ind)
  fprintf(1,['Weighted      : %7.4f   %7.4f   %7.4f   %7.4f\n'],...
    res(:,3))
  fprintf(1,['Weighted deep : %7.4f   %7.4f   %7.4f   %7.4f\n'],...
    res_deep(:,3))
else
  disp(['Weighted      : NaN'])
  disp(['Weighted deep : NaN'])
end

disp(' ')
ind = find(~isnan(res(:,6)));
fprintf(1,['Salinity      :  <2.0d     <1.0d     <0.5d     <0.2d\n']);
if ~isempty(ind)
  fprintf(1,['Weighted      : %7.4f   %7.4f   %7.4f   %7.4f\n'],...
    res(:,6))
  fprintf(1,['Weighted deep : %7.4f   %7.4f   %7.4f   %7.4f\n'],...
    res_deep(:,6))
else
  disp(['Weighted      : NaN'])
  disp(['Weighted deep : NaN'])
end

disp(' ')
ind = find(~isnan(res(:,12)));
fprintf(1,['Ox Aanderaa   :  <2.0d     <1.0d     <0.5d     <0.2d\n']);
if ~isempty(ind)
  fprintf(1,['Weighted      : %7.2f   %7.2f   %7.2f   %7.2f\n'],...
    res(:,12))
  fprintf(1,['Weighted deep : %7.2f   %7.2f   %7.2f   %7.2f\n'],...
    res_deep(:,12))
else
  disp(['Weighted      : NaN'])
  disp(['Weighted deep : NaN'])
end

disp(' ')
ind = find(~isnan(res(:,9)));
fprintf(1,['Ox GEOMAR     :  <2.0d     <1.0d     <0.5d     <0.2d\n']);
if ~isempty(ind)
  fprintf(1,['Weighted      : %7.2f   %7.2f   %7.2f   %7.2f\n'],...
    res(:,9))
  fprintf(1,['Weighted deep : %7.2f   %7.2f   %7.2f   %7.2f\n'],...
    res_deep(:,9))
else
  disp(['Weighted      : NaN'])
  disp(['Weighted deep : NaN'])
end

if ~exist('plots','dir')
  mkdir plots
end
figure(2)
suplabel([op.deplname,' Temperature differences -glider + ',comp_type],'t');
orient landscape
%print('-djpeg',['plots/t_offsets_',comp_type,'.jpg'])

figure(4)
suplabel([op.deplname,' Salinity differences -glider + ',comp_type],'t');
orient landscape
%print('-djpeg',['plots/s_offsets_',comp_type,'.jpg'])

figure(6)
suplabel([op.deplname,' Oxygen (GEOMAR) differences -glider + ',comp_type],'t');
orient landscape
%print('-djpeg',['plots/og_offsets_',comp_type,'.jpg'])

figure(8)
suplabel([op.deplname,' Oxygen (Aanderaa) differences -glider + ',comp_type],'t');
orient landscape
%print('-djpeg',['plots/oa_offsets_',comp_type,'.jpg'])

figure(9)
%suplabel([op.deplname,' TS diagr glider(k)  ',comp_type,'(r)'],'t');
orient landscape
if ~exist('plots','dir')
  mkdir plots
end
%print('-djpeg',['plots',filesep,'ts_',comp_type,'.jpg'])

figure(10)
%suplabel([op.deplname,' TS diagr glider(k)  ',comp_type,'(r)'],'t');
orient landscape
if ~exist('plots','dir')
  mkdir plots
end
%print('-djpeg',['plots',filesep,'ts_',comp_type,'.jpg'])





function [ax,h]=suplabel(text,whichLabel,supAxes)
% PLaces text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.
%  [ax,h]=suplabel(text,whichLabel,supAxes)
% returns handles to both the axis and the label.
%  ax=suplabel(text,whichLabel,supAxes)
% returns a handle to the axis only.
%  suplabel(text) with one input argument assumes whichLabel='x'
%
% whichLabel is any of 'x', 'y', or 't', specifying whether the 
% text is to be the xlable, ylabel, or title respectively.
%
% supAxes is an optional argument specifying the Position of the 
%  "super" axes surrounding the subplots. 
%  supAxes defaults to [.075 .075 .85 .85]
%  specify supAxes if labels get chopped or overlay subplots
%
% EXAMPLE:
%  subplot(2,2,1);ylabel('ylabel1');title('title1')
%  subplot(2,2,2);ylabel('ylabel2');title('title2')
%  subplot(2,2,3);ylabel('ylabel3');xlabel('xlabel3')
%  subplot(2,2,4);ylabel('ylabel4');xlabel('xlabel4')
%  [ax,h1]=suplabel('super X label');
%  [ax,h2]=suplabel('super Y label','y');
%  [ax,h3]=suplabel('super Title'  ,'t');
%  set(h3,'FontSize',30)
%
% SEE ALSO: text, title, xlabel, ylabel, zlabel, subplot,
%           suptitle (Matlab Central)

% Author: Ben Barrowes <barrowes@alum.mit.edu>
% modified to move axes to the back as described in Mathworks fileexchange
% turn of latex interpreter

if nargin < 3
 supAxes=[.08 .08 .84 .84];
 ah=findall(gcf,'type','axes');
 if ~isempty(ah)
  supAxes=[inf,inf,0,0];
  leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
  axBuf=.04;
  set(ah,'units','normalized')
  ah=findall(gcf,'type','axes');
  for ii=1:length(ah)
   if strcmp(get(ah(ii),'Visible'),'on')
    thisPos=get(ah(ii),'Position');
    leftMin=min(leftMin,thisPos(1));
    bottomMin=min(bottomMin,thisPos(2));
    leftMax=max(leftMax,thisPos(1)+thisPos(3));
    bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
   end
  end
  supAxes=[leftMin-axBuf,bottomMin-axBuf,leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
 end
end
if nargin < 2, whichLabel = 'x';  end
if nargin < 1, help(mfilename); return; end

if ~isstr(text) | ~isstr(whichLabel)
  error('text and whichLabel must be strings')
end
whichLabel=lower(whichLabel);

ax=axes('Units','Normal','Position',supAxes,'Visible','off');
ch = get(get(ax, 'parent'), 'children');
set(get(ax, 'parent'), 'children', [ch(2:end); ch(1)]);
if strcmp('t',whichLabel)
  h=get(ax,'Title');
  set(h,'Visible','on','interpreter','none')
  title(text);
elseif strcmp('x',whichLabel)
  h=get(ax,'XLabel');
  set(h,'Visible','on')
  xlabel(text);
elseif strcmp('y',whichLabel)
  h=get(ax,'YLabel');
  set(h,'Visible','on')
  ylabel(text);
end
