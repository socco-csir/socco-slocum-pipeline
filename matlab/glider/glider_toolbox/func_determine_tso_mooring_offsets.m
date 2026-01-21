function [] = func_determine_mooring_tso_offsets(reference_file)
% function [] = func_determine_mooring_tso_offsets()
% 
% GEOMAR SVN $Id: func_determine_tso_mooring_offsets.m 303 2017-03-06 08:52:37Z gkrahmann@geomar.de $
%
% compare the glider T,S,O data against mooring data 
%
% version 2  last change 24.06.2015

% G.Krahmann, GEOMAR,  Jun 2014

% run once with applied offsets and once without  GK, 24.06.2015  1-->2


%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load the mooring reference data
%
load(reference_file);
alat = nlat;
alon = nlon;
at = nt;
as = ns;
ao = no;
ap = np;
atim = ntim;


%
% load glider data and rearrange
%
data = load([op.deplname,'_gridded'],'latitude','longitude','main_datenum','ctd_temperature',...
    'aanderaa_oxygen_calculated_undelayed_filtered',...
    'geomar_oxygen_calculated_undelayed_filtered','good_salinity','pressure');
main_datenum = data.main_datenum;
latitude = data.latitude;
longitude = data.longitude;
temperature = data.ctd_temperature;
temperature_offset = data.ctd_temperature - op.t_offset;
oxygen_a = data.aanderaa_oxygen_calculated_undelayed_filtered;
oxygen_a_offset = data.aanderaa_oxygen_calculated_undelayed_filtered - op.oa_offset;
if ~isfield(data,'geomar_oxygen_calculated_undelayed_filtered')
  oxygen_g = nan*oxygen_a;
  oxygen_g_offset = nan*oxygen_a;
else
  oxygen_g = data.geomar_oxygen_calculated_undelayed_filtered;
  oxygen_g_offset = data.geomar_oxygen_calculated_undelayed_filtered - op.og_offset;
end
pressure = data.pressure;
salinity = data.good_salinity;
salinity_offset = data.good_salinity - op.s_offset;
lat1 = nmean(latitude);
lon1 = nmean(longitude);
tim1 = nmean(main_datenum);


%
% loop over different moorings
%
for k=1:length(alat)

  %
  % extract single mooring
  %
  ntim = atim{k};
  nlat = alat{k};
  nlon=  alon{k};
  np = ap{k};
  nt = at{k};
  ns = as{k};
  no = ao{k};

  %
  % reduce mooring data to glider deployment period
  %
  ind = find( ntim>=nmin(main_datenum(:)) & ntim<=nmax(main_datenum(:)));
  ntim = ntim(ind);
  nt = nt(:,ind);
  ns = ns(:,ind);
  no = no(:,ind);
  np = np(:,ind);

  %
  % determine distance between glider and mooring
  % and loop over several maximal distances
  %
  dis_gl_mo = sqrt( (lat1-nlat(1)).^2 + cosd(nlat(1))^2*(lon1-nlon(1)).^2 )*60;
  distances = [1,2,5];
  for ndis = [1:length(distances)];

    ind = find(dis_gl_mo<distances(ndis));

    for m=1:length(ind)
    
      [dtim,indtim] = nmin( abs( tim1(ind(m))-ntim ) );

      if dtim<=1/24

        for k=1:size(np,1)
      
          moor_press = round(np(k,indtim)); 
          indpress = find(pressure==moor_press);
          if length(indpress)==1
            res_p(ndis,k,m) = np(k,indtim);
            res_t(ndis,k,m) = temperature(indpress,ind(m)) - nt(k,indtim);
            res_s(ndis,k,m) = salinity(indpress,ind(m)) - ns(k,indtim);
            res_oa(ndis,k,m) = oxygen_a(indpress,ind(m)) - no(k,indtim);
            res_og(ndis,k,m) = oxygen_g(indpress,ind(m)) - no(k,indtim);
            res_t_offset(ndis,k,m) = temperature_offset(indpress,ind(m)) - nt(k,indtim);
            res_s_offset(ndis,k,m) = salinity_offset(indpress,ind(m)) - ns(k,indtim);
            res_oa_offset(ndis,k,m) = oxygen_a_offset(indpress,ind(m)) - no(k,indtim);
            res_og_offset(ndis,k,m) = oxygen_g_offset(indpress,ind(m)) - no(k,indtim);
          else
            res_p(ndis,k,m) = nan;
            res_t(ndis,k,m) = nan;
            res_s(ndis,k,m) = nan;
            res_oa(ndis,k,m) = nan;
            res_og(ndis,k,m) = nan;
            res_t_offset(ndis,k,m) = nan;
            res_s_offset(ndis,k,m) = nan;
            res_oa_offset(ndis,k,m) = nan;
            res_og_offset(ndis,k,m) = nan;
            res_n(ndis,k,m) = nan;
          end

        end
      end
    end
  end

  if ~exist('res_t')
    disp('glider was not close to mooring')
    disp(' ')
  else

    res_p = nans(res_p,nan,0,'==');
    res_t = nans(res_t,nan,0,'==');
    res_s = nans(res_s,nan,0,'==');
    res_oa = nans(res_oa,nan,0,'==');
    res_og = nans(res_og,nan,0,'==');

    %
    % display results
    %
    disp('------------------------------------------------------------------------ ')
    disp('------------------------------------------------------------------------ ')
    disp(['Offsets Glider minus Mooring set #',int2str(k)])
    disp(' ')
    disp('Temperature')
    disp('depth   <= 1nm               <=2nm                <=5nm                    ')
    indgood = find(~isnan(nmean(nmean(res_t,3),1)));
    for n=indgood
      fprintf(1,'%s \n',[int2str0(nmean(res_p(1,n,:)),4,' '),'m   ',...
                         num2str(nmean(res_t(1,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_t(1,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_t(2,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_t(2,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_t(3,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_t(3,n,:)),'%5.4f'),...
      ])
      mres1(n) = nmean(res_t(1,n,:));
      mres2(n) = nmean(res_t(2,n,:));
      mres3(n) = nmean(res_t(3,n,:));
      sres1(n) = nstd(res_t(1,n,:));
      sres2(n) = nstd(res_t(2,n,:));
      sres3(n) = nstd(res_t(3,n,:));
    end
    fprintf(1,'%s \n',['        ',int2str0(sum(~isnan(res_s(1,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(2,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(3,indgood(1),:))),3,' values')]);
    mres1 = nans(mres1,nan,0,'==');
    mres2 = nans(mres2,nan,0,'==');
    mres3 = nans(mres3,nan,0,'==');
    warning off
    w1 = meanweight(mres1,1./sres1,0);
    w2 = meanweight(mres2,1./sres2,0);
    w3 = meanweight(mres3,1./sres3,0);
    warning on
    disp('weighted average')
    fprintf(1,'%s \n',[num2str(w1,'%+6.4f'),...
      '              ',             num2str(w2,'%+6.4f'),...
      '              ',             num2str(w3,'%+6.4f'),...
      ])
  
    disp(' ')
    disp('------------------------------------------------------------------------ ')
    clear mres1 mres2 mres3 sres1 sres2 sres3
    disp('Salinity')
    disp('depth   <= 1nm               <=2nm                <=5nm                    ')
    indgood = find(~isnan(nmean(nmean(res_s,3),1)));
    for n=indgood
      fprintf(1,'%s \n',[int2str0(nmean(res_p(1,n,:)),4,' '),'m   ',...
                         num2str(nmean(res_s(1,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_s(1,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_s(2,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_s(2,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_s(3,n,:)),'%+6.4f'),' +/- ',num2str(nstd(res_s(3,n,:)),'%5.4f'),...
      ])
      mres1(n) = nmean(res_s(1,n,:));
      mres2(n) = nmean(res_s(2,n,:));
      mres3(n) = nmean(res_s(3,n,:));
      sres1(n) = nstd(res_s(1,n,:));
      sres2(n) = nstd(res_s(2,n,:));
      sres3(n) = nstd(res_s(3,n,:));
    end
    fprintf(1,'%s \n',['        ',int2str0(sum(~isnan(res_s(1,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(2,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(3,indgood(1),:))),3,' values')]);
    mres1 = nans(mres1,nan,0,'==');
    mres2 = nans(mres2,nan,0,'==');
    mres3 = nans(mres3,nan,0,'==');
    warning off
    w1 = meanweight(mres1,1./sres1,0);
    w2 = meanweight(mres2,1./sres2,0);
    w3 = meanweight(mres3,1./sres3,0);
    warning on
    disp('weighted average')
    fprintf(1,'%s \n',['        ',num2str(w1,'%+6.4f'),...
      '              ',             num2str(w2,'%+6.4f'),...
      '              ',             num2str(w3,'%+6.4f'),...
      ])
  
    disp(' ')
    clear mres1 mres2 mres3 sres1 sres2 sres3
    disp('------------------------------------------------------------------------ ')
    disp('Oxygen Aanderaa')
    disp('depth   <= 1nm               <=2nm                <=5nm                    ')
    indgood = find(~isnan(nmean(nmean(res_oa,3),1)));
    for n=indgood
      fprintf(1,'%s \n',[int2str0(nmean(res_p(1,n,:)),4,' '),'m   ',...
                         num2str(nmean(res_oa(1,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_oa(1,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_oa(2,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_oa(2,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_oa(3,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_oa(3,n,:)),'%5.4f'),...
      ])
      mres1(n) = nmean(res_oa(1,n,:));
      mres2(n) = nmean(res_oa(2,n,:));
      mres3(n) = nmean(res_oa(3,n,:));
      sres1(n) = nstd(res_oa(1,n,:));
      sres2(n) = nstd(res_oa(2,n,:));
      sres3(n) = nstd(res_oa(3,n,:));
    end
    fprintf(1,'%s \n',['        ',int2str0(sum(~isnan(res_s(1,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(2,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(3,indgood(1),:))),3,' values')]);
    disp(' ')
    mres1 = nans(mres1,nan,0,'==');
    mres2 = nans(mres2,nan,0,'==');
    mres3 = nans(mres3,nan,0,'==');
    warning off
    w1 = meanweight(mres1,1./sres1,0);
    w2 = meanweight(mres2,1./sres2,0);
    w3 = meanweight(mres3,1./sres3,0);
    warning on
    disp('weighted average')
    fprintf(1,'%s \n',[num2str(w1,'%+6.4f'),...
      '              ',             num2str(w2,'%+6.4f'),...
      '              ',             num2str(w3,'%+6.4f'),...
      ])
    clear mres1 mres2 mres3 sres1 sres2 sres3
    disp('------------------------------------------------------------------------ ')
    disp('Oxygen GEOMAR')
    disp('depth   <= 1nm               <=2nm                <=5nm                    ')
    indgood = find(~isnan(nmean(nmean(res_og,3),1)));
    for n=indgood
      fprintf(1,'%s \n',[int2str0(nmean(res_p(1,n,:)),4,' '),'m   ',...
                         num2str(nmean(res_og(1,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_og(1,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_og(2,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_og(2,n,:)),'%5.4f'),...
      '   ',             num2str(nmean(res_og(3,n,:)),'%+7.3f'),' +/- ',num2str(nstd(res_og(3,n,:)),'%5.4f'),...
      ])
      mres1(n) = nmean(res_og(1,n,:));
      mres2(n) = nmean(res_og(2,n,:));
      mres3(n) = nmean(res_og(3,n,:));
      sres1(n) = nstd(res_og(1,n,:));
      sres2(n) = nstd(res_og(2,n,:));
      sres3(n) = nstd(res_og(3,n,:));
    end
    fprintf(1,'%s \n',['        ',int2str0(sum(~isnan(res_s(1,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(2,indgood(1),:))),3,' '),'                  ',...
                       int2str0(sum(~isnan(res_s(3,indgood(1),:))),3,' values')]);
    disp(' ')
    mres1 = nans(mres1,nan,0,'==');
    mres2 = nans(mres2,nan,0,'==');
    mres3 = nans(mres3,nan,0,'==');
    warning off
    w1 = meanweight(mres1,1./sres1,0);
    w2 = meanweight(mres2,1./sres2,0);
    w3 = meanweight(mres3,1./sres3,0);
    warning on
    disp('weighted average')
    fprintf(1,'%s \n',[num2str(w1,'%+6.4f'),...
      '              ',             num2str(w2,'%+6.4f'),...
      '              ',             num2str(w3,'%+6.4f'),...
      ])
  end 
end
