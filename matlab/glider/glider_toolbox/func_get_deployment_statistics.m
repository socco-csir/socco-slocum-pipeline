function [] = func_get_deployment_statistics()
% function [] = func_get_deployment_statistics()
% 
% GEOMAR SVN $Id: func_get_deployment_statistics.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% read _yos.mat file and extract some statistics
%
% version 1  last change 20.02.2013

% G.Krahmann, GEOMAR Feb 2013

%
% look for and load processing parameters
%
if exist('./processing_parameters.m') | exist('.\processing_parameters.m')
  processing_parameters
else
  error(['could not find file ''processing_parameters.m'' in current folder ',pwd])
end


%
% load data
%
load([op.deplname,'_yos']);
data = load([op.deplname,'_1vector']);


disp(' ')
disp(['Deployment : ',op.deplname])
disp(['Yos        : ',int2str(length(m_pressure))])

coulomb1 = nmin(data.m_coulomb_amphr);
coulomb2 = nmax(data.m_coulomb_amphr);
coul = coulomb2-coulomb1;
if isnan(coul)
  coul = '-';
else
  coul = int2str(round(coul));
end

tim1 = webbtime2mattime(nmin(m_present_time{1}));
tim2 = webbtime2mattime(nmax(m_present_time{end}));
dv1 = datestr(datevec(tim1),'DD.mm.YYYY');
dv2 = datestr(datevec(tim2),'DD.mm.YYYY');
 
dx = 0;
ca2 = 0;
ca1 = 1000000;
di2 = 0;
di1 = 1000000;
for n=1:length(m_pressure)
  la = n_lat{n};
  lo = n_lon{n};
  good = find(~isnan(la));
  if ~isempty(good)
    la1 = la(good(1));
    lo1 = lo(good(1));
    la2 = la(good(end));
    lo2 = lo(good(end));
    dx = dx + sqrt((la2-la1)^2 + cosd(la1)^2*(lo2-lo1)^2);
  end
end
dx = dx*1.852*60;
dx2 = 0;
for n=3:2:length(m_pressure)
  la = n_lat{n-2};
  lo = n_lon{n-2};
  good = find(~isnan(la));
  if ~isempty(good)
    la1 = la(good(1));
    lo1 = lo(good(1));
  else
    la1 = nan;
    lo1 = nan;
  end
  la = n_lat{n};
  lo = n_lon{n};
  good = find(~isnan(la));
  if ~isempty(good)
    la2 = la(good(1));
    lo2 = lo(good(1));
  else
    la2 = nan;
    lo2 = nan;
  end
  dx2 = dx2 + sqrt((la2-la1)^2 + cosd(la1)^2*(lo2-lo1)^2);
end
dx2 = dx2*1.852*60;
if dx2>dx
  dx = dx2;
end
dx = int2str(dx);

if ~exist('m_iridium_call_num')
  ca = 'x';
  di = 'x';
else
keyboard
  ca2 = nmax(m_iridium_call_num);
  ca1 = nmin(m_iridium_call_num);
  di2 = nmax(m_iridium_dialed_num);
  di1 = nmin(m_iridium_dialed_num);
  ca = int2str(ca2-ca1);
  di = int2str(di2-di1);
end

disp(['Days       : ',int2str(round(tim2-tim1))])
disp(['Energy     : ',coul,' Ahr'])
disp(['Distance   : ',dx,' km'])
disp(['Dialed     : ',di])
disp(['Call       : ',ca])

disp(['| ',depl(end+[-1,0]),' || ',dv1,' - ',dv2,' || ',int2str(round(tim2-tim1)),...
	' || ',int2str(length(data)),' || ',dx,' km || ',coul,' Ah || ',ca,' / ',di,' || ?'])
%keyboard
