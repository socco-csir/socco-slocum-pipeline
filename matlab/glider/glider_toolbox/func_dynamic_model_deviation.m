function [res]=dyn_optim(par)
% function [res]=dyn_optim(par);
% 
% GEOMAR SVN $Id: func_dynamic_model_deviation.m 703 2020-07-28 14:50:02Z gkrahmann@geomar.de $
%
% glider flight model cost function to be optimized 
%
% input  : par                     - parameter vector to be optimized
%
% output : res                     - cost value
%
% globals: val                     - structure with all relevant non-changing data
%          plotnew                 - refresh the plot (this is not done each optimization step)
%          handle_model_graph      - handle to the line plot of the flight model
%          w_glider_model          - resulting modeled vertical speed
%          nnn                     - dive number (required when drag is time dependent)
%          facs                    - scaling factors for the parameter vector to be roughly 1
%          fixed_params            - structure of values for parameters that are not to be optimized
%          param_list              - structure-fieldnames of parameters to be optimized
%          attack                  - resulting angle of attack
%          pitch                   -
%
% version 2.1.0,  last change 31.01.2023

% G.Krahmann, GEOMAR, April 2013

% handle backward movement                               GK, 31.01.2017  1-->2
% change from CSIRO seawater to TEOS library             GK, 31.01.2023  2-->2.1.0

global val plotnew handle_model_graph w_glider_model nnn facs fixed_params param_list attack pitch no_plot


% scale the parameters that are optimized to their original size
if ~isempty(par)
  par = par./facs;
end

% loop over the parameter list and figure out which ones do NOT need to be optimized
% because a value is already given in the structure  fixed_params
count = 1;
for n=1:length(param_list)
  if isempty(getfield(fixed_params,param_list{n}))
    if strcmp('C_D_0',param_list{n})
      C_D_0 = par(count);
    elseif strcmp('epsilon',param_list{n})
      epsilon = par(count);
    elseif strcmp('m_g',param_list{n})
      m_g = par(count);
    elseif strcmp('alpha_T',param_list{n})
      alpha_T = par(count);
    elseif strcmp('temperature_filter',param_list{n})
      temperature_filter = par(count);
    end
    count = count+1;
  else
    if strcmp('C_D_0',param_list{n})
      if length(fixed_params.C_D_0)==1
        C_D_0 = fixed_params.C_D_0;
      else
        C_D_0 = fixed_params.C_D_0(nnn);
      end
    elseif strcmp('epsilon',param_list{n})
      epsilon = fixed_params.epsilon;
    elseif strcmp('m_g',param_list{n})
      if length(fixed_params.m_g)==1
        m_g = fixed_params.m_g;
      else
        m_g = fixed_params.m_g(nnn);
      end
    elseif strcmp('alpha_T',param_list{n})
      alpha_T = fixed_params.alpha_T;
    elseif strcmp('temperature_filter',param_list{n})
      temperature_filter = fixed_params.temperature_filter;
    end
  end
end

S = 0.10;
a_w = 3.7;
a_h = 2.4;
a_h = a_w * (val.A_h/S) + 1.2;
C_D_1w = 0.78;
C_D_1h = 2.1;
g = gsw_grav(val.latitude);
T_0 = 16;

V_g = val.V_g;
theta = val.pitch;
rho = val.rho;
shear_lift = val.shear_lift;
P = val.P;
T = val.T;
w = val.w;
delta_V_bp = val.delta_V_bp;
%fin = val.fin;
good = val.good;

T = val.T;

alpha = 0;
for m=1:10
  alpha = (C_D_0 + (C_D_1w+C_D_1h)*alpha.^2)./((a_w+a_h)*tan(theta+alpha));
  bad = find(abs(alpha+theta)>pi);
  % catch runaway cases because of too large C_D_0 in the optimization
  if ~isempty(bad)
    alpha(bad) = sign(theta(bad))*pi-theta(bad);
  end
end
pitch = theta;
attack = alpha;
gamma = alpha+theta;

% add this only to the real drag, not to the AoA calculation
%C_D_0 = C_D_0 + abs(fin/200);  % 200 looks ok
%C_D_0 = C_D_0 + fin.^2/8000;

F_g = m_g * g;
F_b = g * rho .* ( V_g * (1 - epsilon*P +  alpha_T*( T - T_0 ) ) + delta_V_bp );
t1 = F_b - F_g + shear_lift;
t2 = 1/2 * rho * S .* ( C_D_0 + ( C_D_1w + C_D_1h )*alpha.^2 ) .*...
     (sin(gamma).^2+cos(gamma).^2)./sin(gamma);
sig = sign( t1 ./ t2 );
U = sig .* sqrt( abs( t1 ./ t2 ) );

%u_g = U .* cos( gamma );
w_glider_model = U .* sin( gamma );

% calculate cost function
%res = rms(w-w_glider_model);
good = find(imag(U)==0);
res = sum((w(good)-w_glider_model(good)).^2);
if isnan(res)
w(good)
res
end

% add some penalties
res = res + length(find(imag(U)~=0));
if epsilon<0
  res = res * (1-epsilon*1e8);
end
if C_D_0<0
  res = res * (1-C_D_0);
end
if alpha_T<0
  res = res * (1-alpha_T*1e4);
end
if temperature_filter<0
  res = res * (1-temperature_filter*1e3);
end
if alpha_T>1e-3
  res = res * (1+alpha_T*1e2);
end
if m_g>72
  res = res * (1+m_g-72);
end
if m_g<52
  res = res * (1-m_g+52);
end

if isnan(res) & imag(res)==0
  disp('res is NaN')
  keyboard
end



if plotnew==floor(plotnew/20)*20
  set(handle_model_graph,'ydata',[w_glider_model]);
  drawnow
if 0
figure(2)
clf
plot(g*rho.*V_g.*epsilon.*P)
hold on
plot(g*rho.*V_g.*alpha_T.*(T-T_0),'r')
plot(g*rho.*delta_V_bp,'g')
drawnow
end
end
plotnew = plotnew + 1;
