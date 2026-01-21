function out = opt_evalcal(temp,phase_diff,cal_method,foilcoef,p,sal,phasefunstr,beta,pcfactor,O2_ref,ind)
%%
% function out = opt_evalcal(temp,phase_diff,cal_method,foilcoef,p,sal,phasefunstr,beta,pcfactor,O2_ref,ind)
% 
% GEOMAR SVN $Id: opt_evalcal.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% =========================================================================
%
% function to evaluate optode oxygen calibration using fitting method from
% Uchida and additional fit on 'BPhase'-'RPhase'('phase_diff') of optode. 'phase_diff' calibration
% coefficients ('beta') and calibration function ('phasefunstr') are
% calculated with the function 'optdoinsitufit.m'
%
% INPUT:
% ------
%   temp        ... temperature (degreeC)
%   phase_diff  ... optode phase difference (BPhase-RPhase) (degree)
%   cal_method  ... calibration method ( = 'aanderaa' or 'uchida' etc.)
%   foilcoef    ... Uchida refit or Aanderaa fit of (Aanderaa) multipoint calibration
%                   (surface)
%   p           ... pressure (dbar)
%   sal         ... salinity (PSU)
%   phasefunstr ... string of phase calibration function
%                   (calibration of BPhase)
%                       X(:,1) - temperature
%                       X(:,2) - phase_diff
%                       X(:,3) - Pressure
%   beta        ... calibration coefficients for phasefunstr
%                       beta(:,1) - mean values
%                       beta(:,2) - confidence interval (95%)
%   pcfactor    ... pressure correction factor
%   O2_ref      ... (optional) reference oxygen; if given, then compared to
%                   calibrated optode oxygen
%   ind         ... (optional) indices, on which calibrated oxygen shall be
%                   compared with 'O2_ref'
%                   
%       ----> temp, phase_diff, p, sal, (O2_ref) all must have same size
%
% OUTPUT:
% -------
%   out.O2fit           ... calibrated oxygen [mumol/kg]
%   out.rms_kgpc        ... (optional - only if 'O2_ref' is given)
%                           error of 'O2fit' to 'O2_ref' (mumol/kg) for all
%                           values
%   out.rms_kgpc_ind    ... (optional - only if 'O2_ref' and 'ind' is given)
%                           error of 'O2fit' to 'O2_ref' (mumol/kg) for
%                           given indices
%
% version 1.0       last change 31.01.2012
%
% Author:
%   J.Hahn, IfM GEOMAR, Kiel, Jan. 2012
%
%   initial coding                  J.Hahn      Jan. 2012   v1.0
%
% =========================================================================

%% calibrate 'phase_diff'

% % prepare foil coefficients

foilstr=num2str(foilcoef(:,1),'%.8g');
foilstr_row = ['[' foilstr(1,:)];
for i=2:size(foilstr,1)
    foilstr_row = [foilstr_row ';' foilstr(i,:)];
end
foilstr_row = [foilstr_row ']'];

% % prepare 'phase_diff' calibration function

eval(['phasefun=@(beta,X)' phasefunstr ';'])

% % determine calibrated oxygen

% eval(['O2fit_l=optcalcO2(temp,phasefun(beta(:,1),[temp phase_diff p]),' foilstr_row ',''uchida'',sal,1013.25);']) % surface saltwater fit
eval(['O2fit_l=optcalcO2(temp,phasefun(beta(:,1),[temp phase_diff p]),' foilstr_row ',''',cal_method,''',sal,1013.25);']) % surface saltwater fit
O2fit_kg = molar2molal(O2fit_l,temp,sal,p);
O2fit_kg_pc = optprescorr(O2fit_kg,p,pcfactor);

%% calculate error

if exist('O2_ref','var')
    
%     O2_ref_l = optreverseprescorr(molal2molar(O2_ref,temp,sal,p),p,pcfactor);
%     out.rms_l = rms(O2_ref_l-O2fit_l);      % uniform rms (in l,surf)
    out.rms_kgpc = rms(O2_ref-O2fit_kg_pc); % uniform rms (in kg,pc)
    
    if exist('ind','var')
    %     out.rms_l_ind = rms(O2_ref_l(ind)-O2fit_l(ind));      % uniform rms for indeces (in l,surf)
        out.rms_kgpc_ind = rms(O2_ref(ind)-O2fit_kg_pc(ind)); % uniform rms for indeces (in kg,pc)
    end
end

out.O2fit = O2fit_kg_pc;
