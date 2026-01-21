function [oxygen]=optcalcO2(temp,phase,foilcoef,modeltype,sal,P_atm,P_dbar,pcfactor)
% function [oxygen]=optcalcO2(temp,phase,foilcoef,modeltype,sal,P_atm,P_dbar,pcfactor)
% 
% GEOMAR SVN $Id: optcalcO2.m 303 2017-03-06 08:52:37Z gkrahmann@geomar.de $
%
% apply polynomial fit/foil coefficients to optode phase and temperature 
% measurements according to specified model and calculate O2 in umol/l (or
% as desired)
% model: 
% 'aanderaa'   - 20 coefficient matrix 3rd order in t and 4th order
%                (via O2conc calculation only) (phase = DPhase/CalPhase!)
%                No pressure required
% '3x4b'       - 14 coefficient matrix by Craig Neill (via pO2 calculation)
%                (phase = BPhase/TCPhase)
% '5x5b'       - 21 coefficient matrix by Craig Neill (via pO2 calculation)
%                (phase = BPhase/TCPhase)
% 'uchida'     - 7 coefficients following Stern-Volmer equation derived
%                approach after Uchida et al. 2008
% 'uchidasimple'- 6 coefficients following Stern-Volmer equation derived
%                approach after Uchida et al. 2008, simplified after GO
%                SHIP Hydrography Manual (Version 1, 2010)
% 'any number' - upper polynomial matrix in bphase & t up to degree n
%                (1x1b, 2x2b, 3x3b, 4x4b) (via pO2 calculation)
% (if no modeltype is given size/shape of foilcoef is used as indicator)
%
% salinity can be of same length as temp/phase or singular value (default:
% sal = 0)
% atmospheric/gas pressure can be of same length as temp/phase or singular 
% value (default: P_atm = 1013.25 mbar) (=atmospheric pressure or 
% pressure from GTD if available)
% hydrostatic pressure can be of same length as temp/phase or singular
% value (default: P_dbar = 0 dbar)
% pcfactor (optional) gives value of pressure correction factor for optode
% foil (default: pcfactor = 3.2) (3.2% per 1000 dbar)
%
% relevant subfunctions at end of m-file
%
% part of optcalc-toolbox
% Henry Bittig, IFM-GEOMAR
% 01.11.2010
% revised 31.03.2011

% check and clean up input variables
if nargin<3
    disp('too few input variables')
    return
else
    if size(temp)~=size(phase)
        disp('variable sizes don''t match')
        return
    end
end

% remember shape of input
shape=size(temp);

if nargin<4
    %modeltype='aanderaa';
    [modeltype,foilcoef]=getmodeltype(foilcoef);
else
    modeltype=lower(modeltype);
    [n,m]=size(foilcoef);
    if m>n
        disp('foilcoef flipped')
        foilcoef=foilcoef';
    end
end

if nargin<5
    sal=zeros(shape);
else
    if find(~(size(sal)==shape)) %size(sal)~=shape
        if length(sal)==1
            sal=ones(shape).*sal;
            disp('expand salinity to match variable size')
        else
            disp('salinity size doesn''t match')
            return
        end
    end
end

if nargin<6
    P_atm=ones(shape).*1013.25; % mbar
else
    if size(P_atm)~=shape
        if length(P_atm)==1
            P_atm=ones(shape).*P_atm;
            %disp('expand atmospheric pressure to match variable size')
        else
            disp('atmospheric pressure size doesn''t match')
            return
        end
    end
end

if nargin<7
    P_dbar=zeros(shape); % dbar
else
    if size(P_dbar)~=shape
        if length(P_dbar)==1
            P_dbar=ones(shape).*P_dbar;
            %disp('expand hydrostatic pressure to match variable size')
        else
            disp('hydrostatic pressure size doesn''t match')
            return
        end
    end
end

if nargin<8
    pcfactor=3.2; % 3.2% per 1000 dbar as default for 3830 / 4330 standard foil
else
    if length(pcfactor)~=1
        disp('singular value for pcfactor expected!')
        return
    end
end

%make column vectors
temp=temp(:);
phase=phase(:);
%foilcoef=foilcoef(:);
sal=sal(:);
P_atm=P_atm(:);
P_dbar=P_dbar(:);

%do actual calculation according to specified modeltype
if strcmp(modeltype,'unknown')
    disp('unknown model: foilcoef doesn''t fit any expected size!')
    return
elseif (strcmp(modeltype,'aanderaa')|strcmp(modeltype,'3830'))
    %check number of foil coefficients
    if size(foilcoef,1)~=20
        if find(size(foilcoef)==[5 4])
            foilcoef=reshape(foilcoef',20,1);
        elseif find(size(foilcoef)==[4 5])
            foilcoef=reshape(foilcoef,20,1);
        else
            disp('number of foil coefficient doesn''t match model type!')
            return
        end
    end
    % oxygen concentration in umol/l from Aanderaa coefficients
    % freshwater (phase should be dphase!)
    oxy_conc0_fresh=optodefun3830(foilcoef,[temp,phase]);
    % salinity correction
    oxy_conc0=O2freshtosal(oxy_conc0_fresh,temp,sal);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
% elseif      strcmp(modeltype,'3830')
%    %check number of foil coefficients
%    if size(foilcoef,1)~=20
%        disp('number of foil coefficient doesn''t match model type!')
%        return
%    end
%    % pO2 in mbar from individual coefficients (3830 model polynom)
%    % (phase should be bphase)
%    pO2 = optodefun3830(foilcoef,[temp,phase]);
%    % oxygen concentration (salinity corrected) in umol/l
%    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
elseif  strcmp(modeltype,'3x4b')
    if size(foilcoef,1)~=14
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (3x4b model polynom)
    % (phase should be bphase)
    pO2 = optodefun3x4b(foilcoef,[temp,phase]);
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
elseif  strcmp(modeltype,'5x5b')
    if size(foilcoef,1)~=21
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (5x5b model polynom)
    % (phase should be bphase)
    pO2 = optodefun5x5b(foilcoef,[temp,phase]);
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
elseif  strcmp(modeltype,'uchida')
    if size(foilcoef,1)~=7
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (5x5b model polynom)
    % (phase should be bphase)
    pO2 = optodefunUchida(foilcoef,[temp,phase]);
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
elseif  strcmp(modeltype,'uchidasimple')
    if size(foilcoef,1)~=6
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (5x5b model polynom)
    % (phase should be bphase)
    pO2 = optodefunUchidasimple(foilcoef,[temp,phase]);
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
elseif  strcmp(modeltype,'3n')
    if size(foilcoef,1)~=10
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (5x5b model polynom)
    % (phase should be bphase)
    pO2 = optodefun3x3b(foilcoef,[temp,phase]);
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
elseif ~isnan(str2double(modeltype))
    n=str2double(modeltype);
    if size(foilcoef,1)~=(n+1)*(n+2)/2
        disp('number of foil coefficient doesn''t match model type!')
        return
    end
    % pO2 in mbar from individual coefficients (mixed polynom)
    % (phase should be bphase)
    pO2 = diag(polyval2(foilcoef',temp,phase));
    % oxygen concentration (salinity corrected) in umol/l
    oxy_conc0=pO2toO2conc(pO2,temp,sal,P_atm);
    % pressure correction to saltwater concentration
    oxy_conc0=optprescorr(oxy_conc0,P_dbar,pcfactor);
else
    disp('unknown model! Check input and size of foil coefficients.')
    return
end

% oxygen saturation (salinity corrected and pressure corrected (input)) 
% in % (use intake values if applicable)
%oxy_sat=O2conctoO2sat(oxy_conc0,temp,sal,P_atm);
% oxygen concentration in umol/kg (salinity corrected and pressure
% corrected (input))
oxy_conc=molar2molal(oxy_conc0,temp,sal,P_dbar);

% set actual output: 
oxygen=reshape(oxy_conc0,shape);

%figh=figure;
%axh=scatter3(temp,phase,oxygen,'filled');
%set(figh,'color',[1 1 1])
%set(gcf,'Position',[500 350 560 420]) % original size
%title(['foil coefficient evaluation (' modeltype '-model)'])
%view(110,20)

%function out=scaledT(in)
%% calculate scaled temperature
%out=log((298.15-in)./(273.15+in));

%function pw=watervapor(T,S)
%% calculating pH2O / atm after Weiss and Price 1980
%% T in �C
%pw=(exp(24.4543-(67.4509*(100./(T+273.15)))-(4.8489*log(((273.15+T)./100)))-0.000544.*S));

%function out=O2solubility(T,S)
%% calculate oxygen solubilty / umol/l
%sca_T = scaledT(T);
%out=((exp(2.00856+3.224.*sca_T+3.99063.*sca_T.^2+4.80299.*sca_T.^3+0.978188.*sca_T.^4+...
%    1.71069.*sca_T.^5+S.*(-0.00624097-0.00693498.*sca_T-0.00690358.*sca_T.^2-0.00429155.*sca_T.^3)...
%    -0.00000031168.*S.^2))./0.022391903);

%function O2sal=O2freshtosal(O2fresh,T,S)
%% apply salinity correction to oxygen concentration / umol/l
%sca_T = scaledT(T);
%O2sal=O2fresh.*exp(S.*(-0.00624097-0.00693498*sca_T-0.00690358*sca_T.^2-0.00429155*sca_T.^3)-3.11680e-7*S.^2);

%function O2conc_sal=pO2toO2conc(pO2,T,S,P_atm)
%% calculate O2 concentration in umol/l from pO2 / mbar with salinity and
%% atm. pressure / mbar (or GTD pressure, respectively) correction
%
%% pH2O / atm
%pH2O = watervapor(T,S);
%% atm. pressure / atm
%atm_press=P_atm/1013.25; 
%% theta0
%th0=1-(0.999025+0.00001426.*T-0.00000006436.*T.^2);
%% O2 solubility / umol/l
%oxy_sol=O2solubility(T,S);
%% pressure corrected O2 solubility / umol atm / l
%oxy_sol_pc=oxy_sol.*atm_press.*(((1-pH2O./atm_press).*(1-th0.*atm_press))./((1-pH2O).*(1-th0)));
%% oxygen concentration in umol/l (salinity corrected)
%O2conc_sal=(pO2.*oxy_sol_pc)./((atm_press-pH2O).*0.20946.*1013.25);

%function O2sat_sal=pO2toO2sat(pO2,T,S,P_atm)
%% calculate O2 saturation in % from pO2 / mbar with salinity and
%% atmospheric/GTD pressure / mbar correction
%
%pH2O = watervapor(T,S); % pH2O / atm
%atm_press=P_atm/1013.25; % atm. pressure / atm
%th0=1-(0.999025+0.00001426.*T-0.00000006436.*T.^2); % theta0
%oxy_sol=O2solubility(T,S); % O2 solubility / umol/l
%% pressure corrected O2 solubility / umol atm / l
%oxy_sol_pc=oxy_sol.*atm_press.*(((1-pH2O./atm_press).*(1-th0.*atm_press))./((1-pH2O).*(1-th0)));
%% oxygen concentration in umol/l (salinity corrected)
%O2conc_sal=(pO2.*oxy_sol_pc)./((atm_press-pH2O).*0.20946.*1013.25);
%% saturation
%O2sat_sal=((O2conc_sal)./oxy_sol_pc).*100;

%function out_kg=molar2molal(in_l,T,S,P_dbar)
%% convert from molar unit x/l to molal unit x/kg using seawater densitiy at
%% pressure level P / db(pressure level of intake: 11 dbar or surface/uw 
%% box: 0 dbar)
%
%if nargin<4
%	P_dbar=0; %surface/box
%end
%
%dens_ss=sw_pden(S,T,P_dbar,0); % potential density at S, T, P with 0 dbar as reference
%out_kg=in_l./(dens_ss./1000); %umol/l -> umol/kg
