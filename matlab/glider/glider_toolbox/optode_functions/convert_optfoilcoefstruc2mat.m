function foilcoef = convert_optfoilcoefstruc2mat(C)
%% function
%    foilcoef = convert_optfoilcoefstruc2mat(C)
% 
% GEOMAR SVN $Id: convert_optfoilcoefstruc2mat.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% function to convert format of optode Aanderaa foil coefficients:
%
%       from structure C --> matrix
%
% INPUT:    C ... structure of foil coefficient file
%
% OUTPUT:   foilcoef ... foil coefficients as matrix
%
% version 1.0         last change 22.04.2014
%
% Author:
%   J.Hahn, GEOMAR, Kiel, Apr. 2014
%   G.Krahmann
%
%   initial coding                  J.Hahn      Apr. 2014
%   merging from 2 functions of J.Hahn

if isfield(C,'C0Coeff0')
  foilcoef = [C.C0Coeff0 C.C0Coeff1 C.C0Coeff2 C.C0Coeff3 ...
              C.C1Coeff0 C.C1Coeff1 C.C1Coeff2 C.C1Coeff3 ...
              C.C2Coeff0 C.C2Coeff1 C.C2Coeff2 C.C2Coeff3 ...
              C.C3Coeff0 C.C3Coeff1 C.C3Coeff2 C.C3Coeff3 ...
              C.C4Coeff0 C.C4Coeff1 C.C4Coeff2 C.C4Coeff3 ]';
else
  foilcoef = [C.FoilCoefA(1)  C.FoilCoefA(2)  C.FoilCoefA(3)  C.FoilCoefA(4)  ...
              C.FoilCoefA(5)  C.FoilCoefA(6)  C.FoilCoefA(7)  C.FoilCoefA(8)  ...
              C.FoilCoefA(9)  C.FoilCoefA(10) C.FoilCoefA(11) C.FoilCoefA(12) ...
              C.FoilCoefA(13) C.FoilCoefA(14) C.FoilCoefB(1)  C.FoilCoefB(2)  ...
              C.FoilCoefB(3)  C.FoilCoefB(4)  C.FoilCoefB(5)  C.FoilCoefB(6)  ...
              C.FoilCoefB(7) ...
               ]';
end
