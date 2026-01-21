function foilcoef = convert_optfoilcoefstruc2mat_4330(C)
%% function
%    foilcoef = convert_optfoilcoefstruc2mat_4330(C)
% 
% GEOMAR SVN $Id: convert_optfoilcoefstruc2mat_4330.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% function to convert format of optode Aanderaa foil coefficients:
%
%       from structure C --> matrix
%
% INPUT:    C ... structure of foil coefficient file
%
% OUTPUT:   foilcoef ... foil coefficients as matrix
%
% version 1.0         last change 11.04.2014
%
% Author:
%   J.Hahn, GEOMAR, Kiel, Apr. 2014
%
%   initial coding                  J.Hahn      Apr. 2014

foilcoef = [C.FoilCoefA(1)  C.FoilCoefA(2)  C.FoilCoefA(3)  C.FoilCoefA(4)  ...
            C.FoilCoefA(5)  C.FoilCoefA(6)  C.FoilCoefA(7)  C.FoilCoefA(8)  ...
            C.FoilCoefA(9)  C.FoilCoefA(10) C.FoilCoefA(11) C.FoilCoefA(12) ...
            C.FoilCoefA(13) C.FoilCoefA(14) C.FoilCoefB(1)  C.FoilCoefB(2)  ...
            C.FoilCoefB(3)  C.FoilCoefB(4)  C.FoilCoefB(5)  C.FoilCoefB(6)  ...
            C.FoilCoefB(7) ...
             ]';
