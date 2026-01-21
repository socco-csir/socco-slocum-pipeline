function foilcoef = convert_optfoilcoefstruc2mat_3830(C)
%% function
%    foilcoef = convert_optfoilcoefstruc2mat_3830(C)
% 
% GEOMAR SVN $Id: convert_optfoilcoefstruc2mat_3830.m 186 2016-07-04 14:13:01Z gkrahmann@geomar.de $
%
% function to convert format of optode Aanderaa foil coefficients:
%
%       from structure C --> matrix
%
% INPUT:    C ... structure of foil coefficient file
%
% OUTPUT:   foilcoef ... foil coefficients as matrix
%
% version 1.0         last change 08.04.2014
%
% Author:
%   J.Hahn, GEOMAR, Kiel, Apr. 2014
%
%   initial coding                  J.Hahn      Apr. 2014

foilcoef = [C.C0Coeff0 C.C0Coeff1 C.C0Coeff2 C.C0Coeff3 ...
            C.C1Coeff0 C.C1Coeff1 C.C1Coeff2 C.C1Coeff3 ...
            C.C2Coeff0 C.C2Coeff1 C.C2Coeff2 C.C2Coeff3 ...
            C.C3Coeff0 C.C3Coeff1 C.C3Coeff2 C.C3Coeff3 ...
            C.C4Coeff0 C.C4Coeff1 C.C4Coeff2 C.C4Coeff3 ]';
