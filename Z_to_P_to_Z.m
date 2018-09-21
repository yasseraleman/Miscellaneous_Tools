function OutMap = Z_to_P_to_Z(Map,convtype)
%
% Syntax :
% Res = Z_to_P_to_Z(Map,convtype)
%
% This scrips converts Z Maps to p-values maps and viceversa
%
%
% Input Parameters:
%         Map                   : P or Z Map
%       convtype                : Conversion type(from Z to P: 'z2p', from P to Z: 'p2z')
%
% Output Parameters:
%      OutMap                   : Output Map
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0

switch convtype
    case 'z2p'
        OutMap = 1-normcdf(abs(Map),0,1);
    case 'p2z'
        Z = @(p) -sqrt(2) * erfcinv(p*2);
        OutMap = Z(Map);
end
return;