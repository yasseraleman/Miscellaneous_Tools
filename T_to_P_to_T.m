function OutMap = T_to_P_to_T(T,convtype,n);
%
% Syntax :
% OutMap = T_to_P_to_T(T,convtype)
%
% This scrips converts T Maps to p-values maps and viceversa
%
%
% Input Parameters:
%         Map                   : P or T Map
%       convtype                : Conversion type(from T to P: 't2p', from P to T: 'p2t')
%         n                     : Degrees of Freedom
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
    case 't2p'
        OutMap=2*(1-tcdf(abs(T),n-2))
    case 'p2t'
        OutMap= tinv(1-T/2,n-2);
end
return;