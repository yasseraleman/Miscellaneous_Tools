function varargout = cohensD(varargin);
%
% Syntax :
%  [cohensd] = cohensD(X,Y)
%
% This script computes the effect size (Cohen's d ) between variables X and Y. 
%
% Input Parameters:
%       X                     : First variable.
%       Y                     : Second variable.
%
% Output Parameters:
%      cohensd                : Cohen's d value.
%
% See also:
%__________________________________________________
% Authors:Luis Marcos
% LIM, HUGGM
% March 3rd 2017
% Version $1.0

%% ============================ Checking input parameters ============================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    X = varargin{1};
    Y = varargin{2};
    sampType = 'indep';
    if nargin > 3
        error('To many outputs');
        return;
    elseif nargin == 3
         sampType = varargin{3};
    end
end
%% ===================== End of checking input parameters ============================== %

%% ================================== Main Program ===================================== %
if length(X) == length(Y)
    switch sampType
        case 'indep'
            spooled = sqrt((sum((X - repmat(mean(X),[size(X,1) 1])).^2) + sum((Y - repmat(mean(Y),[size(Y,1) 1])).^2))/(size(X,1)+size(Y,1)-2));
            cohensd = (mean(X)-mean(Y))./spooled;
        case 'related'
            cohensd = (mean(X)-mean(Y))./ std(X-Y);
    end
else
    n1 = length(X);
    n2 = length(Y);
    s1 = std(X);
    s2 = std(Y);
    sp = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
    cohensd = (mean(X)-mean(Y))./ sp;
end

%% =========================== End of Main Program ===================================== %
% Outputs
varargout{1} = cohensd;

end
