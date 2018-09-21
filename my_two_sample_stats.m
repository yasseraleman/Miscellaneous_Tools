function varargout = my_two_sample_stats(X,Y,varargin)
%
% Syntax :
%  [statsCad] = my_two_sample_stats(X,Y,varargin)
%
% This script computes different statistical tests between two samples
%
% Input Parameters:
%       X                    : Sample 1
%       Y                    : Sample 2
%
%
% Output Parameters:
%      statsCad               : tats from comparison between data columns.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 22th 2014
% Version $1.0

if nargin ~=0
    method = 'correlation';
    statparam = 'parametric';
    cohensMode = 'indep';
    tailvar = 'both';
end


% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'method'
                    method = varargin{2}; % Methods (correlation,hipotest,cohensd,regression)
                case 'statparam' % Type of analysis (parametric or nonparametric)
                    statparam=varargin{2};
                case 'cohensmode'
                    cohensMode = varargin{2};
                case 'tailvar' % Tails (both, left or right)
                    tailvar=varargin{2};
            
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

%% ================================ Removing Outliers ================================== %
% % % % plotOp = 0;
% % % % noutliers = 2;
% % % % [X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
% % % % X(outliers_idx) = [];
% % % % Y(outliers_idx) = [];
%% ========================= End of Removing Outliers ================================== %

%% ================================ Main Program ======================================= %

switch lower(method)
    case 'correlation'
        switch lower(statparam)
            case 'parametric'
                switch lower(tailvar)
                    case 'right'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','right');
                    case 'left'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','left');
                    case 'both'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','both');
                end
                statsCad = [strvcat('Pearson Correlation Value','pValue') repmat(' , ',[2 1]) num2str([rcorr;rpval])];

            case 'nonparametric'
                switch lower(tailvar)
                    case 'right'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','right','Type','Spearman');
                    case 'left'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','left','Type','Spearman');
                    case 'both'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','both','Type','Spearman');
                end
                        statsCad = [strvcat('Spearman Correlation Value','pValue') repmat(' , ',[2 1]) num2str([rcorr;rpval])];
       end
        
    case 'hipotest'
        switch lower(statparam)
            case 'parametric'
                switch lower(tailvar)
                    case 'right'
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','right');
                    case 'left'
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','left');
                    case 'both'
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','both');
                end
                tValues   = STATS.tstat;
                ttestpvalValues = tpval;
                statsCad = [strvcat('T-Value','pValue') repmat(' , ',[2 1]) num2str([STATS.tstat;tpval])];
            case 'nonparametric'
                switch lower(tailvar)
                    case 'right'
                        [upval,H,STATS] = ranksum(X(:),Y(:),'tail','right');
                    case 'left'
                        [upval,H,STATS] = ranksum(X(:),Y(:),'tail','left');
                    case 'both'
                        [upval,H,STATS] =ranksum(X(:),Y(:),'tail','both');
                end
                statsCad = [strvcat('U-ZValue','pValue') repmat(' , ',[2 1]) num2str([STATS.zval;upval])];
        end
    case 'cohensd'
        switch cohensMode
            case 'related'
                cohens_D = cohensD(X,Y,'related');
            case 'indep'
                cohens_D = cohensD(X,Y,'indep');
        end
       statsCad = ['Cohens D , ' num2str(cohens_D)];

    case 'regression'
        % Regression
        [B,BINT,R,RINT,STATS] = regress(Y,[ones(length(X),1) X]);
        statsCad = [strvcat('Regression Slope','Regression Intercept','R-Square Value') repmat(' , ',[3 1]) num2str([B(2);B(1);STATS(1)])];
    case 'allstats'
        statsCad = '';
        switch lower(statparam)
            case 'parametric'
                switch lower(tailvar)
                    case 'right'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','right');
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','right');
                    case 'left'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','left');
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','left');
                    case 'both'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','both');
                        [H,tpval,CI,STATS] = ttest2(X(:),Y(:),'tail','both');
                end
                statsCad = [strvcat(statsCad,'Pearson Correlation Value','corrpValue','T-Value','ttestpValue')];
                valValues = [rcorr;rpval;STATS.tstat;tpval];
                
            case 'nonparametric'
                switch lower(tailvar)
                    case 'right'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','right','Type','Spearman');
                        [upval,H,STATS] = ranksum(X(:),Y(:),'tail','right');
                    case 'left'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','left','Type','Spearman');
                        [upval,H,STATS] = ranksum(X(:),Y(:),'tail','left');
                    case 'both'
                        [rcorr,rpval] = corr(X(:),Y(:),'tail','both','Type','Spearman');
                        [upval,H,STATS] =ranksum(X(:),Y(:),'tail','both');
                end
                statsCad = [strvcat(statsCad,'Spearman Correlation Value','corrpValue','STATS.zval','upValue')];
                valValues = [rcorr;rpval;STATS.zval;upval];
        end
        [B,BINT,R,RINT,STATS] = regress(Y,[ones(length(X),1) X]);
        statsCad = [strvcat(statsCad,'Regression Slope','Regression Intercept','R-Square Value')];
        valValues = [valValues;B(2);B(1);STATS(1)];
        
        switch cohensMode
            case 'related'
                cohens_D = cohensD(X,Y,'related');
            case 'indep'
                cohens_D = cohensD(X,Y,'indep');
        end
       statsCad = [strvcat(statsCad,'Cohens D')];
       valValues = [valValues;cohens_D];
       statsCad = [statsCad repmat(' , ',[length(valValues) 1]) num2str(valValues)];
end
%% ========================== End of Main Program ====================================== %
% Outputs
varargout{1} = statsCad;
return;