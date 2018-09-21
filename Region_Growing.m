function varargout = Region_Growing(varargin);
%
% Syntax :
%     [labels] = Region_Growing(imMatrix);
%
% This function performs a region growing algorithm. Clusters of points
% with label equal to zero will be identified and labeled.
%
% Input Parameters:
%        imMatrix                       : Input matrix(Nx1,NxM or NxMxO).
% Output Parameters:
%         labels                        : Labels
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0


%% ================================ Main Program ====================================== %%
imMatrix = varargin{1};

labels = imMatrix;
ind = find(imMatrix ~= 0);
labels(ind) = -1;
cont = 1;
matDim = size(imMatrix);
dim = length(matDim);
if dim ==2
    if find(matDim == 1)
        dim = 1;
    end
end
switch dim
    case 1 % Region growing algorithm for one dimensional data
        Npoints = length(imMatrix);
        ind = find(imMatrix == 0);
        
        while ~isempty(ind)
            labels(ind(1)) = cont;
            Neigh = [ind(1)-1;ind(1)+1];
            Neigh(Neigh<1|Neigh>Npoints) = [];
            ind2rem = find(labels(Neigh)~=0);
            Neigh(ind2rem) = [];
            if isempty(Neigh)
                cont = cont +1
            else
                labels(Neigh) = cont;
            end
            ind = find(labels == 0);
        end

% % % % % % % % %     
        
    case 2
    case 3
end
ind = find(labels == -1);
labels(ind) = 0;
%% ====================== End of main program  =================00===================== %%
% Outputs
varargout{1} = labels;
return;
