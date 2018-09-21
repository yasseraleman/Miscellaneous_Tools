function I = Iso_Rem(T,Nhood);
%
%This function removes isolated points from mask. 
%
% Input Parameters:
%   T            : Mask
%   Nhood        : Minimun number of neighbors. 
% Output Parameters:
%   I            : Mask without isolated points  
%__________________________________________________________________________
% Authors:  Yasser Alem?n G?mez 
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

%=========================Main program====================================%  
warning off
if length(size(T))==3
    I = zeros(size(T)+2);
    I(2:end-1,2:end-1,2:end-1) = T;
    clear T
    ind = find(I>0);
    [x,y,z] = ind2sub(size(I), ind);
    s = size(x,1);
    sROI = zeros(size(I));
    [X, Y, Z] = meshgrid(-1:1,-1:1,-1:1);
    X = X(:);Y = Y(:);Z = Z(:);
    Neib = [X Y Z];clear X Y Z;
    pos = find((Neib(:,1)==0)&(Neib(:,2)==0)&(Neib(:,3)==0));
    Neib(pos,:) = [];
    for i =1:26
        M = Neib(i,:);
        S = [x y z] + M(ones(s,1),:);
        ind2 = sub2ind(size(I),S(:,1),S(:,2),S(:,3));
        sROI(ind) = sROI(ind) + I(ind2);
    end
    ind = find(sROI<Nhood);
    I(ind) =0;
    I = I(2:end-1,2:end-1,2:end-1);
elseif length(size(T))==2
    I = zeros(size(T)+2);
    I(2:end-1,2:end-1) = T;
    clear T
    ind = find(I>0);
    [x,y] = ind2sub(size(I), ind);
    s = size(x,1);
    sROI = zeros(size(I));
    [X, Y] = meshgrid(-1:1,-1:1);
    X = X(:);Y = Y(:);
    Neib = [X Y];clear X Y;
    pos = find((Neib(:,1)==0)&(Neib(:,2)==0));
    Neib(pos,:) = [];
    for i =1:8
        M = Neib(i,:);
        S = [x y] + M(ones(s,1),:);
        ind2 = sub2ind(size(I),S(:,1),S(:,2));
        sROI(ind) = sROI(ind) + I(ind2);
    end
    ind = find(sROI<Nhood);
    I(ind) =0;
    I = I(2:end-1,2:end-1);
end

%========================End of main program==============================%
return;