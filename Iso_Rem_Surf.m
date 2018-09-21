function [I,IB] = Iso_Rem_Surf(T,Nhood);
%
% Syntax :
%     [I,IB] = Iso_Rem_Surf(T,Nhood);
%
% This function removes boundary points from a mask if the number of neightbor 
% points is less than a predefined neighborhood. 
%
% Input Parameters:
%   T            : White Matter  Mask
%   Nhood        : Minimun number of neighbors.
%
% Output Parameters:
%   I            : White Matter Mask without isolated points
%__________________________________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

warning off
ind = find(T);
T = zeros(size(T));
T(ind) = 1;
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
indbt = 0;
for i =1:26
    M = Neib(i,:);
    S = [x y z] + M(ones(s,1),:);
    ind2 = sub2ind(size(I),S(:,1),S(:,2),S(:,3));
    sROI(ind) = sROI(ind) + I(ind2);
    indb = find(I(ind2)==0);
    indt = sub2ind(size(I),x(indb),y(indb),z(indb));
    indbt = [indbt; indt];
end
ind = indbt(2:end,1);
indb = unique(ind);
IB = zeros(size(I));
IB(indb) = 1;
ind = find(sROI<Nhood);
I(ind) = 0;
IB(ind) = 0;
I = I(2:end-1,2:end-1,2:end-1);
IB = IB(2:end-1,2:end-1,2:end-1);
return;