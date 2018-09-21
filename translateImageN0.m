function r = translateImageN0(f,di,dj,dk)

% (subfunction)
% Translation of a ND image in 3D space
%
% Author: Pedro Antonio Valdes-Hernandez
% Neuroimaging Department
% Cuban Neuroscience Center
% March 1 2001
% Version $4.0

%=========================Main program====================================%
S = size(f);
r = zeros(S);
N = S(1); M = S(2); P = S(3);
if di > 0
    iindr =  di+1:N;
    iindf =  1:N-di;
elseif di < 0
    iindr =  1:N+di;
    iindf = -di+1:N;
else
    iindr = 1:N;
    iindf = 1:N;
end
if dj > 0
    jindr =  dj+1:M;
    jindf =  1:M-dj;
elseif dj < 0
    jindr =  1:M+dj;
    jindf = -dj+1:M;
else
    jindr = 1:M;
    jindf = 1:M;
end
if dk > 0
    kindr =  dk+1:P;
    kindf =  1:P-dk;
elseif dk < 0
    kindr =  1:P+dk;
    kindf = -dk+1:P;
else
    kindr = 1:P;
    kindf = 1:P;
end
r(iindr,jindr,kindr,:) = f(iindf,jindf,kindf,:);
%========================End of main program==============================%
return;