function varargout = Branch_Labelling(varargin);
%
% Syntax :
%      [reordsLines] = Branch_Labelling(Surf, sLines);
%
% This function computes all sulcal lines using a curvature map.
%
% Input Parameters:
%        Surf                    : Matlab surface variable. The surface
%                                  structure variable must containg the
%                                  Is field with curvature values.
%        sLines                  : Sulcal Lines Branches
%
% Output Parameters:
%      reordsLines               : Lines with labeled branches.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
% % % % %  load('/media/COSAS/scripts/Gyral_Crowns_and_Sulcal_Lines_Extraction/matlab.mat');
if nargin < 2
    error('Two Inputs are mandatory');
    return
end
Surf = varargin{1};
sLines = varargin{2};

% Surface Checking
Surf = Surface_Checking(Surf);

if nargin > 2
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end

%% ========================= End of Checking Inputs ==================== %%

%% ======================= Main Program ================================ %%

% ======================= Removing Isolated clusters ================== %%
tempLines = sLines(:,1:2);  % Keeping only the selected edges
[~,~,auxOrder] = unique(tempLines(:));   % Locating the order in the original edge matrix
Nedge = size(tempLines,1);               % Number o edges
reorgtempLines = reshape(auxOrder,[Nedge 2]);  % Reorganicing edges to avoid memory problems in the graph creating step
reorgtempLines = sort(reorgtempLines')';

Npoints = max(reorgtempLines(:));
graphMat = sparse(Npoints,Npoints); % Creating empty Graph Matrix
indg = sub2ind(size(graphMat),[reorgtempLines(:,1) reorgtempLines(:,2)],[reorgtempLines(:,2) reorgtempLines(:,1)]); % Establishing connections
graphMat(indg) = 1;
graphMat = Label_Graph_Components(graphMat); % Labelling Networks

if max(graphMat(:)) > 1
    temp = nonzeros(graphMat(:));
    temp2 = accumarray(temp, temp*0+1);
    [~,maxgraph] = max(temp2);
    [X,Y]  = find(graphMat == maxgraph); % Locating pairs of connected points
    auxPairs = unique(sort([X Y]')','rows'); % Auxiliar organization
    indr = ismember(reorgtempLines, auxPairs,'rows'); % Indexes to reorganice
    tempLines  = tempLines(indr,:); % Each network edge
    ind2rem = find(ismember(sort(sLines(:,1:2)')',sort(tempLines')','rows') == 0); % Edge To remove
    sLines(ind2rem,:) = [];
end
% ================= End of  Removing Isolated clusters ================ %%

tempVar = sLines(:,1:2);
tempVar = accumarray(tempVar(:),tempVar(:)*0+1);
inde  = find(tempVar == 1); % Endpoints
indb  = find(tempVar > 2);  % Junction points

MaxCurvPath = Compute_Max_Curvature_Path(Surf, sLines(:,1:2));
sLines(:,3) = ismember(sort(sLines(:,1:2)')',sort(MaxCurvPath')','rows');

% ====================== Reorganicing Main Line ======================== %
indMainline = find(sLines(:,3) == 1);
reordsLines = sLines(indMainline,1:2);  % Keeping only the edges from main line
tempVar = accumarray(reordsLines(:),reordsLines(:)*0+1);
endpMainline  = find(tempVar == 1);
start_point = endpMainline(1);
end_point = endpMainline(2);

% Building Graph
Graph = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices)); % Empty Graph
ind = sub2ind(size(Graph),[sLines(:,1);sLines(:,2)],[sLines(:,2);sLines(:,1)]); % Branch index
Graph(ind) = 1; % Curvature weight (the inverse of the sum of curvature).

% Endpoint - junction points pairs
distances = graphallshortestpaths(Graph);

if ~isempty(indb)
    D = distances(inde,indb); % Distance between endpoints and branches
    [minVal,location] = min(D,[],2); % Finding endpoint-branch pair
    end_points = inde;
    start_points = indb(location);
    ind2rem = ismember(inde,endpMainline);
    start_points(ind2rem) = [];
    end_points(ind2rem) = [];
end



[DIST, nativePath]=graphshortestpath(Graph,start_point,end_point);
nativePath = nativePath(:);

% Remove points with high curvature
nativePath = Remove_High_Curvature_Points_from_Path(Surf, nativePath, indb);

reordsLines = [nativePath(1:end-1) nativePath(2:end) ones(length(nativePath)-1,1)];

indNomainLine = find(sLines(:,3)~=1);

% =================== End of Reorganicing Main Line ==================== %


sLinesTemp = sLines(indNomainLine,:); %removing main line

if ~isempty(sLinesTemp)
    Graphb = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices)); % Empty Graph
    indb = sub2ind(size(Graphb),[sLinesTemp(:,1);sLinesTemp(:,2)],[sLinesTemp(:,2);sLinesTemp(:,1)]); % Branch index
    
    %     curvGraph = surface2graph(Surf,Surf.Is);
    
    % ============= Creating Distance and Curvature Graphs ============ %%
    % -------- Create Empty Graphs
    distGraph = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices)); % Empty Distance Graph
    curvGraph = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices)); % Empty Curvature Graph
    indall = sub2ind(size(curvGraph),[sLines(:,1);sLines(:,2)],[sLines(:,2);sLines(:,1)]); % Edges indexes
    
    % -------- Create Distance Graph
    distGraph = surface2graph(Surf);
    X = Surf.SurfData.vertices([sLines(:,1);sLines(:,2)],1) - Surf.SurfData.vertices([sLines(:,2);sLines(:,1)],1);
    Y = Surf.SurfData.vertices([sLines(:,1);sLines(:,2)],2) - Surf.SurfData.vertices([sLines(:,2);sLines(:,1)],2);
    Z = Surf.SurfData.vertices([sLines(:,1);sLines(:,2)],3) - Surf.SurfData.vertices([sLines(:,2);sLines(:,1)],3);
    distGraph(indall) = sqrt(X.^2 + Y.^2 +Z.^2);
    
    
    % -------  Create Curvature Graph
    curvGraph(indall) = (Surf.Is([sLines(:,1);sLines(:,2)]) + Surf.Is([sLines(:,1);sLines(:,2)]))/2;
    % ======= End of Creating Distance and Curvature Graphs =========== %%
    Nc = length(end_points);
    % ===================== Creating Branch Paths ====================== %
    cont = 0;
    for i  = 1:Nc  % For each cluster
        
        start_point =  start_points(i);
        end_point   =  end_points(i);
        
        
        
        [DISTcurv, nativePathcurv]=graphshortestpath(curvGraph,start_point,end_point); nativePathcurv = nativePathcurv(:);
        [DIST,     ~]=graphshortestpath(distGraph,start_point,end_point); nativePath = nativePath(:);
        
        for j = 1:length(DIST)
            cont = cont + 1;
            if iscell(nativePathcurv)
                branchPath = nativePathcurv{j};
            else
                branchPath = nativePathcurv;
            end
            branchPath = branchPath(:);
            branchPath = Remove_High_Curvature_Points_from_Path(Surf, branchPath);
            
            BranchPaths{cont} = [branchPath(1:end-1) branchPath(2:end) ones(length(branchPath)-1,1)*(cont+1)];
            branchLength(cont) = DIST(j);
        end
    end
    % =================== End of Creating Branch Paths ================== %
    
    
    if cont
        [branchLength,tempOrd] = sort(branchLength,'descend');
        reordsLines = [reordsLines;vertcat(BranchPaths{tempOrd})];
        
        Nedges = accumarray(reordsLines(:,3),reordsLines(:,3)*0+1);
        edge2del = ismember(reordsLines(:,3),find(Nedges <=3));
        reordsLines(edge2del,:) = [];
        
    end
    Nedges = accumarray(reordsLines(:,3),reordsLines(:,3)*0+1);
    [~,tempOrd] = sort(Nedges,'descend');
    newreordsLines = [ 0 0 0];
    for j = 1:length(Nedges)
        ind = find(reordsLines(:,3) == tempOrd(j));
        newreordsLines = [newreordsLines;reordsLines(ind,1:2) ones(length(ind),1)*j];
    end
    newreordsLines(1,:) = [];
    reordsLines = newreordsLines;
end



%% ======================== End of Main Program ========================= %
% Outputs
varargout{1} = reordsLines;


% % % % % % % %% %%%%%%%%%% Testing THINGS
% % % % % % % 
% % % % % % % col = [1 1 1; 0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
% % % % % % % Ncolor = size(col,1);
% % % % % % % re = floor(Nc/Ncolor); col = repmat(col,[re+1 1]);
% % % % % % % 
% % % % % % % 
% % % % % % % % Compute Neighborhoods
% % % % % % % [Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
% % % % % % % Temp = sum(Trip);
% % % % % % % Trip(:,Temp==0) = [];
% % % % % % % temp = Trip(:,3:end);
% % % % % % % indz = find(temp == 0);
% % % % % % % temp(indz) = 1;
% % % % % % % 
% % % % % % % % Compute Normals
% % % % % % % N = patchnormals(Surf.SurfData);
% % % % % % % X = N(temp,1);X(indz) = 0;X = sum(X,2)./sum(logical(X),2);
% % % % % % % Y = N(temp,2);Y(indz) = 0;Y = sum(Y,2)./sum(logical(Y),2);
% % % % % % % Z = N(temp,2);Z(indz) = 0;Z = sum(Z,2)./sum(logical(Z),2);
% % % % % % % N = [X Y Z];
% % % % % % % norma = normm(N);
% % % % % % % N = N./[norma norma norma];
% % % % % % % hfig = figure('numbertitle','off','name','New Figure','Color',[1 1 1],'InvertHardcopy','off');
% % % % % % % 
% % % % % % % branchlabels = unique(reordsLines(:,3));
% % % % % % % nLabels = length(branchlabels);
% % % % % % % for i = 1:nLabels
% % % % % % %     indbranch = find(reordsLines(:,3) == branchlabels(i));
% % % % % % %     nativePath = [reordsLines(indbranch,1);reordsLines(indbranch(end),2)];
% % % % % % %     Y = abs(acos(dot(N(nativePath,:),repmat(N(start_point,:),[length(nativePath) 1]),2))*180/pi);X = [1:length(Y)];X = X(:);Y = Y(:);
% % % % % % %     % Detecting peaks
% % % % % % %     if length(Y)>6
% % % % % % %         [~,localiz] = findpeaks(Y);
% % % % % % %         indcorrect = find([Y(localiz-1)-Y(localiz+1)] <0);
% % % % % % %         Y(localiz(indcorrect)) = [Y(localiz(indcorrect)-1)+Y(localiz(indcorrect)+1)]/2;
% % % % % % %         FO = fit(X, Y,'smoothingspline');
% % % % % % % 
% % % % % % %         X = [1:length(Y)/(length(Y)*5):length(Y)]';
% % % % % % %         Y = FO(X);
% % % % % % %         subplot(2,2,1,'Color',[0 0 0])
% % % % % % %         hold on
% % % % % % %         plot(X,Y,'Color',col(i,:));
% % % % % % %     else
% % % % % % %         subplot(2,2,1,'Color',[0 0 0]);
% % % % % % %         hold on
% % % % % % %         plot(Y,'Color',col(i,:));
% % % % % % %     end
% % % % % % % 
% % % % % % %     subplot(2,2,2,'Color',[0 0 0]);
% % % % % % %     hold on
% % % % % % %     plot(Surf.Is(nativePath),'Color',col(i,:));
% % % % % % %     cont(i,1) = sum(Surf.Is(nativePath));
% % % % % % % 
% % % % % % %     if i == 1
% % % % % % %         subplot(2,2,3,'Color',[0 0 0]);
% % % % % % %         Plot_Surf(Surf,'FigID',hfig);
% % % % % % %         Xp = [Surf.SurfData.vertices(reordsLines(:,1),1) Surf.SurfData.vertices(reordsLines(:,2),1)]';
% % % % % % %         Yp = [Surf.SurfData.vertices(reordsLines(:,1),2) Surf.SurfData.vertices(reordsLines(:,2),2)]';
% % % % % % %         Zp = [Surf.SurfData.vertices(reordsLines(:,1),3) Surf.SurfData.vertices(reordsLines(:,2),3)]';
% % % % % % %         subplot(2,2,3);
% % % % % % %         line(Xp,Yp,Zp,'Color',[1 1 1],'Linewidth',2);
% % % % % % %     end
% % % % % % % 
% % % % % % %     subplot(2,2,4,'Color',[0 0 0]);
% % % % % % %     Plot_Surf(Surf,'FigID',hfig);
% % % % % % %     Xp = [Surf.SurfData.vertices(reordsLines(indbranch,1),1) Surf.SurfData.vertices(reordsLines(indbranch,2),1)]';
% % % % % % %     Yp = [Surf.SurfData.vertices(reordsLines(indbranch,1),2) Surf.SurfData.vertices(reordsLines(indbranch,2),2)]';
% % % % % % %     Zp = [Surf.SurfData.vertices(reordsLines(indbranch,1),3) Surf.SurfData.vertices(reordsLines(indbranch,2),3)]';
% % % % % % %     line(Xp,Yp,Zp,'Color',col(i,:),'Linewidth',2);
% % % % % % %     hlink = linkprop([subplot(2,2,3);subplot(2,2,4)], {'CameraPosition','CameraUpVector'});
% % % % % % %     a = 1;
% % % % % % % end

return;