function varargout = Points_Rotation(Coord,varargin)
%
% Syntax :
%     Surfout = Points_Rotation(Coord,'angX', angX,'angY', angY, 'angZ', angZ );
%
% This function rotates points coordinates. 
%
% Input Parameters:
%        Coord                   : Points Coordinates
%
% Output Parameters:
%        rCoord                  : Rotated coordinates
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2011
% Version $1.0

%% ============================= Checking Inputs ======================= %%

if nargin<1 % the indispensable input arguments are not provided
    error('At least one input is mandatory');
    return;
else
    angX = 0;
    angY = 0;
    angZ = 0;
end


% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}

                case 'angX'
                    % Limits in X axis (Left-Right axis)
                    angX=varargin{2};
                case 'angY'
                    % Limits in X axis (Anterior-Posterior axis)
                    angY=varargin{2};
                case 'angZ'
                    angZ=varargin{2};
                                    case 'rotMat'
                     rotMat=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

if nargout > 1
    error('To Many Output Parameters');
    return;
end


%% ========================= End of Checking Inputs ==================== %%

%% ======================== Main Program =============================== %%

angX = angX*pi/180;
angY = angY*pi/180;
angZ = angZ*pi/180;
Rx = [1 0 0 0; 0 cos(angX) -sin(angX) 0; 0 sin(angX) cos(angX) 0; 0 0 0 1];
Ry = [cos(angY) 0 sin(angY) 0; 0 1 0 0; -sin(angY) 0 cos(angY) 0; 0 0 0 1];
Rz = [cos(angZ) -sin(angZ) 0 0; sin(angZ) cos(angZ) 0 0; 0 0 1 0; 0 0 0 1];
Mat = Rx*Ry*Rz;
t = [Coord ones(size(Coord,1),1)]*Mat';
rCoord = t(:,1:3);
%% ======================== End of Main Program ======================== %%

% Outputs
varargout{1} = rCoord;

return;