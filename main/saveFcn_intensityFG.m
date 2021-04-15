function saveFcn_intensityFG(subjID, deviantType, expName, saveType, varargin)
%% Function to save out the parameters and results of intensityFG_main
%
% USAGE: saveFcn_intensityFG(subjID, deviantType, expName, varargin)
% 
% The function generates a unique file name and saves out the variables
% passed to it when called.
% The file name is constructed from the experimental function name, the
% task type, the subject ID and the current time.
%
% Mandatory inputs:
% subjID        - Numeric value, one of 1:99, subject ID.
% deviantType   - Char array, one of {'figure', 'background'}. Defines the
%               task.
% expName       - Char array, name of the experimental function.
% saveType      - Char array, one of {'temporary', 'final'}. Appended to
%               the save file name. 
%
% Optional inputs:
% All optional inputs passed via varargin are saved out.
%
%


%% Input checks

if nargin < 4
    error('Input arguments "subjID", "deviantType", "expName" and "saveType" are mandatory!');
end
if ~isnumeric(subjID) || ~ismember(subjID, 1:99)
    error('Input arg "subjID" should be an integer in range 1:99!');
end
if ~ischar(deviantType) || ~ismember(deviantType, {'figure', 'background'})
    error('Input arg "deviantType" should be either "figure" or "background"!');
end
if ~ischar(expName)
    error('Input arg "expName" should be a character array!');
end
if ~ischar(saveType) || ~ismember(saveType, {'temporary', 'final'})
    error('Input arg "saveType" should be one of {"temporary", "final"}!');
end


%% Sort passed variables into fields of a struct for saving out
% IMPORTANT! In the loop, "i+4" depends on the number of mandatory args,
% needs to be changed if arguments are changed

output = struct;
for i = 1:length(varargin)
    output.(inputname(i+4)) = varargin{i};
end


%% Get unique name 

% provide unique filename string to ensure previous data parameters don't
% accidently get written over 
tmp = clock;
hrStr = num2str(tmp(4));
minStr = num2str(tmp(5));

if length(hrStr) < 2; hrStr = ['0', hrStr]; end  % force the hour/minute string to have 4 digits in total
if length(minStr) < 2; minStr = ['0', minStr]; end  % force the hour/minute string to have 4 digits in total

dateStr = [num2str(tmp(1)), '_', num2str(tmp(2)), '_', num2str(tmp(3)), '_', hrStr, minStr];

saveFileName = [expName, '_subj', num2str(subjID), '_', deviantType, '_', saveType, '_', dateStr, '.mat'];


%% Save

save(saveFileName,...
    'subjID', 'deviantType', 'expName', 'saveType',...
    'output'); 


return