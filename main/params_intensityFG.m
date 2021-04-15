function params = params_intensityFG
%% Parameters for intensityFG_main
%
% USAGE: params_intensityFG
%
% To be used with intensityFG_main
%

params = struct;


%% Parameters for overall experiment

params.fs = 44100;  % sampling rate
params.N = 5;  % Number of each of the four stimulus types per block
params.dbscl = -45;  % dB scaling of output signal
params.trialNo = params.N*4;  % number of trials per block
params.blockNo = 6;  % Number of blocks:
params.iti = 1.8;  % set ITI/ISI in secs
% params.iti = 0.760;  % default value from Darrin & Brigi's experiment

params.detectKey = 'space';  % response key
params.goKey = 'Return';  % start key
params.abortKey = 'escape';  % abort key

%--------------------------------------------------------------------------
%%% These are the STIMULUS TYPE LABELS
%--------------------------------------------------------------------------
params.figStandardLabel = 10;  % Figure in FG, no deviant
params.backStandardLabel = 20;  % No figure in FG, no deviant
params.figDeviantLabel = 30;  % Figure in FG, deviant figure
params.backDeviantLabel = 40;  % No figure in FG, deviant background


%% Parameters for stimulus generation

params.numNotesPerToken = 17;
params.lenNote = 0.120;  % length (in secs) of one note (chord)
params.lenRamp = 0.010;  % length (in secs) of onset/offset ramps on each chord
params.stimLength = params.numNotesPerToken*params.lenNote;  % stimulus length in secs

params.clFig = 8;  % This defines the number of DEVIANTS in both the figure AND background deviant conditions
params.clNotFig = 2;

params.clBack = 0; 
params.clNotBack = 10; 

params.clITD = 0*ones(1, params.numNotesPerToken); 
params.clNotITD = 0*ones(1, params.numNotesPerToken);    

params.clTraj = 0*ones(1, params.numNotesPerToken-1);  % this can also be empty [] to generate a random trajectory
params.clStart = []; % use [] to specify random starting notes
params.baseFreq = 150;
params.upFreq = 5000;
params.maxSemitoneStep = 2;
params.minSemitoneStep = 0.5;

params.deviantIndMin = 7;  % Earliest chord that the deviant can begin
params.deviantIndMax = 11;  % Last chord that the deviant can begin
params.deviantDuration = 6;  % Integer duration in # of chords for the deviant

params.deviantSNR = -25;
params.NONdeviantSNR = 0;  % SNR for the non-deviant stimuli

% For background deviant stimuli
params.deviantNum = params.clFig;  % Number of background tones attributed to the background deviant

% For figure deviant stimuli
params.baseSNR = nan;  % when 'nan', the function calculates an estimate based on the sum(cl)/sum(clNot) ratio


return