function params = params_intensityFG_intro
%% Parameters for intensityFG_intro
%
% USAGE: params = params_intensityFG_intro
%
% To be used with intensityFG_intro
%


params = struct;

params.fs = 44100;
params.plotMe = 0;
params.dbscl = -45;
params.numNotesPerToken = 17;
params.lenNote = 0.120; 
params.lenRamp = 0.010;  % length (in seconds) of onset/offset ramps on each chord

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