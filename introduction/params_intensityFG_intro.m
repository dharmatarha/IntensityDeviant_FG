function params = params_intensityFG_intro
%% Parameters for intensityFG_intro
%
% USAGE: params = params_intensityFG_intro
%
% To be used with intensityFG_intro
%
% Output: 
% params    - Struct with its fields defining the stimulus and
%           experimental procedure parameters


%% Init output var

params = struct;


%% Parameters for overall experiment

params.fs = 44100;  % sampling rate
params.plotMe = 0;  % flag for plotting generated stimulus
params.dbscl = -45;  % dB scaling of output signal


%% Parameters for stimulus generation


params.numNotesPerToken = 17;
params.lenNote = 0.120;  % length (in secs) of one note (chord)
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


%% Parameters for triggers

% Params for triggers, using the serial port
% (with Micromed SD LTM EEG in mind)
params.serial.portName = 'COM4';  % For the stimulus PC at SOTE, Szalardy lab
params.serial.baudRate = 9600;  % safe value, should be supported by everything

% Trigger types
params.trig.format = '%i';  % triggers are written to the serial port as integers
params.trig.blockStart = 50;  % at the very start of each block
params.trig.trialStart = 100;  % at the start of each trial
params.trig.soundOnset = 81; 
params.trig.response = 82;
params.trig.hit = 83;
params.trig.falseAlarm = 84;


return