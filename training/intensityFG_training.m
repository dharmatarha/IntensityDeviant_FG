function intensityFG_training(subjID, trainingType, varargin)
%% Intensity-deviant-detection Figure-ground experiment, training
%
% USAGE: intensityFG_training(subjID, trainingType, feedback='feedback', serialTrigger='trigger') 
% 
% Function for the training part of the Intensity-deviant-detection
% Figure-ground experiment. 
%
% The task is to detect an intensity deviant either in the "figure" or the
% "background" portion of the SFG stimuli. There are four types of trials 
% in the overall experiment:
%   - figure and no deviant; 
%   - figure and deviant figure; 
%   - no figure and no deviant; 
%   - no figure and deviant background.
% 
% During training, the subject is asked to detect deviants in one of the
% following scenarios:
%   (1) Detect deviant figure among figure stimuli (no stimuli with only background)
%   (2) Detect deviant background among background stimuli (no stimuli with figure)
%   (3) Detect deviant figure among figure and background
%   stimuli (all types of stimuli present)
%   (4) Detect deviant background among figure and background
%   stimuli (all types of stimuli present)
%
% Importantly, there is feedback after each response (by default, but can
% be turned off by setting "feedback" arg to "nofeedback").
%
% Main parameters are defined in "params_intensityFG_training".
%
% Mandatory inputs:
% subjID            - Numeric value, one of 1:99. Subject ID.
% trainingType      - Numeric value, one of 1:4. Determines the task: 
%                   1 = deviant figure among figure stimuli;
%                   2 = deviant background among background stimuli;
%                   3 = deviant figure among figure and background stimuli;
%                   4 = deviant background among figure and background stimuli;
%
% Optional inputs:
% feedback          - Char array, one of {'feedback', 'nofeedback'}.
%                   Defaults to 'feedback'. Determines if there is an
%                   auditory feedback after each response or not.
%                   Auditory feedback is loaded from files specified in
%                   params_intensityFG_training. Must be .wav files with
%                   sampling rate matching the sampling rate of the
%                   stimuli.
% serialTrigger     - Char array, one of {'notrigger', 'trigger'}. Controls
%                   if the function sends triggers to the serial port 
%                   for EEG recordings or not. Serial port usage is 
%                   determined by "params.serial" and "params.trig", 
%                   from the output of "params_intensityFG.m". 
%                   Defaults to 'trigger'.
%
% Outputs:
% All outputs are saved out into a subject- and run-specific .mat file. See
% "saveFcn_intensityFG_training" for details of filename generation.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Based mostly on previous code written by:
% - Darrin K. Reed (currently a Senior Research Audiologist at Apple)
% Credit also to:
% - Brigitta Tóth, at RCNS, Budapest
% - Tamas Kurics, Zsuzsanna Kocsis, Botond Hajdu (all ex-RCNS)
%
% Regarding the experiment, please contact the PI:
% Orsolya Szalardy, szalardy.orsolya@med.semmelweis-univ.hu
%
% Written / edited from pieces by Adam Boncz
% 2021.04.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check number of inputs
if ~ismember(nargin, 2:4)    
    error(['Wrong number of input args! Input args "subjID" and ',...
        '"trainingType" are mandatory while input args "feedback" and ',...
        '"serialTrigger" are optional!']);
end

% check mandatory inputs
if ~isnumeric(subjID) || ~ismember(subjID, 1:99)
    error('Input arg "subjID" should be an integer in range 1:99!');
end
if ~isnumeric(trainingType) || ~ismember(trainingType, 1:4)
    error('Input arg "trainingType" should be a numeric value in range 1:4!');
end

% check optional inputs
if ~isempty(varargin)
    for v = 1:numel(varargin)
        if ischar(varargin{v}) && ismember(varargin{v}, {'feedback', 'nofeedback'}) && ~exist('feedback', 'var')
            feedback = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'trigger', 'notrigger'}) && ~exist('serialTrigger', 'var')
            serialTrigger = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "feedback" or "serialTrigger"!');
        end
    end  % for v
end  % if ~isempty
        
% assign default values
if ~exist('feedback', 'var')
    feedback = 'feedback';
end
if ~exist('serialTrigger', 'var')
    serialTrigger = 'trigger';
end

% get logical flag and text for user message from "serialTrigger"
if strcmp(serialTrigger, 'trigger')
    triggerFlag = true;
    triggerText = 'Yes';
elseif strcmp(serialTrigger, 'notrigger')
    triggerFlag = false;
    triggerText = 'No';
end

% get logical flag for "feedback"
if strcmp(feedback, 'feedback')
    feedbackFLAG = 1;
elseif strcmp(feedback, 'nofeedback')
    feedbackFLAG = 0;
end

% user message
disp([char(10), 'Started the Training portion of the Intensity-deviant-detection Figure-ground experiment with inputs:',...
    char(10), 'Subject ID: ', num2str(subjID),...
    char(10), 'Training type: ', num2str(trainingType),...
    char(10), 'Feedback after each response: ', feedback,...
    char(10), 'Triggers on serial port: ', triggerText]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic settings, parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for psychtoolbox installation + unify key names
PsychDefaultSetup(1);
Priority(1);

% load parameters
params = params_intensityFG_training;

% init serial port if triggering was requested
if triggerFlag
    serialObj = serial(params.serial.portName);
    fopen(serialObj);
end

% Set response key code (response key is set in "params" struct
detectKeyCode = KbName(params.detectKey);
% Same for block / task start key and the abort key
goKeyCode = KbName(params.goKey);
abortKeyCode = KbName(params.abortKey);
% if there are multiple keys with the name (e.g. 'Return'), only use the
% first one
if numel(detectKeyCode) ~= 1; detectKeyCode = detectKeyCode(1); end
if numel(goKeyCode) ~= 1; goKeyCode = goKeyCode(1); end
if numel(abortKeyCode) ~= 1; abortKeyCode = abortKeyCode(1); end

% Define a flag for the auditory deviant type: 0 = deviant background, 
% 1 = deviant figure 
if ismember(trainingType, [2 4])
    auditoryFigDeviantFLAG = 0;
    deviantType = 'background';  % used in user messages at the start of task
elseif ismember(trainingType, [1 3])
    auditoryFigDeviantFLAG = 1;
    deviantType = 'figure';
end
% Define a flag for FG stimuli type: 0 = figure only, 1 = background
% only, 2 = both
if trainingType == 1
    FGtypeFLAG = 0;
    FGtype = 'only figure';  % used for user messages at the start of task
elseif trainingType == 2
    FGtypeFLAG = 1;
    FGtype = 'only background';
elseif ismember(trainingType, [3 4])
    FGtypeFLAG = 2;    
    FGtype = 'figure and background';
end

% If there is auditory feedback, prepare the sounds
if feedbackFLAG
    % user message
    disp([char(10), 'There is auditory feedback in this training block, loading and checking feedback sounds...']);
    % load sounds
    [feedbackCorrect, fsTmp1] = audioread(params.feedbackCorrectFile);
    [feedbackWrong, fsTmp2] = audioread(params.feedbackWrongFile);
    % check samplings and length
    if fsTmp1 ~= fsTmp2 || fsTmp1 ~= params.fs
        error('Sampling rates of the feedback sounds and the stimuli are incompatible!');
    end
    if size(feedbackCorrect, 1) ~= size(feedbackWrong, 1)
        error('Feedback sounds have different lengths!');
    end
    % get sound length in secs
    feedbackSoundLength = size(feedbackCorrect, 1)/fsTmp1;
    % user message
    disp(['Feedback sounds are OK. ',...
        'Their sampling rate is ', num2str(fsTmp1), ' Hz and they are ',... 
        num2str(feedbackSoundLength), ' secs long.']);
end

% restrict keys for KbCheck
RestrictKeysForKbCheck([detectKeyCode, goKeyCode, abortKeyCode]);

% Force costly mex functions into memory to avoid latency later on
GetSecs; WaitSecs(0.1); KbCheck();

% Store the name of the current experimental function, used later when
% results are saved out
expName = mfilename(pwd);

% user message
disp([char(10), 'Basic parameters set / loaded', char(10)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init PsychPortAudio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the sound driver - hopefully WASAPI if on Win 10
InitializePsychSound;

% % This part is for selecting a very specific sound device if the default
% % does not work out-of-the-box
% % get correct audio device
device = [];  % system default is our default as well
% targetString = 'ESI Juli@: ICE1724 (hw:2,0)';  % target device name
% % we only change audio device in the lab, when we see the correct audio
% % card
% tmpDevices = PsychPortAudio('GetDevices');
% for i = 1:numel(tmpDevices)
%     if strcmp(tmpDevices(i).DeviceName, targetString)
%         device = tmpDevices(i).DeviceIndex;
%     end
% end

% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;

% open PsychPortAudio device for playback
pahandle = PsychPortAudio('Open', device, mode, reqLatencyClass, params.fs, nrChannels);

% get and display device status
pahandleStatus = PsychPortAudio('GetStatus', pahandle);
disp([char(10), 'PsychPortAudio device status: ']);
disp(pahandleStatus);

% initial start & stop of audio device to avoid potential initial latencies
tmpSound = zeros(2, params.fs/10);  % silence
tmpBuffer = PsychPortAudio('CreateBuffer', [], tmpSound);  % create buffer
PsychPortAudio('FillBuffer', pahandle, tmpBuffer);  % fill the buffer of audio device with silence
PsychPortAudio('Start', pahandle, 1);  % start immediately
PsychPortAudio('Stop', pahandle, 1);  % stop when playback is over

% scale and buffer response sounds if there is feedback after responses
if feedbackFLAG
    % scale the output signal.  First normalize to RMS == 1 and then scale
    % by 'params.dbscl'.
    feedbackCorrect = feedbackCorrect ./ sqrt(mean(feedbackCorrect.^2, 1));  % broadcasting
    feedbackCorrect = feedbackCorrect .* (10^(params.dbscl/20)); 
    feedbackWrong = feedbackWrong ./ sqrt(mean(feedbackWrong.^2, 1));  % broadcasting
    feedbackWrong = feedbackWrong .* (10^(params.dbscl/20));
    % fill special buffers in advance with response sounds
    feedbackCorrectBuffer = PsychPortAudio('CreateBuffer', [], feedbackCorrect');
    feedbackWrongBuffer = PsychPortAudio('CreateBuffer', [], feedbackWrong');
    % user message
    disp('Scaled and buffered feedback sounds.');
end  % if feedbackFLAG


% user message
disp([char(10), 'PsychPortAudio initialized', char(10)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the stimulus, prepare for trial start

% user message
disp([char(10), 'Done with general settings, preparing the training block...']);

% set flag for aborting experiment
abortFlag = 0;

% user message
disp([char(10), 'Generating stimuli and preallocating / appending response variables...']);

% Generate the stimulus
[sigOut, trialTypes, ~, deviantIndices, stimParams] = genSNRdeviantStimulusBlock_training(params, FGtypeFLAG);

% Scale the output signal.  First normalize to RMS == 1 and then scale
% by 'params.dbscl'.
sigOut = sigOut ./ sqrt(mean(sigOut.^2, 2));  % broadcasting
sigOut = sigOut .* (10^(params.dbscl/20)); 

% Store the deviant note indices and the general stimulus params
allDeviantIndices = deviantIndices;
allStimParams = stimParams;

% Preallocate response variables
outData.rt = nan(params.trialNo, 1);  % Numeric, response time relative to stimulus onset
outData.rtDeviant = nan(params.trialNo, 1);  % Numeric, response time relative to deviant onset
outData.accuracy = nan(params.trialNo, 1);  % Binary [1 = correct, 0 = wrong]

% Collect trial / stimlus type information for block
outData.trialTypes = trialTypes;  % stimulus / trial type, numeric code

% Collect deviant onsets (in time, relative to audio onset)
outData.deviantOnset = (deviantIndices(:, 1)-1) * params.lenNote + 0.001;  % Numeric, deviant onset time relative to stimulus onset

% Send AUDITORY targetMask (i.e. indicies of deviant/target trials)
% for the i-th block to the output arrays 
if auditoryFigDeviantFLAG
    deviantMask = trialTypes == params.figDeviantLabel; 
    targetDeviantID = params.figDeviantLabel;  % use for quick calculation of correctness
    nontargetDeviantID = params.backDeviantLabel;  % use for quick calculation of correctness
else
    deviantMask = trialTypes == params.backDeviantLabel;
    targetDeviantID = params.backDeviantLabel;  % use for quick calculation of correctness
    nontargetDeviantID = params.figDeviantLabel;  % use for quick calculation of correctness
end

% Collect mask for trials with target
outData.targetPresented = deviantMask;  % Binary [1 = target, 0 = no target]

% user message
disp('Done, now buffering all audio in advance...');    

% % fill buffers with audio for whole block
buffer = [];
for tmpIdx = 1:params.trialNo
    audioData = squeeze(sigOut(tmpIdx, :, :));
    buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
end

% user message
disp('Done.');    

% user message
disp([char(10), 'Ready to start the block!']);

% Define again the task and the response key for user, just to be safe
disp([char(10), char(10)]);
disp('_______________________TRAINING BLOCK!_______________________');
disp('-------------------------------------------------------------');
disp(['THE TASK IS TO DETECT DEVIANTS IN THE    ', upper(deviantType), '']);
disp('-------------------------------------------------------------');
disp(['STIMULI SET FOR TRAINING:                ', upper(FGtype), '']);
disp('-------------------------------------------------------------');
disp(['RESPONSE KEY FOR DEVIANT DETECTION:      ', upper(params.detectKey), '']);
disp('-------------------------------------------------------------');
disp([char(10), char(10)]);    


% start block after user pressed goKey
if strcmp(params.goKey, 'Return')
    startKey = 'ENTER'; 
else
    startKey = upper(params.goKey); 
end
disp([char(10), char(10)]);
disp('-------------------------------------------------------------');
disp(['-----------    PRESS ', startKey, ' TO START THE BLOCK!    -----------']);
disp('-------------------------------------------------------------');
disp([char(10), char(10)]);

% wait for key press
KbReleaseWait;
while 1
    [keyDown, secs, keyCode] = KbCheck;
    if keyDown && find(keyCode) == goKeyCode
        blockStartTime = secs + 3;
        break;
    elseif keyDown && find(keyCode) == abortKeyCode
        abortFlag = 1;
        break;
    end          
end

% check for abort
if abortFlag
    disp([char(10), char(10), 'USER REQUESTED ABORT!', char(10), char(10)]);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Stop', pahandle, 1, 1);
    PsychPortAudio('Close');
    if triggerFlag
        fclose(serialObj);
    end    
    return;
end   

% block start trigger
if triggerFlag
    fprintf(serialObj, params.trig.format, params.trig.blockStart + trainingType);  % block type (task type) trigger = block start trigger + training type number
end  

% user message
disp('Starting...');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trial loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trialIdx = 1:params.trialNo

    % make sure the subject is not pressing the detection key already
    KbReleaseWait;   
    
    % trial start time is determined by block start for first trial,
    % otherwise it depends on feedback:
    % if there is was a feedback sound, trial start is feedback sound start + feedback length + iti 
    % otherwise it is last trial audio onset + stimulus length + iti
    if trialIdx == 1
        trialStartTime = blockStartTime;
    else
        if feedbackFLAG && feedbackOnset ~= 0
            trialStartTime = feedbackOnset + feedbackSoundLength + params.iti;  % trial start is relative to feedback sound onset
        else
            trialStartTime = oldTrialStartTime + params.stimLength + params.iti;  % trial start is relative to previous stimulus onset
        end  % if feedbackFLAG && feedbackOnset ~= 0
    end  % if trialIdx == 1  
    
    % fill audio buffer with next stimuli
    PsychPortAudio('FillBuffer', pahandle, buffer(trialIdx));

    % user message
    disp([char(10), 'Trial ', num2str(trialIdx)]);
    
    % trial type trigger
    if triggerFlag
        fprintf(serialObj, params.trig.format, trialTypes(trialIdx));  % trial type trigger, from var "trialTypes" 
    end    

    % wait till we are 50 ms from the start of the playback
    while trialStartTime - GetSecs >= 0.05
        WaitSecs(0.01);
    end

    % blocking playback start for precision
    audioOnset = PsychPortAudio('Start', pahandle, 1, trialStartTime, 1);        

    % stimulus onset trigger
    if triggerFlag
        fprintf(serialObj, params.trig.format, params.trig.soundOnset);
    end      
    
    % user message
    disp(['Audio started at ', num2str(audioOnset-trialStartTime), ' secs relative to requested start time']);
    if trialIdx ~= 1
        disp(['Audio start relative to last audio onset was ', num2str(audioOnset-oldTrialStartTime), ' secs ',...
            '(should be ', num2str(params.stimLength + params.iti), ' secs)']);
    end
    if deviantMask(trialIdx)
        disp('There is a TARGET in the trial');
        disp(['Deviant onset is at ', num2str(outData.deviantOnset(trialIdx)), ' secs from audio onset']);
    end

    % Listen to response until we are only 200 ms before next stimulus
    respFlag = 0;
    while (audioOnset+params.stimLength + params.iti) - GetSecs >= 0.2
        [keyDown, secs, keyCode] = KbCheck;
        % if subject detected deviant
        if keyDown && find(keyCode) == detectKeyCode
            % response trigger
            if triggerFlag
                fprintf(serialObj, params.trig.format, params.trig.response);
            end   
            % set response flag, store RT
            respFlag = 1;
            outData.rt(trialIdx) = secs-audioOnset;
            outData.rtDeviant(trialIdx) = secs-(audioOnset + outData.deviantOnset(trialIdx));
            break;
        % if experiment is aborted
        elseif keyDown && find(keyCode) == abortKeyCode
            abortFlag = 1;
            break;
        end
    end

    % if abort was requested, quit
    if abortFlag
        disp([char(10), char(10), 'USER REQUESTED ABORT!', char(10), char(10)]);
        Priority(0);
        RestrictKeysForKbCheck([]);
        PsychPortAudio('Stop', pahandle, 1, 1);
        PsychPortAudio('Close');
        if triggerFlag
            fclose(serialObj);
        end        
        return;
    end        

    % check accuracy 
    accFlag = 0;
    if respFlag
        if deviantMask(trialIdx)
            outData.accuracy(trialIdx) = 1;
            accFlag = 1;
        else
            outData.accuracy(trialIdx) = 0;
        end
    elseif ~respFlag
        if ~deviantMask(trialIdx)
            outData.accuracy(trialIdx) = 1;
            accFlag = 1;
        else
            outData.accuracy(trialIdx) = 0;
        end
    end 

    % if there is feedback in training and there was a response, play
    % feedback sound
    feedbackOnset = 0;
    if feedbackFLAG && respFlag
        % stop playback if it is still ongoing
        PsychPortAudio('Stop', pahandle, 0);
        % play 'correct' / wrong sound immediately
        if accFlag
            PsychPortAudio('FillBuffer', pahandle, feedbackCorrectBuffer);
            feedbackOnset = PsychPortAudio('Start', pahandle, 1);
        else
            PsychPortAudio('FillBuffer', pahandle, feedbackWrongBuffer);
            feedbackOnset = PsychPortAudio('Start', pahandle, 1);
        end  % if accFlag
    end  % if feedbackFLAG && respFlag
         
    % Hit / False alarm trigger
    if triggerFlag
        % if response was a HIT
        if respFlag && accFlag
            WaitSecs(0.1);  % workaround for buggy trigger recognition by Micromed SD LTM EEG, need break between response trigger and hit/false alarm trigger
            fprintf(serialObj, params.trig.format, params.trig.hit);
        % if response was a FALSE ALARM
        elseif respFlag && ~accFlag
            WaitSecs(0.1);  % workaround for buggy trigger recognition by Micromed SD LTM EEG, need break between response trigger and hit/false alarm trigger
            fprintf(serialObj, params.trig.format, params.trig.falseAlarm);
        end 
    end  % if triggerFlag    
    
    % audio onset becomes the last trial start
    oldTrialStartTime = audioOnset;

    % user messages
    if respFlag && accFlag
        disp([char(10), 'Subject pressed detection key. HIT!'])
        disp(['Subject pressed detection key at ', num2str(outData.rt(trialIdx)*1000), ' msec relative to stimulus onset, '...
            'and ', num2str(outData.rtDeviant(trialIdx)*1000), ' msec relative to deviant onset']);
    elseif respFlag && ~accFlag
        disp([char(10), 'Subject pressed detection key. FALSE ALARM!'])
        disp(['Subject pressed detection key at ', num2str(outData.rt(trialIdx)*1000), ' msec relative to stimulus onset']);
    end


end  % for trialIdx



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block-level feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user message
disp([char(10), char(10), 'Done with training, providing overall feedback:', char(10)]);

% Auditory task       
% Determine TARGET-DEVIANT detection rates
tmpIDmask = outData.trialTypes == targetDeviantID;
tmpCorr = sum(outData.accuracy(tmpIDmask));
tmpTotal = sum(tmpIDmask);
percentCorrectDetect = round(100*tmpCorr/tmpTotal);
fprintf(['For **target deviant** stimuli... ' , num2str(percentCorrectDetect), '%% accuracy. \n'])

% Determine NON-TARGET-DEVIANT false alarm rates
% Only if trainingType was 3 or 4, that is, there were non-target deviants!
if ismember(trainingType, 3:4)
    tmpIDmask = outData.trialTypes == nontargetDeviantID;
    tmpCorr = sum(~outData.accuracy(tmpIDmask));
    tmpTotal = sum(tmpIDmask);
    falseAlarmNonTargetDeviant = round(100*tmpCorr/tmpTotal);
    fprintf(['False alarm for non-target deviant stimuli... ' , num2str(falseAlarmNonTargetDeviant), '%%. \n'])
end

% Determine TARGET-STANDARD false alarm rates
tmpIDmask = outData.trialTypes == (targetDeviantID-20);  % ------- NOTE  -----------> This statement is very specific to the condition label values specified in 'params_intensityFG_training.m'!!!!!!!!!!!!!!!!!
tmpCorr = sum(~outData.accuracy(tmpIDmask));
tmpTotal = sum(tmpIDmask);
falseAlarmTargetStandard = round(100*tmpCorr/tmpTotal);
fprintf(['False alarm for target standard conditions... ' , num2str(falseAlarmTargetStandard), '%%. \n'])

% Determine NON-TARGET-STANDARD false alarm rates
% Only if trainingType was 3 or 4, that is, there were non-target standards!
if ismember(trainingType, 3:4)
    tmpIDmask = outData.trialTypes == (nontargetDeviantID-20);  % ------- NOTE  -----------> This statement is very specific to the condition label values specified in 'params_intensityFG_training.m'!!!!!!!!!!!!!!!!!
    tmpCorr = sum(~outData.accuracy(tmpIDmask));
    tmpTotal = sum(tmpIDmask);
    falseAlarmNonTargetStandard = round(100*tmpCorr/tmpTotal);
    fprintf(['False alarm for non-target standard conditions... ' , num2str(falseAlarmNonTargetStandard), '%%. \n'])
end

% Determine OVERALL ACCURACY rates
tmpCorr = sum(outData.accuracy);
tmpTotal = params.trialNo;
percentALLaccuracy = round(100*tmpCorr/tmpTotal);
fprintf(['Across ALL conditions... ' , num2str(percentALLaccuracy), '%% accuracy. \n'])

% Determine OVERALL FALSE ALARM rates
tmpIDmask = outData.trialTypes == (nontargetDeviantID-20) | outData.trialTypes == (targetDeviantID-20) | outData.trialTypes == nontargetDeviantID;  % ------- NOTE  -----------> This statement is very specific to the condition label values specified in 'params_intensityFG_training.m'!!!!!!!!!!!!!!!!!
tmpCorr = sum(~outData.accuracy(tmpIDmask));
tmpTotal = sum(tmpIDmask);
falseAlarmALL = round(100*tmpCorr/tmpTotal);
fprintf(['Overall false alarm rate for non-target-deviant conditions... ' , num2str(falseAlarmALL), '%%. \n'])

% Determine D-PRIME... this applies the d' correction mentioned in 
% Doelling and Poeppel (2015) PNAS E6233-E6242 described in their
% Methods and Materials section under 'Behavioral Analysis'.
% Ultimately, it is based off of Macmillan and Kaplan (1985)
% Psychol. Bulletin, vol 98, pp. 185-199 
if percentCorrectDetect == 100; pHit = 1-1/(2*sum(deviantMask)); 
elseif percentCorrectDetect == 0; pHit = 1/(2*sum(deviantMask));
else pHit = percentCorrectDetect/100; end

if falseAlarmTargetStandard == 0; pFA = 1/(2*sum(deviantMask)); 
elseif falseAlarmTargetStandard == 100; pFA = 1-1/(2*sum(deviantMask)); 
else pFA = falseAlarmTargetStandard/100; end

dprimeResult = dprime(pHit, pFA);
fprintf(['The d'' value for relevant "target" conditions... ' , num2str(dprimeResult), '. \n'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save & cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFcn_intensityFG_training(subjID, trainingType, deviantType, expName,... 
    params, outData, allDeviantIndices, allStimParams);

Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Stop', pahandle, 1, 1);
PsychPortAudio('Close');      
if triggerFlag
    fclose(serialObj);
end
                
disp([char(10), char(10)]);
fprintf('-----------------------------------------------\n')
fprintf('Experiment completed.  Thanks for your time!\n')
fprintf('-----------------------------------------------\n')


return
