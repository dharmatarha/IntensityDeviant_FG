function intensityFG_intro
%% Function for introducing the four stimuli categories of the Intensity-deviant-detection FG task
%
% USAGE: intensityFG_intro
%
% Simple command line interface for generating and playing stimuli of given
% types one-by-one. 
% Pressing the keys "1" to "4" generates and plays a stimulus with:
% "1"   - FG stimulus with Figure, no deviant
% "2"   - FG stimulus without Figure (only Background), no deviant
% "3"   - FG stimulus with Figure, deviant Figure
% "4"   - FG stimulus without Figure (only Background), deviant Background
%
% Pressing "escape" aborts the task.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PsychDefaultSetup(1);
Priority(1);

% load all parameters
params = params_intensityFG_intro;


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

% user message
disp([char(10), 'PsychPortAudio initialized', char(10)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample sounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp([char(10), char(10)]);
disp('------------------------------------------------------');
disp('------------   STIMULUS INTRODUCTION   ---------------');
disp('------------------------------------------------------');
disp('PRESS "1" FOR NON-DEVIANT FIGURE');
disp('PRESS "2" FOR NON-DEVIANT BACKGROUND');
disp('PRESS "3" FOR DEVIANT FIGURE');
disp('PRESS "4" FOR DEVIANT BACKGROUND');
disp('PRESS "ESCAPE" TO ABORT');
disp([char(10), char(10)]);

abortFlag = 0;
while 1
    
    request = [];
    
    [keyDown, ~, keyCode] = KbCheck;
    if keyDown
        if find(keyCode) == KbName('1')
            request = 'figStandardLabel';
        elseif find(keyCode) == KbName('2')
            request = 'backStandardLabel';
        elseif find(keyCode) == KbName('3')
            request = 'figDeviantLabel';            
        elseif find(keyCode) == KbName('4')
            request = 'backDeviantLabel';        
        elseif find(keyCode) == KbName('escape')
            abortFlag = 1;  
        end
    else
        WaitSecs(0.05);
    end
    
    if abortFlag 
        Priority(0);
        PsychPortAudio('Stop', pahandle, 1, 1);
        PsychPortAudio('Close', pahandle);
        break;
    end
    
    if ~isempty(request)
        % Create the (random) indices for the deviant
        tmp = randi([params.deviantIndMin, params.deviantIndMax],1);
        deviantInd = tmp:(tmp+params.deviantDuration-1);

        switch request
            case 'figStandardLabel'  % Generate NON-DEVIANT FIGURE
                [sigOut, ~] = osullivanTokenIntensityDeviantFig(params.fs, params.clFig, params.clNotFig,...
                    round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                    params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                    params.maxSemitoneStep, params.clStart, [1 1], params.plotMe,...
                    params.baseSNR, deviantInd, params.NONdeviantSNR);   

            case 'backStandardLabel'  % Generate NON-DEVIANT BACKGROUND
                [sigOut, ~] = osullivanTokenIntensityDeviantBack(params.fs, params.clBack, params.clNotBack, ...
                    round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                    params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                    params.maxSemitoneStep, params.clStart, [1 1], params.plotMe,...
                    params.deviantNum, deviantInd, params.NONdeviantSNR);

            case 'figDeviantLabel'  % Generate DEVIANT FIGURE
                [sigOut, ~] = osullivanTokenIntensityDeviantFig(params.fs, params.clFig, params.clNotFig,...
                    round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                    params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                    params.maxSemitoneStep, params.clStart, [1 1], params.plotMe,...
                    params.baseSNR, deviantInd, params.deviantSNR);   

            case 'backDeviantLabel'  % Generate DEVIANT BACKGROUND
                [sigOut, ~] = osullivanTokenIntensityDeviantBack(params.fs, params.clBack, params.clNotBack, ...
                    round(params.lenNote*params.fs),  round(params.lenRamp*params.fs),  params.numNotesPerToken,...
                    params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                    params.maxSemitoneStep, params.clStart, [1 1], params.plotMe,...
                    params.deviantNum, deviantInd, params.deviantSNR); 

            otherwise
                error('DKR: Invalid stimulus type.')
        end   % switch request
        
        % Scale the output signal.  First normalize to RMS == 1 and then scale
        % by 'params.dbscl'.
        sigOut = sigOut ./ sqrt(mean(sigOut.^2, 1));  % broadcasting
        sigOut = sigOut .* (10^(params.dbscl/20)); 

        % buffer into PsychPortAudio and play
        buffer = PsychPortAudio('CreateBuffer', [], sigOut');
        PsychPortAudio('FillBuffer', pahandle, buffer);
        PsychPortAudio('Start', pahandle, 1);
        KbReleaseWait;
        
    end  % if ~isempty(request)
    
    
end