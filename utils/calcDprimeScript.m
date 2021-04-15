% first load a participant's .mat file with the 'paramsOut.outData' structure

expType = 'figure';  % Either 'figure' or 'background' depending on if the listener was to detect
col =  size(paramsOut.outData.PresentedCondID,2);  % Specify which block to analyze... this is setup to select the last block

switch lower(expType)
    case 'figure'
        % For 'FIGURE Deviant detection'
        mskDev = paramsOut.outData.PresentedCondID(:,col) == 30;
        mskStd = paramsOut.outData.PresentedCondID(:,col) == 10;
%         mskStd = paramsOut.outData.PresentedCondID(:,col) == 20 | paramsOut.outData.PresentedCondID(:,col) == 10 | paramsOut.outData.PresentedCondID(:,col) == 40;
    case 'background'
        % % For 'BACKGROUND Deviant detection'
        mskDev = paramsOut.outData.PresentedCondID(:,col) == 40;
        mskStd = paramsOut.outData.PresentedCondID(:,col) == 20;
%         mskStd = paramsOut.outData.PresentedCondID(:,col) == 20 | paramsOut.outData.PresentedCondID(:,col) == 10 | paramsOut.outData.PresentedCondID(:,col) == 30;
    otherwise
        error('DKR: Invalid condition type.')
end


hit = sum(paramsOut.outData.Responses(mskDev,col));
fa = sum(~paramsOut.outData.Responses(mskStd,col));
totalStd = sum(mskStd);
totalDev = sum(mskDev);

percentCorrect = hit/totalDev
falseAlarm = fa/totalStd

% Determine D-PRIME... this applies the d' correction mentioned in 
% Doelling and Poeppel (2015) PNAS E6233-E6242 described in their
% Methods and Materials section under 'Behavioral Analysis'.
% Ultimately, it is based off of Macmillan and Kaplan (1985)
% Psychol. Bulletin, vol 98, pp. 185-199 
if percentCorrect == 1; pHit = 1-1/(2*totalDev); 
elseif percentCorrect == 0; pHit = 1/(2*totalDev);
else pHit = percentCorrect; end

if falseAlarm == 0; pFA = 1/(2*totalDev); 
elseif falseAlarm == 1; pFA = 1-1/(2*totalDev); 
else pFA = falseAlarm; end


dprimeResult = dprime(pHit, pFA)
