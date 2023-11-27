% PREPROCESSING SINGLE-PULSE TMS-EEG DATA
%
% Timo van Hattem
% Updated: 13-3-2023
% Adjusted by Emile d'Angremont, 13-9-2023
cd('/scratch/anw/edangremont/TMS-EEG/')

%% Clean workspace
clear
close all
clc 
%% Set path
%addpath('/data/anw/anw-gold/NP/projects/data_TIPICCO/TMS_EEG/tvh/eeglab2023.0/');
%addpath(genpath('/data/anw/anw-gold/NP/projects/data_TIPICCO/TMS_EEG/tvh/eeglab2023.0/FastICA_25/'));
addpath(genpath('/scratch/anw/edangremont/TMS-EEG/data/'));
addpath('/scratch/anw/edangremont/TMS-EEG/code/'); 
addpath('/scratch/anw/edangremont/TMS-EEG/code/eeglab2023.0')
addpath(genpath('/scratch/anw/edangremont/TMS-EEG/code/eeglab2023.0/FastICA_25/')); 
fprintf('Paths added!\n')

%% Initialize variables
ppn = 'TC923'; %input subject number
br = 'rDLPFC'; %input brain region of interest
con = 'SP'; %input condition

%set input and output paths
DATAIN = ['/scratch/anw/edangremont/TMS-EEG/data/processed/', ppn, '/', br, '/'];
DATAOUT = '/scratch/anw/edangremont/TMS-EEG/data/processed/';

%% Import data
eeglab;
EEG = pop_loadset('filename', [ppn, '_rawepochs_', br, '_', con, '.set'], 'filepath', DATAIN);

%% Remove TMS artefact 
if strcmpi(con, 'SP')
    EEG = pop_tesa_removedata(EEG, [-5 15]);
elseif strcmpi(con, 'D2')
    EEG = pop_tesa_removedata(EEG, [-7 15]);
elseif strcmpi(con, 'D10')
    EEG = pop_tesa_removedata(EEG, [-15 15]);
end

%% Interpolate 
EEG = pop_tesa_interpdata(EEG, 'cubic', [1 1]);

%% Downsample
EEG = pop_resample(EEG, 1000);

%% Automated removal bad epochs
EEG = pop_jointprob(EEG,1,1:size(EEG,1),3,3,0,1,0,[],0); % This looks function rejects epochs in which there is strong acitivty in all electrodes? should we save these parameters?
% EEG.rejepoch = find(~ismember([EEG.rawepochs.eventurevent], [EEG.epoch.eventurevent]));
EEG.rejepoch = setdiff([EEG.rawepochs.eventurevent],[EEG.epoch.eventurevent]); % NB this does something different. Do we want the event or the eventurevent?
fprintf('Number of rejected epochs: %d\n', length(EEG.rejepoch));

%% Manual check data
% trial rejection
pop_eegplot(EEG,1,1,1)
%% Rejected epochs?
EEG.rejepoch = find(~ismember([EEG.rawepochs.eventurevent], [EEG.epoch.eventurevent]));
fprintf('Number of rejected epochs including manual: %d\n', length(EEG.rejepoch));
%% bad channel rejection
EEG = pop_select(EEG, 'nochannel', {'TP9'}); %manual input, create GUI?
EEG.rejchan = find(~ismember([EEG.allchan.urchan], [EEG.chanlocs.urchan]));
EEG.rejchan = setdiff({EEG.allchan.labels},{EEG.chanlocs.labels});

%% Remove data
if strcmpi(con, 'SP')
    EEG = pop_tesa_removedata(EEG, [-5 15]);
elseif strcmpi(con, 'D2')
    EEG = pop_tesa_removedata(EEG, [-7 15]);
elseif strcmpi(con, 'D10')
    EEG = pop_tesa_removedata(EEG, [-15 15]);
end

%% ICA round 1
EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
EEG = pop_tesa_compselect(EEG,'compCheck','on','remove','on','saveWeights','off', 'figSize','medium','plotTimeX',[-200 500],'plotFreqX',[1 100],'freqScale', 'log', 'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[15 30],'tmsMuscleFeedback','off','blink','off','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','off','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',-0.31,'muscleFreqIn',[7 70],'muscleFreqEx',[48 52],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',4,'elecNoiseFeedback','off' );

%% Interpolate 
EEG = pop_tesa_interpdata(EEG, 'cubic', [10 10]);

%% Filtering
EEG = pop_tesa_filtbutter(EEG, 1, 80, 4, 'bandpass');
EEG = pop_tesa_filtbutter(EEG, 48, 52, 4, 'bandstop' );
%EEG_erp = pop_tesa_filtbutter(EEG, 1, 45, 4, 'bandpass');

%% SOUND
%EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 5);

%% Manual trial rejection 
pop_eegplot(EEG,1,1,1)
%% Reject epochs?
% EEG = pop_select(EEG, 'notrial', [27]); %manual input
% EEG.rejepoch = find(~ismember([EEG.rawepochs.eventurevent], [EEG.epoch.eventurevent]));

%% Remove data around TMS-puls and add constant amplitude
if strcmpi(con, 'SP')
    EEG = pop_tesa_removedata(EEG, [-5 15]);
elseif strcmpi(con, 'D2')
    EEG = pop_tesa_removedata(EEG, [-7 15]);
elseif strcmpi(con, 'D10')
    EEG = pop_tesa_removedata(EEG, [-15 15]);
end

%% ICA round 2
EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
EEG = pop_tesa_compselect(EEG,'compCheck','on','remove','on','saveWeights','off', 'figSize','medium','plotTimeX',[-200 500],'plotFreqX',[1 100],'freqScale', 'log', 'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[15 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'AF7','AF8'},'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','on','muscleThresh',-0.31,'muscleFreqIn',[7 70],'muscleFreqEx',[48 52],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );

%% Interpolate data around TMS-pulse
EEG = pop_tesa_interpdata(EEG, 'cubic', [10 10]); 

%% Number of epochs and channels after processing
fprintf('N_epochs (%s): %d\n', con, length(EEG.epoch));
fprintf('N_channels (%s): %d\n', con, length(EEG.chanlocs));

%% Interpolate removed channels
EEG = pop_interp(EEG, EEG.allchan, 'spherical');

%% Rereference data to average of all electrodes
EEG = pop_reref(EEG, []);
%pop_eegplot(EEG,1,1,1)

%% Save dataset
mkdir([DATAOUT, '/', ppn, '/', br])
pop_saveset(EEG, 'filename', [ppn, '_preprocessed_', br, '_', con, '.set'], 'filepath', [DATAOUT, '/', ppn, '/' br]);

