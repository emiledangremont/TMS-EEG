% PREPROCESSING TMS-EEG DATA
%
% Timo van Hattem
% Updated: 13-3-2023
% Adjusted by Emile d'Angremont, 13-9-2023

cd('/scratch/anw/edangremont/TMS-EEG/code/')
%% Clean workspace
clear
close all
clc

%% Set path
%addpath('/data/anw/anw-gold/NP/projects/data_TIPICCO/TMS_EEG/tvh/eeglab2023.0/');
%addpath(genpath('/data/anw/anw-gold/NP/projects/data_TIPICCO/TMS_EEG/tvh/eeglab2023.0/FastICA_25/'));
addpath(genpath('/scratch/anw/edangremont/TMS-EEG/data/'));
addpath('/scratch/anw/edangremont/TMS-EEG/code/'); 
addpath('/scratch/anw/edangremont/TMS-EEG/code/eeglab2023.0'); 
% addpath(genpath('/scratch/anw/edangremont/TMS-EEG/code/eeglab2023.0/FastICA_25/')); 
fprintf('Paths added!\n')

%% Initialize variables
ppn = 'TC923'; %input subject number % make gui out of this
br = 'rDLPFC'; %input brain region of interest

% set input and output paths
DATAIN = ['/scratch/anw/edangremont/TMS-EEG/data/raw/', ppn, '/', br, '/'];
%DATAIN = ['/scratch/anw/tvanhattem/analysis_tvh/TMSEEG_data/convert/', ppn, '/', br, '/'];
DATAOUT = '/scratch/anw/edangremont/TMS-EEG/data/processed/';

%% Import data
eeglab;
EEG = loadcurry([DATAIN, dir(fullfile(DATAIN, '*.cdt')).name], 'KeepTriggerChannel', 'False', 'CurryLocations', 'False');
%EEG = pop_biosig(); %for loading convert file
fprintf('N_events start processing: %d\n', length(EEG.event)); % manually add to excel file?

%% Load channel locations
EEG = pop_chanedit(EEG,'lookup','/scratch/anw/edangremont/TMS-EEG/code/eeglab2023.0/plugins/dipfit/standard_BEM/elec/standard_1020.elc');

%% Remove unused electrodes
EEG = pop_select(EEG, 'nochannel', 63:68); % these are not used
EEG.allchan = EEG.chanlocs;

%% Automated removal bad electrodes step 1
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'Highpass', 'off','ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'BurstCriterion','off','WindowCriterion','off'); % should these parameters be saved somewhere?
% EEG.rejchan = find(~ismember([EEG.allchan.urchan], [EEG.chanlocs.urchan]));
rejchan_one = setdiff({EEG.allchan.labels}, {EEG.chanlocs.labels}); % I think this is superior to the previous line
fprintf('Rejected channel step 1: %s\n',rejchan_one{:}); % manually copy into excel file?

%% Automated removal bad electrodes step 2
EEG = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1), 'threshold', 4, 'norm', 'on', 'measure', 'kurt'); % same for these parameters
% EEG.rejchan = find(~ismember([EEG.allchan.urchan], [EEG.chanlocs.urchan]));
rejchan_two = setdiff(setdiff({EEG.allchan.labels}, {EEG.chanlocs.labels}),rejchan_one);
fprintf('Rejected channel step 2: %s\n',rejchan_two{:}); % manually copy into excel file?
EEG.rejchan = [rejchan_one rejchan_two];

% manual rejection of electrodes should be added here
%% Fix latency of events in D2 and D10 condition (from marker on conditioning pulse to marker on test pulse)
EEG.oldeventlatency = [EEG.event.latency];
for i = 1:size(EEG.event,2)
    if EEG.event(i).type == 3
        EEG.event(i).latency = EEG.event(i).latency + 21;
        EEG.urevent(i).latency = EEG.urevent(i).latency + 21;
    elseif EEG.event(i).type == 5
        EEG.event(i).latency = EEG.event(i).latency + 101;
        EEG.urevent(i).latency = EEG.urevent(i).latency + 101;
    end
end

%% Epoch segmentation
EEG = pop_epoch(EEG, {'1', '3', '5'}, [-1.5 1.5]);

%% Baseline correction
EEG = pop_rmbase(EEG, [-800 -110]);
%pop_eegplot(EEG,1,1,1);

%% Seperate epochs on condition
EEG_SP = pop_select(EEG, 'trial', find([EEG.event.type] == 1));
EEG_D2 = pop_select(EEG, 'trial', find([EEG.event.type] == 3));
EEG_D10 = pop_select(EEG, 'trial', find([EEG.event.type] == 5));

%% Save additional information for later checks
EEG_SP.rawepochs = EEG_SP.epoch;
EEG_SP.rawurevents = EEG_SP.urevent;

EEG_D2.rawepochs = EEG_D2.epoch;
EEG_D2.rawurevents = EEG_D2.urevent;

EEG_D10.rawepochs = EEG_D10.epoch;
EEG_D10.rawurevents = EEG_D10.urevent;

%% Manually check for true presence of TMS-pulse at given marker SP
EEG_pulsecheck_SP = epoch2continuous(EEG_SP);
EEG_pulsecheck_SP = tesa_findpulse(EEG_pulsecheck_SP, 'CZ', 'refract', 10, 'rate', 2e4, 'tmsLabel', 'SP'); % CZ, but PZ can be used if CZ was already filtered out
EEG_SP.pulseinfo = EEG_pulsecheck_SP.event;
fprintf('Is number of single pulses detected equal to %d?\n',EEG_SP.trials);
% pop_eegplot(EEG_pulsecheck_SP,1,1,1);
% EEG_SP = pop_select(EEG_SP, 'notrial', [44:51])

%% Manually check for true presence of TMS-pulse at given marker D2
EEG_pulsecheck_D2 = epoch2continuous(EEG_D2);
EEG_pulsecheck_D2 = tesa_findpulse(EEG_pulsecheck_D2, 'CZ', 'refract', 2, 'rate', 2e4, 'paired', 'yes', 'ISI', 2);
EEG_D2.pulseinfo = EEG_pulsecheck_D2.event;
fprintf('Is number of test pulses detected equal to %d?\n',EEG_D2.trials); % dit automatiseren (latencies vergelijken en epoch verwijderen waar niet overeen)
% pop_eegplot(EEG_pulsecheck_D2,1,1,1);
% EEG_D2 = pop_select(EEG_D2, 'notrial', [43:44]);

%% Manually check for true presence of TMS-pulse at given marker D10
EEG_pulsecheck_D10 = epoch2continuous(EEG_D10);
EEG_pulsecheck_D10 = tesa_findpulse(EEG_pulsecheck_D10, 'CZ', 'refract', 10, 'rate', 2e4, 'paired', 'yes', 'ISI', 10);
EEG_D10.pulseinfo = EEG_pulsecheck_D10.event;
fprintf('Is number of test pulses detected equal to %d?\n',EEG_D10.trials);
% pop_eegplot(EEG_pulsecheck_D10,1,1,1);
% EEG_D10 = pop_select(EEG_D10, 'notrial', [21]);

%% Number of raw epochs per condition
fprintf('N_rawepochs (SP/D2/D10): %d/%d/%d\n', length(EEG_SP.epoch),... % klopt dit nu?
    length(EEG_D2.epoch), length(EEG_D10.epoch));

%% Seperate epochs for conditions and save
mkdir([DATAOUT, '/', ppn, '/', br])
pop_saveset(EEG_SP, 'filename', [ppn, '_rawepochs_', br, '_SP.set'], 'filepath', [DATAOUT, ppn, '/', br]);
pop_saveset(EEG_D2, 'filename', [ppn, '_rawepochs_', br, '_D2.set'], 'filepath', [DATAOUT, ppn, '/', br]);
pop_saveset(EEG_D10, 'filename', [ppn, '_rawepochs_', br, '_D10.set'], 'filepath', [DATAOUT, ppn, '/', br]);