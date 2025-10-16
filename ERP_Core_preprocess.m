% ERP Core - Whole Brain analysis
% ----------------------------------------
% This script uses EEGLAB to preprocess the data.
% Data are saved and reused in other script for analyses. 
% Cyril Pernet 16 Octobre 2015

clear variables
InputDataset   = '/indirect/staff/cyrilpernet/multiverse_analyses/ERP_CORE_BIDS_Raw_Files';
OutputLocation = '/indirect/staff/cyrilpernet/multiverse_analyses/GAEtesting';
ERP_Core_check_install

% task we process
% Note ERN was removed as there is unequal design matrix depending on error
% responses - while this can be solved merging all conditions, for testing
% purposes here, we simply did not process that task
% N2Pc could be analyzed but since the interpretation relies on
% hemispheric differences, it does not reflect dirreclty what the GAE could
% be representing (given we that we use here channel neighbouring for the
% graph - adding hemispheric connection can/could change results, TBC).
TaskLabel       = {'MMN','N170','N400','P3'};

% task parameters
epoch_window    = repmat([-0.2 0.8],length(TaskLabel),1);
baseline_window = repmat([-200 0],length(TaskLabel),1);

% subject list
all_sub = dir(fullfile(InputDataset,'sub-*'));
SubjectLabel = arrayfun(@(x) x.name, all_sub, 'UniformOutput', false)';
sublist = find(ismember({all_sub.name}', SubjectLabel))'; % labels to num

% start eeglab and check plug-ins
rng('default');
current_folder = pwd;
if ~exist('OutputLocation','dir')
    mkdir(OutputLocation)
end

% analysis parameters
high_pass    = 0.5;
tmpfilter    = 'off';
ICAname      = 'picard';

% edit participants.tsv checking the same subjects are present
participants = readtable(fullfile(InputDataset,'participants.tsv'), 'FileType', 'text', ...
    'Delimiter', '\t', 'TreatAsEmpty', {'N/A','n/a'}); N = size(participants,1);
for p=length(participants.participant_id):-1:1
    name_match(:,p) = arrayfun(@(x) strcmpi(x.name,participants.participant_id{p}),all_sub);
end

if ~isempty(find(sum(name_match,1)==0)) %#ok<EFIND>
    participants(find(sum(name_match,1)==0),:) = []; %#ok<FNDSB>
    warning('mismatch between files and participants.tsv -%g subject(s)',N-size(participants,1))
    writetable(participants, fullfile(InputDataset,'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');
end

% edit events.tsv files
% should we correct epoching +26ms for stimuli from events.tsv files? as opposed to eeg channels

% edit events.tsv files for meaningful epoching for N170
if any(contains(TaskLabel,'N170'))
    for sub = 1:size(all_sub,1)
        root   = fullfile(all_sub(sub).folder,[all_sub(sub).name filesep 'ses-N170' filesep 'eeg']);
        file   = [all_sub(sub).name,'_ses-N170_task-N170_events.tsv'];
        if exist(fullfile(root,file),'file')
            events = readtable(fullfile(root,file), 'FileType', 'text', ...
                'Delimiter', '\t', 'TreatAsEmpty', {'N/A','n/a'});
            for s = size(events,1):-1:1
                if events.value(s) <= 40
                    event{s} = 'faces';
                elseif (events.value(s) >= 41) && (events.value(s) < 101)
                    event{s} = 'cars';
                elseif (events.value(s) >= 101) && (events.value(s) < 141)
                    event{s} = 'scrambled_faces';
                else
                    event{s} = 'scrambled_cars';
                end
            end
            t = table(events.onset,events.duration,events.sample,events.trial_type,event',events.value,...
                'VariableNames',{'onset', 'duration', 'sample', 'trial_type', 'event', 'value'});
            writetable(t, fullfile(root,file), 'FileType', 'text', 'Delimiter', '\t');
            clear event events t
        end
    end
end

% loop by TaskLabel
for t = 1:length(TaskLabel)

    %% IMPORT
    outdir = fullfile(OutputLocation);
    if ~exist(outdir,'dir')
        mkdir(outdir)
    end

    if strcmpi(TaskLabel{t},'N170')
        [STUDY, ALLEEG] = pop_importbids(InputDataset, 'bidsevent','on','bidschanloc','on', ...
            'bidstask',TaskLabel{t},'eventtype', 'event', 'outputdir' ,outdir, 'studyName',TaskLabel{t}, 'subjects', sublist);
    else
        [STUDY, ALLEEG] = pop_importbids(InputDataset, 'bidsevent','on','bidschanloc','on', ...
            'bidstask',TaskLabel{t},'eventtype', 'value', 'outputdir' ,outdir, 'studyName',TaskLabel{t}, 'subjects', sublist);
    end

    if t == 1 % also export metadata
        addpath([fileparts(which('pop_importbids.m')) filesep 'JSONio']);
        json = jsonread([InputDataset filesep 'dataset_description.json']);
        json.DatasetType = 'Derivative';
        json.Authors = 'Cyril Pernet';
        json.SourceDatasets = "https://osf.io/9f5w7/files/osfstorage";
        jsonwrite(fullfile(outdir,'dataset_description.json'),json,'prettyprint','on');
        % ignore extra files
        lines = {'*.study', '*.mat'};
        fid = fopen([outdir filesep '.bidsignore'], 'w');
        if fid == -1
            error('Cannot open .bidsignore for writing.');
        else
            for i = 1:length(lines)
                fprintf(fid, '%s\n', lines{i});
            end
            fclose(fid);
        end
    end

    if length(ALLEEG) == 1 %#ok<ISCL>
        ALLEEG = eeg_checkset(ALLEEG, 'loaddata');
    end
    ALLEEG = pop_select( ALLEEG, 'nochannel',{'HEOG_left','HEOG_right','VEOG_lower'});
    STUDY = pop_statparams(STUDY, 'default');

    [STUDY,~,AvgChanlocs] = std_prepare_neighbors(STUDY, ALLEEG, 'force', 'on');
    % remove connections 8-9/3 ie P7-P9/F7, 26-27/19 ie P8-P10/F8 and 7-25/22 ie P3-P4/Cz
    pairs(1,:) = [3 8];   pairs(2,:) = [3 9];
    pairs(3,:) = [19 26]; pairs(4,:) = [19 27];
    pairs(5,:) = [7 22];  pairs(6,:) = [25 22];
    for p=1:6
        AvgChanlocs.channeighbstructmat(pairs(p,1),pairs(p,2)) = 0;
        AvgChanlocs.channeighbstructmat(pairs(p,2),pairs(p,1)) = 0;
    end
    save(fullfile(outdir, [TaskLabel{t} '-AvgChanlocs.mat']),'AvgChanlocs')

    %% Pre-processing
    % for each subject, downsample, clean 50Hz, remove bad channels,
    % interpolate, re-reference to the average, run ICA to remove
    % eye and muscle artefacts, delete bad segments

    EEG = ALLEEG;
    for s=1:size(ALLEEG,2)
        try
            % downsample
            EEGTMP = eeg_checkset(EEG(s), 'loaddata');
            if EEGTMP.srate ~= 250
                EEGTMP = pop_resample(EEGTMP, 250);
            end
            % line freq removal
            EEGTMP = pop_zapline_plus(EEGTMP,'noisefreqs','line',...
                'coarseFreqDetectPowerDiff',4,'chunkLength',30,...
                'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);
            % remove bad channels
            EEGTMP = pop_clean_rawdata(EEGTMP,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
                'LineNoiseCriterion',4,'Highpass',[high_pass-0.25 high_pass+0.25] ,...
                'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
                'Distance','Euclidian','WindowCriterionTolerances','off' );
            % interpolate missing channels and reference
            [~,idx] = setdiff({AvgChanlocs.expected_chanlocs.labels},{EEGTMP.chanlocs.labels});
            stats.interpolated_channels= idx;
            if ~isempty(idx)
                EEGTMP = pop_interp(EEGTMP, AvgChanlocs.expected_chanlocs(idx), 'sphericalKang');
            end

            % ICA cleaning
            if ~strcmpi(tmpfilter,'off')
                tmpeeg = pop_eegfiltnew(EEGTMP,tmpfilter,0); %#ok<UNRCH>
                if strcmpi(ICAname,'picard')
                    tmpeeg = pop_runica(tmpeeg, 'icatype',ICAname,'maxiter',500,'mode','standard','concatcond','on', 'options',{'pca',EEGTMP.nbchan-1});
                else
                    tmpeeg = pop_runica(tmpeeg, 'icatype',ICAname,'concatcond','on', 'options',{'pca',tmpeeg.nbchan-1});
                end
                % project solution to the 0.5Hz filtered data
                EEGTMP.icasphere   = [];
                EEGTMP.icasphere   = tmpeeg.icasphere;
                EEGTMP.icaweights  = [];
                EEGTMP.icaweights  = tmpeeg.icaweights;
                EEGTMP.icachansind = [];
                EEGTMP.icachansind = tmpeeg.icachansind;
                EEGTMP.icaact      = [];
                EEGTMP.icawinv     = [];
                EEGTMP             = eeg_checkset(EEGTMP); % re-compute EEG.icawinv
                clear tmpeeg
            else
                if strcmpi(ICAname,'picard')
                    EEGTMP = pop_runica(EEGTMP, 'icatype',ICAname,'maxiter',500,'mode','standard','concatcond','on', 'options',{'pca',EEGTMP.nbchan-1});
                else
                    EEGTMP = pop_runica(EEGTMP, 'icatype',ICAname,'concatcond','on', 'options',{'pca',EEGTMP.nbchan-1}); %#ok<UNRCH>
                end
            end
            % clean data with ICs using IClabels
            EEGTMP = pop_iclabel(EEGTMP, 'default');
            EEGTMP = pop_icflag(EEGTMP,[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
            if isempty(EEGTMP.reject)
                stats.removed_ica_components = 0;
            else
                stats.removed_ica_components = sum(EEGTMP.reject.gcompreject);
            end
            EEGTMP = pop_subcomp(EEGTMP,[],0);

            % clean data using ASR - just the bad segment
            EEGTMP = pop_clean_rawdata(EEGTMP,'FlatlineCriterion','off','ChannelCriterion','off',...
                'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
                'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
                'WindowCriterionTolerances',[-Inf 7] );
            stats.percentage_removed_data = 1-(EEGTMP.pnts/EEG(s).pnts);

            % re-reference
            EEGTMP = pop_reref(EEGTMP,[],'interpchan','off');
            EEGTMP = pop_saveset(EEGTMP,'savemode','resave');
            EEG    = eeg_store(EEG, EEGTMP, s); % does EEG(s) = EEGTMP but with extra checks
            jsonwrite(fullfile(STUDY.filepath,[filesep EEGTMP.subject 'task_' TaskLabel{t} '_stats.json']), stats, 'prettyprint','on');
        catch pipe_error
            error_report{s} = pipe_error.message; %#ok<SAGROW>
        end
    end

    % Save study
    if exist('error_report','var')
        mask = cellfun(@(x) ~isempty(x), error_report); % which subject/session
        if all(mask)
            save(fullfile(OutputLocation,'error_report_preprocessing'),error_report);
            error('there has been a preprocessing issue with all included datasets, cannot proceed');
        else
            STUDY = std_rmdat(STUDY, EEG, 'datinds', find(mask));
            EEG(mask) = [];
        end
    end
    ALLEEG = EEG;

    % Extract data epochs (windowing as per ERP core github)
    if strcmpi(TaskLabel{t},'ERN')
        eventTypes = {'111','112','121','122','211','212','221','222'};
        EEG = pop_epoch(ALLEEG,eventTypes,epoch_window(t,:) ,'epochinfo','yes');
    elseif strcmpi(TaskLabel{t},'MMN')
        eventTypes = {'80','70'};
        EEG = pop_epoch(ALLEEG,eventTypes,epoch_window(t,:) ,'epochinfo','yes');
    elseif strcmpi(TaskLabel{t},'N170')
        eventTypes = {'faces','cars','scrambled_faces','scrambled_cars'};
        EEG = pop_epoch(ALLEEG,eventTypes, epoch_window(t,:) ,'epochinfo','yes');
    elseif strcmpi(TaskLabel{t},'N2pc')
        eventTypes = {'111','112','121','122','211','212','221','222'};
        EEG = pop_epoch(ALLEEG,eventTypes,epoch_window(t,:) ,'epochinfo','yes');
    elseif strcmpi(TaskLabel{t},'N400')
        eventTypes = {'111','112','121','122','211','212','221','222'};
        EEG = pop_epoch(ALLEEG,eventTypes, epoch_window(t,:) ,'epochinfo','yes');
    elseif strcmpi(TaskLabel{t},'P3')
        eventTypes = {'11','12','13','14','15','21','22','23','24','25',...
            '31','32','33','34','35','41','42','43','44','45','51','52','53','54','55'};
        EEG = pop_epoch(ALLEEG,eventTypes, epoch_window(t,:) ,'epochinfo','yes');
    end
    EEG    = eeg_checkset(EEG);
    EEG    = pop_saveset(EEG, 'savemode', 'resave');
    if any(strcmpi(TaskLabel{t},{'ERN','N170'}))
        [STUDY, EEG] = std_editset(STUDY, EEG, 'commands',{{'remove',4}},'updatedat','on','rmclust','on');
    elseif any(strcmpi(TaskLabel{t},{'P3'}))
        [STUDY, EEG] = std_editset(STUDY, EEG, 'commands',{{'remove',4},{'remove',35}},'updatedat','on','rmclust','on');
    end

    % Create study design
    STUDY  = std_checkset(STUDY, EEG);
    STUDY  = std_makedesign(STUDY, EEG, 1, 'name',TaskLabel{t}, ...
        'delfiles','off','defaultdesign','off','variable1','type','values1',{});

    % Precompute ERP measures
    [STUDY, EEG] = std_precomp(STUDY, EEG, {}, 'savetrials','on','interp','on','recompute','on',...
        'erp','on','erpparams', {'rmbase' baseline_window(t,:)}, ...
        'spec','off','ersp','off','itc','off');

    % output preprocessed files
    for s=size(EEG,2):-1:1
        old = fullfile(EEG(s).filepath,EEG(s).filename(1:end-4));
        EEG(s).setname = 'preprocessed';
        EEG(s).filename = [EEG(s).filename(1:end-7) 'desc-preprocessed_eeg.set'];
        EEG(s) = pop_saveset(EEG(s), 'filename', [EEG(s).filename(1:end-7) 'desc-preprocessed_eeg.set'], 'filepath', EEG(s).filepath);
        STUDY.datasetinfo(s).filename = EEG(s).filename;
        delete([old '.set']);
        delete([old '.fdt']);
    end
    STUDY  = std_checkset(STUDY, EEG);
    pop_savestudy(STUDY,EEG,'savemode','resave');
end