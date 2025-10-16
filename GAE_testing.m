% ERP Core - Whole Brain analysis
% ----------------------------------------
% The workflow below uses EEGLAB and LIMO MEEG to 
% perform that statistical analyses - the 1st level 
% analysis is done with and without trial weighting (see
% https://doi.org/10.52294/ApertureNeuro.2022.2.SEOO9435)
% and the group level analysis with and without
% subject weighting. Subject weighting is done using
% the Graph Based Auto Encoder. Results are all saved
% and reused in another script for analyses. One for
% analysising statistical outputs and the other for xAI
% % of the GAE model.

clear variables
eeglab('nogui')
outdir = '/indirect/staff/cyrilpernet/multiverse_analyses/GAEtesting';

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
analysis_window = repmat([-200 600],length(TaskLabel),1);

% start eeglab and check plug-ins
rng('default');
current_folder = pwd;

% analysis parameters
estimation   = {'OLS','WLS'};
nboot        = 1000;
tfce         = 0;

% -------------------
%% Statistics
% -----------------
for t = 1 % 1:length(TaskLabel)
    clear STUDY EEG
    [STUDY,EEG] = pop_loadstudy(fullfile(outdir,[TaskLabel{t} '.study']));
    AvgChanlocs = load(fullfile(outdir,[char(TaskLabel{t}) '-AvgChanlocs.mat']));
    AvgChanlocs = AvgChanlocs.AvgChanlocs;
    for est=1:2
        [STUDY, files] = std_limo(STUDY, EEG, 'method',estimation{est},...
            'measure','daterp', 'chanloc',AvgChanlocs,...
            'timelim',analysis_window(t,:),'erase','on',...
            'splitreg','off','interaction','off');

        if isempty(STUDY.filepath) % this seems to happen no unknown reason
            STUDY.filepath = outdir;
        end
        STUDY  = std_checkset(STUDY, EEG);
        pop_savestudy(STUDY,EEG,'savemode','resave')

        % add contrasts - which is study specific and run 2nd level
        if strcmpi(TaskLabel{t},'MMN')

            resultdir = fullfile([outdir filesep 'derivatives'],...
                ['LIMO_MMN' filesep estimation{est}]);
            mkdir(resultdir);
            cd(resultdir);
            limo_random_select('paired t-test',AvgChanlocs,...
                'LIMOfiles',fullfile(files.LIMO,"Beta_files_MMN_MMN_GLM_Channels_Time_OLS.txt"), ...
                'parameter',[1 2], 'analysis_type',...
                'Full space analysis', 'type','Channels','nboot',0,'tfce',tfce);
            limo_get_effect_size('Paired_Samples_Ttest_parameter_1_2.mat')
   
            resultdir = fullfile([outdir filesep 'derivatives'],...
                ['LIMO_MMN' filesep 'weighted_' estimation{est}]);
            mkdir(resultdir);
            cd(resultdir);
            limo_random_select('paired t-test',AvgChanlocs,...
                'LIMOfiles',fullfile(files.LIMO,"Beta_files_MMN_MMN_GLM_Channels_Time_OLS.txt"), ...
                'parameter',[1 2], 'analysis_type','Full space analysis', ...
                'method','weighted','type','Channels','nboot',0,'tfce',tfce,...
                'saveGAE','yes');
            limo_get_effect_size('Paired_Samples_Ttest_parameter_1_2.mat')

        elseif strcmpi(TaskLabel{t},'N170')

            resultdir = fullfile([outdir filesep 'derivatives'],...
                ['LIMO_N170' filesep estimation{est}]);
            mkdir(resultdir);
            cd(resultdir);

            % there are two analyses
            % faces vs cars
            % faces-scrambled vs cars-scrambled
            mkdir('Cars_vs_Faces'); cd('Cars_vs_Faces');
            limo_random_select('paired t-test',AvgChanlocs,...
                'LIMOfiles',files.Beta, 'parameter',[2 1], 'analysis_type',...
                'Full scalp analysis', 'type','Channels','nboot',nboot,'tfce',tfce);
            limo_get_effect_size('Paired_Samples_Ttest_parameter_2_1.mat')
            % ERPs (use limo_add_plots to visualize)
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                1, AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Cars')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                2, AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Faces')
            Diff = limo_plot_difference('ERPs_Faces_single_subjects_Weighted mean.mat',...
                'ERPs_Cars_single_subjects_Weighted mean.mat',...
                'type','paired','fig',0,'name','ERP_Difference');
            save('ERP_difference','Diff')

            cd(resultdir); mkdir('Cars_vs_Faces_controlled'); cd('Cars_vs_Faces_controlled'); clear data
            limo_random_select('Repeated Measures ANOVA',AvgChanlocs,...
                'LIMOfiles',{fullfile(files.LIMO,['Beta_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt'])}, ...
                'analysis type','Full scalp analysis', 'parameters',{[1 2],[3 4]}, ...
                'factor names',{'Scrambling','Category'},'type','Channels',...
                'nboot',nboot,'tfce',tfce,'skip design check','yes');
            limo_get_effect_size('Rep_ANOVA_Main_effect_1_Scrambling.mat')
            limo_get_effect_size('Rep_ANOVA_Main_effect_2_Category.mat')
            limo_get_effect_size('Rep_ANOVA_Interaction_Factors_12.mat')
            % Param avg (use limo_add_plots to visualize)
            % we can also do double diff ERP if needed (do as above twice)
            CarDiff  = [1 0 -1 0];
            FaceDiff = [0 1 0 -1];
            [~,~,LFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']));
            [~,R,BFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['Beta_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']));
            clear con1_files con2_files
            for s=length(LFiles):-1:1
                name = limo_get_subname(R{s});
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, CarDiff);
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, FaceDiff);
                con1_files{s,:} = fullfile(R{s},[name '_desc-con_1.mat']);
                con2_files{s,:} = fullfile(R{s},[name '_desc-con_2.mat']);
            end
            writecell(con1_files,fullfile(files.LIMO,'con1_files.txt'))
            writecell(con2_files,fullfile(files.LIMO,'con2_files.txt'))

            limo_central_tendency_and_ci(fullfile(files.LIMO,'con1_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_Cars')
            limo_central_tendency_and_ci(fullfile(files.LIMO,'con2_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_Faces')
            Diff = limo_plot_difference('Con_Faces_single_subjects_mean.mat',...
                'Con_Cars_single_subjects_mean.mat',...
                'type','paired','fig',0,'name','Con_diff');
            save('Parameter_difference','Diff')

            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                1 , AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Cars')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                3, AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Cars_control')
            set1      = load('ERPs_Cars_single_subjects_Weighted mean.mat');
            set2      = load('ERPs_Cars_control_single_subjects_Weighted mean.mat');
            Data.data = set1.Data.data - set2.Data.data;
            Data.limo = set1.Data.limo;
            save('ERPs_Cars_diff_single_subjects_Weighted mean.mat','Data')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                2, AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Faces')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N170_N170_GLM_Channels_Time_' estimation{est} '.txt']), ...
                4, AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_Faces_control')
            set1      = load('ERPs_Faces_single_subjects_Weighted mean.mat');
            set2      = load('ERPs_Faces_control_single_subjects_Weighted mean.mat');
            Data.data = set1.Data.data - set2.Data.data;
            Data.limo = set1.Data.limo;
            save('ERPs_Faces_diff_single_subjects_Weighted mean.mat','Data')
            Diff = limo_plot_difference('ERPs_Faces_diff_single_subjects_Weighted mean.mat',...
                'ERPs_Cars_diff_single_subjects_Weighted mean.mat',...
                'type','paired','fig',0,'name','ERP_Faces_Cars_Difference');
            save('ERP_difference','Diff')

        elseif strcmpi(TaskLabel{t},'N400')

            resultdir = fullfile([outdir filesep 'derivatives'],...
                ['LIMO_N400' filesep estimation{est}]);
            mkdir(resultdir);
            cd(resultdir);
            mkdir(resultdir);
            cd(resultdir);

            % 111 prime word, related word pair, list 1
            % 112 prime word, related word pair, list 2
            % 121 prime word, unrelated word pair, list 1
            % 122 prime word, unrelated word pair, list 2
            % 211 target word, related word pair, list 1
            % 212 target word, related word pair, list 2
            % 221 target word, unrelated word pair, list 1
            % 222 target word, unrelated word pair, list 2
            related   = [0 0 0 0 1 1 0 0];
            unrelated = [0 0 0 0 0 0 1 1];
            [~,~,LFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['LIMO_files_N400_N400_GLM_Channels_Time_' estimation{est} '.txt']));
            [~,R,BFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['Beta_files_N400_N400_GLM_Channels_Time_' estimation{est} '.txt']));

            clear con1_files con2_files
            for s=length(LFiles):-1:1
                name    = limo_get_subname(R{s}); % name
                s_value = find(arrayfun(@(x) strcmpi(x.subject,name),STUDY.limo.subjects)); % ensure match
                cond    = unique(STUDY.limo.subjects(s_value).cat_file);
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, related(cond));
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, unrelated(cond));
                con1_files{s,:} = fullfile(R{s},[name '_desc-con_1.mat']);
                con2_files{s,:} = fullfile(R{s},[name '_desc-con_2.mat']);
            end
            writecell(con1_files,fullfile(files.LIMO,'con1_files.txt'))
            writecell(con2_files,fullfile(files.LIMO,'con2_files.txt'))

            clear data
            for N=size(con1_files,1):-1:1
                data{1,N} = con1_files{N};
                data{2,N} = con2_files{N};
            end
            limo_random_select('paired t-test',AvgChanlocs,...
                'LIMOfiles',data, 'analysis_type',...
                'Full scalp analysis', 'type','Channels','nboot',nboot,'tfce',tfce);
            limo_get_effect_size('Paired_Samples_Ttest_parameter_1_2.mat')
            % Param avg (use limo_add_plots to visualize)
            limo_central_tendency_and_ci(fullfile(files.LIMO,'con1_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_related')
            limo_central_tendency_and_ci(fullfile(files.LIMO,'con2_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_unrelated')
            Diff = limo_plot_difference('Con_unrelated_single_subjects_mean.mat',...
                'Con_related_single_subjects_mean.mat',...
                'type','paired','fig',0,'name','Con_diff');
            save('Parameter_difference','Diff')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N400_N400_GLM_Channels_Time_' estimation{est} '.txt']),...
                'con_1', AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_related')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_N400_N400_GLM_Channels_Time_' estimation{est} '.txt']),...
                'con_2', AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_unrelated')
            Diff = limo_plot_difference('ERPs_unrelated_single_subjects_Weighted mean.mat',...
                'ERPs_related_single_subjects_Weighted mean.mat',...
                'type','paired','fig',0,'name','ERP_diff');
            save('ERP_difference','Diff')

        elseif strcmpi(TaskLabel{t},'P3')

            resultdir = fullfile([outdir filesep 'derivatives'],...
                ['LIMO_P3' filesep estimation{est}]);
            mkdir(resultdir);
            cd(resultdir);
            mkdir(resultdir);
            cd(resultdir);

            % 11: Stimulus - block target A, trial stimulus A,
            % 22: Stimulus - block target B, trial stimulus B,
            % 33: Stimulus - block target C, trial stimulus C,
            % 44: Stimulus - block target D, trial stimulus D,
            % 55: Stimulus - block target E, trial stimulus E,
            % 12, 13, 14, 15
            % 21, 23, 24, 25
            % 31, 32, 34, 35
            % 41, 42, 43, 45
            % 51, 52, 53, 54
            distractor = [0 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 0];
            target     = distractor==0;
            [~,~,LFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['LIMO_files_P3_P3_GLM_Channels_Time_' estimation{est} '.txt']));
            [~,R,BFiles] = limo_get_files([],[],[],...
                fullfile(files.LIMO,['Beta_files_P3_P3_GLM_Channels_Time_' estimation{est} '.txt']));

            clear con1_files con2_files
            for s=length(LFiles):-1:1
                name    = limo_get_subname(R{s}); % name
                s_value = find(arrayfun(@(x) strcmpi(x.subject,name),STUDY.limo.subjects)); % ensure match
                cond    = unique(STUDY.limo.subjects(s_value).cat_file);
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, distractor(cond));
                limo_contrast(fullfile(R{s},[name '_desc-Yr.mat']), BFiles{s}, LFiles{s}, 'T', 1, target(cond));
                con1_files{s,:} = fullfile(R{s},[name '_desc-con_1.mat']);
                con2_files{s,:} = fullfile(R{s},[name '_desc-con_2.mat']);
            end
            writecell(con1_files,fullfile(files.LIMO,'con1_files.txt'))
            writecell(con2_files,fullfile(files.LIMO,'con2_files.txt'))

            clear data
            for N=size(con1_files,1):-1:1
                data{1,N} = con2_files{N};
                data{2,N} = con1_files{N};
            end
            limo_random_select('paired t-test',AvgChanlocs,...
                'LIMOfiles',data, 'analysis_type',...
                'Full scalp analysis', 'type','Channels','nboot',nboot,'tfce',tfce);
            limo_get_effect_size('Paired_Samples_Ttest_parameter_2_1.mat')

            % Param avg (use limo_add_plots to visualize)
            limo_central_tendency_and_ci(fullfile(files.LIMO,'con1_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_distractors')
            limo_central_tendency_and_ci(fullfile(files.LIMO,'con2_files.txt'),...
                1, AvgChanlocs, 'mean', 'Trimmed mean', [],'Con_targets')
            Diff = limo_plot_difference('Con_targets_single_subjects_mean.mat',...
                'Con_distractors_single_subjects_mean.mat',...
                'type','paired','fig',0,'name','Con_diff');
            save('Parameter_difference','Diff')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_P3_P3_GLM_Channels_Time_' estimation{est} '.txt']),...
                'con_1', AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_distractors')
            limo_central_tendency_and_ci(fullfile(files.LIMO,['LIMO_files_P3_P3_GLM_Channels_Time_' estimation{est} '.txt']),...
                'con_2', AvgChanlocs, 'Weighted mean', 'Trimmed mean', [], 'ERPs_targets')
            Diff = limo_plot_difference('ERPs_targets_single_subjects_Weighted mean.mat',...
                'ERPs_distractors_single_subjects_Weighted mean.mat',...
                'type','paired','fig',0,'name','ERP_diff');
            save('ERP_difference','Diff')
        end
    end
    clear STUDY ALLEEG EEG
end
cd(current_folder)