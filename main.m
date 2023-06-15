%% Quadratic-regularized least-squares on real data
%% Initialize
clear all
close all
if ~ismac
    fprintf('\n')
    fprintf('Setting MaxNumCompThreads to 5')
    fprintf('\n')
    maxNumCompThreads(4)
end
addpath redblue
addpath DrosteEffect-BrewerMap-5b84f95/
addpath bluewhitered/
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Set1'))

% Make plots tabs in figure window
set(0,'DefaultFigureWindowStyle','docked')

%% User Parameters
overall_data_folder     = 'data_folder';
sub_data_folder         = 'raw_data';
n                       = 500;     % length of filter for breathing
s                       = 5750;    % length of filter for stimulus
%trial_select            = '3CHO';   % scalar or vector of trial(s); or, odorant name
% Trials just below are the full set of odorant/stimulus trials
%trial_select            = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,40,41,42,43,44,45];
lambda_breathing        = [0 (10.^(1:6))];       % lambda for breathing in optimization
lambda_stimulus         = [(1e-10) (10.^(1:6))];       % lambda for stimulus in optimization
dff_flag                = 1;       % Flag to indicate whether to do deltaF/F
roi_solve_flag          = 0;
roi_glo_filtering_flag  = 1;  % this will be calibrated for n/s longer, both for stim & non-stim
roi_glo_zscore_flag     = 0;
filter_stimulus_flag    = 0; % 1-gaussian, 2-triangular
num_ROIs_to_use         = 20;
stim_extend_flag        = 3; % 1-extend by 2sec, 2-shorten to just one sample (2ms),
                             % 3-same as 2, but moved to peak of breath
                             % after onset of stim
trial_blacklist         = [1 5 9 10 15 20 40 41 42 43 44 45]; % Blacklist of trials to exclude if odorant
                              %     name is used for

if dff_flag==0
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
    fprintf('Not performing deltaF/F0 processing')
    fprintf('\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
end
if roi_glo_filtering_flag==0
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
    fprintf('Not performing FILTERING ON ROI/GLO')
    fprintf('\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
end
if roi_glo_zscore_flag==0
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
    fprintf('Not performing ZSCORE ON ROI/GLO')
    fprintf('\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
end
if filter_stimulus_flag~=0
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
    fprintf('Stimulus is being filtered')
    fprintf('\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')
end


% Flags for individual model fitting
run_individual_fit = 0; % Flag on whether or not to learn filters for trials

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Breathing
fc = 10;
fs = 500;

%%%% ROIs/Glomeruli
filter_select = 'Lowpass'; % 'Lowpass' or 'Bandpass', filtering for data (not breathing)
fl = 1.5;
fh = 70; % will act as cutoff frequency if 'Lowpass' selected
% fl = 0.5; % Usual settings
% fh = 20;
dff_stim_mean_start_t = 1500; % 3sec -> this will be calibrated for n/s longer
dff_stim_mean_end_t = 2000;   % 4sec

%%%% Stimulus
gausswin_alpha_stim = 5;
gausswin_N_stim = 4000;

% Plotting parameters
time_range_to_plot       = [40 43]; % in seconds
colorbar_m = [];
%c_axis = [-1 1];
c_axis = [];      % If c_axis unspecified/empty, normalize c_axis for trial
plot_roi_flag = 0; % Plot roi traces
plot_glo_flag = 0; % Plot glo traces
ignore_raw_for_c_normalization = 1; % Set to 1 if not including raw in normalization of corr

% For graphs with nodes/edges
weight_threshold = 0.04; % Determines which weights to care about
%weight_threshold = 0;

%% Load Data

% Load trial stimulus information
data_excel_doc = readtable('2017-10-27_GAD2-cre-tdTomato+AAV1.Syn.GCaMP6f.xlsx','Format','auto');
data_excel_doc(:,[1 3:6 8:10]) = [];
trial_stim_label = data_excel_doc(contains(data_excel_doc.Var2,'roi_1027_0'),:);
[stim_label,~,stim_index] = unique(trial_stim_label.Var7);

% Check that user has defined trial_select - if not, ask user for trial_select
if ~exist('trial_select','var')
    trial_select = input('Please define trial_select: ','s');
end

% If trial_select is a string, then select trials that have that odorant
if isstr(trial_select)
    odorant_name = trial_select;
    trial_select = find(strcmp(trial_stim_label.Var7,trial_select));
    for ii=1:length(trial_blacklist)
        if ~isempty(find(trial_blacklist(ii)==trial_select))
            fprintf('Removing trial %i\n', trial_blacklist(ii))
            trial_select(find(trial_blacklist(ii)==trial_select))=[];
        end
    end
else
    odorant_name = num2str(trial_select);
end

% Plot stimulus information
figure('NumberTitle', 'off', 'Name','Stim Dist')
    h_hist = histogram(stim_index);
    xticks(1:max(stim_index))
    xticklabels(stim_label)
    set(gca,'FontSize',16)
    xlabel('Stimulus','FontWeight','bold')
    ylabel('# Trials','FontWeight','bold')
    title('Distribution of Stimuli Among Trials','FontWeight','bold','FontSize', 20)
    box off
    h_hist.FaceAlpha = 1;

% Trial-specific data
% Initialize cells so shape is specified
num_trials = length(trial_select);
trial_str  = cell(num_trials,1);
T          = cell(num_trials,1);

% Determine starting index for regression (i.e. when observations start,
% y_i)
b_starting_index = max([n,s])+1;

% Loop through trials
for ii=1:num_trials
    trial_str{ii} = strcat(overall_data_folder,'/',...
                           sub_data_folder,'/roi_1027_', ...
                           num2str(trial_select(ii), '%03.f'), '.csv');
    T{ii}         = readtable(trial_str{ii});
    if (b_starting_index>2000) % If b_starting>4sec
        num_timepoints_to_pad = b_starting_index-2000; % How many zeros to add
        T{ii}((num_timepoints_to_pad+1):(end+num_timepoints_to_pad),:) = ...
            T{ii}; % Shift real values
        T{ii}{1:num_timepoints_to_pad,:} = 0;  % Preprend zeros
    else
        num_timepoints_to_pad = 0;
    end
    original_number_of_ROIs = length(find(cellfun(@(x) ...
            contains(x,'roi'),T{ii}.Properties.VariableNames))')
    T{ii}(:,1+((1+num_ROIs_to_use):original_number_of_ROIs)) = [];
end

%% Show Stimuli for Selected Trials
figure('NumberTitle', 'off', 'Name','Stim for Trials')
    ylim([min(trial_select)-0.5 max(trial_select)+0.5])
    for tt=1:num_trials
        text(0.03,trial_select(tt),trial_stim_label.Var7{trial_select(tt)},...
             'FontSize',20)
    end
    set(gca,'FontSize',20)
    xticks([])
    save_yticks = yticks;
    yticks(min(trial_select):max(trial_select))
    set(gca,'YDir','reverse')

%% Extract Waveforms

% Initialize cells so shape is specified
breathing = cell(num_trials,1);
stimulus  = cell(num_trials,1);

T_roi_indices = cell(num_trials,1);
roi_names     = cell(num_trials,1);
raw_roi_data  = cell(num_trials,1);
raw_glo_data  = cell(num_trials,1);

% Loop through trials
for ii=1:num_trials

    % Breathing/stim
    breathing{ii} = T{ii}.breath;
    stimulus{ii}  = double(T{ii}.stim>0.2);

    % ROIs
    T_roi_indices{ii} = find(cellfun(@(x) contains(x,'roi'),T{ii}.Properties.VariableNames))';
    roi_names{ii}     = T{ii}.Properties.VariableNames(T_roi_indices{ii})';
    raw_roi_data{ii}  = nan(length(roi_names{ii}),length(breathing{ii}));
    for jj = 1:length(roi_names{ii})        % ROIs
        current_roi             = str2double(roi_names{ii}{jj}(4:5));
        raw_roi_data{ii}(current_roi,:) = T{ii}.(roi_names{ii}{jj});
    end

    % Glomeruli
    raw_glo_data{ii} = nan(length(roi_names{ii})/2,length(breathing{ii}));
    for jj = 1:2:length(roi_names{ii})    % Glomeruli
        raw_glo_data{ii}(round(jj/2),:) = mean(raw_roi_data{ii}(jj:(jj+1),:),1);
    end

    % 1: Extend stim length
    if stim_extend_flag==1
        stim_nonzero_indices = find(stimulus{ii}>0);
        final_stim_nonzero_index = max(stim_nonzero_indices);
        stimulus{ii}(final_stim_nonzero_index:(final_stim_nonzero_index+1000)) = 1;
    % 2: Set it so stim only has 1 point
    elseif stim_extend_flag==2
        stim_nonzero_indices = find(stimulus{ii}>0);
        stimulus{ii}(stim_nonzero_indices(2:end))=0;
    % 3: Set it so stim only has 1 point at peak of breathing waveform
    elseif stim_extend_flag==3
        stim_nonzero_indices = find(stimulus{ii}>0);
        stimulus{ii}(stim_nonzero_indices(2:end))=0;

        [PKS,LOCS] = findpeaks(...
            breathing{ii}(stim_nonzero_indices(1):(stim_nonzero_indices+500)),...
            'MinPeakDistance',125,'MinPeakHeight',0.005);
        LOCS(LOCS<50) = [];
        if trial_select(ii)==13 % Trial 13: Have to choose next peak
            br_pk_LOC = LOCS(2)-1+stim_nonzero_indices(1);
        elseif trial_select(ii)==14 % Trial 14: Have to choose next peak
            br_pk_LOC = LOCS(2)-1+stim_nonzero_indices(1);
        else
            br_pk_LOC = LOCS(1)-1+stim_nonzero_indices(1);
        end

        figure('NumberTitle', 'off', 'Name',sprintf('Trial %i', ii))
        plot(raw_glo_data{ii}')
        hold on
        plot(stimulus{ii},'LineWidth',3)
        plot(0.22*zscore(breathing{ii})-0.8,'LineWidth',3)
            line([br_pk_LOC br_pk_LOC],[-1.5 1.5],...
                 'LineWidth',3)
        hold off
        title(sprintf('Trial %i', ii))
        set(gca,'FontSize',16)
        xlim([1800 2800])

        % Set stimulus location to peak breathing
        stimulus{ii}(stim_nonzero_indices(1))=0; % remove old marker
        stimulus{ii}(br_pk_LOC)=1;               % set new marker to pk

        % Add to plot
        hold on
        plot(stimulus{ii}-1.1,'LineWidth',3)
        hold off

    end
end


%% Time
time_generic          = 0.002:0.002:10000;

% Initialize cells so shape is specified
time                  = cell(num_trials,1);
breathing_filter_time = cell(num_trials,1);
stimulus_filter_time  = cell(num_trials,1);
later_time            = cell(num_trials,1);

% Loop through trials
for ii=1:num_trials
    time{ii}                  = time_generic(1:length(breathing{ii}));
    breathing_filter_time{ii} = time{ii}(1:n);
    stimulus_filter_time{ii}  = time{ii}(1:s);
    % later_time - corresponds to b, A*beta, & S*alpha
    later_time{ii}            = time{ii}(max([n,s])+1:end);
end


%% DeltaF/F0
% Initialize cells so shape is specified
dff_roi_data = cell(num_trials,1);
dff_glo_data = cell(num_trials,1);

% Loop through trials
for tt=1:num_trials
    % ROIs
    dff_roi_data{tt} = nan(size(raw_roi_data{tt}));
    for ii=1:size(raw_roi_data{tt},1)
        if dff_flag==1
            if isempty(find(stimulus{tt}>0)) % No stimulus
                dff_roi_data{tt}(ii,:) = (raw_roi_data{tt}(ii,:)/...
                                        mean(raw_roi_data{tt}(ii,(num_timepoints_to_pad+1):end)))-1;
            else % Stimulus
                dff_roi_data{tt}(ii,:) = (raw_roi_data{tt}(ii,:)/...
                    mean(raw_roi_data{tt}(ii,num_timepoints_to_pad+...
                        (dff_stim_mean_start_t:dff_stim_mean_end_t))))-1;
            end
        else
            dff_roi_data{tt}(ii,:) = raw_roi_data{tt}(ii,:);
        end
    end
    % Glomeruli
    dff_glo_data{tt} = nan(size(raw_glo_data{tt}));
    for ii=1:size(raw_glo_data{tt},1)
        if dff_flag==1
            if isempty(find(stimulus{tt}>0)) % No stimulus
                fprintf('Trial %i: No Stimulus\n',trial_select(tt))
                dff_glo_data{tt}(ii,:) = (raw_glo_data{tt}(ii,:)/...
                                        mean(raw_glo_data{tt}(ii,(num_timepoints_to_pad+1):end)))-1;
            else % Stimulus
                fprintf('Trial %i: Stimulus\n',trial_select(tt))
                dff_glo_data{tt}(ii,:) = (raw_glo_data{tt}(ii,:)/...
                    mean(raw_glo_data{tt}(ii,num_timepoints_to_pad+...
                    (dff_stim_mean_start_t:dff_stim_mean_end_t))))-1;
            end
        else
            dff_glo_data{tt}(ii,:) = raw_glo_data{tt}(ii,:);
        end
    end
end


%% Low Pass Filter Breathing + Z-score
[b_filt,a_filt]    = butter(4,fc/(fs/2));
for tt=1:num_trials
    breathing{tt} = filtfilt(b_filt, a_filt, breathing{tt});
    breathing{tt} = zscore(breathing{tt}); % zscore
end

% Plot filter
% Uncomment below to plot breathing filter
% lpf_fv = fvtool(b_filt, a_filt);
% set(lpf_fv, 'MagnitudeDisplay', 'Magnitude squared')
% lpf_fv.Fs = fs;
% set(gca,'FontSize',16)
% xlim([0 50])
% objects = findall(gca);
% objects(2).LineWidth = 3;

%% Low Pass Filter Stimulus
% First plot stimulus
h_stimulus = figure('NumberTitle', 'off', 'Name','Stimulus');
for tt=1:num_trials
    plot(time{tt},stimulus{tt},'LineWidth',3)
    hold on
end
hold off
xlabel('Time (s)')
ylabel('Stimulus State')
set(gca, 'FontSize', 16)


% Filter
if (filter_stimulus_flag~=0)

    if filter_stimulus_flag==1
        % Apply gaussian window to stimulus
        stim_gausswin = gausswin(gausswin_N_stim,gausswin_alpha_stim);
        for tt=1:num_trials
            stimulus{tt} = conv(stimulus{tt}, stim_gausswin, 'same');
            if ~isempty(find(stimulus{tt}>0)) % If stimulus
                stimulus{tt} = stimulus{tt}/max(stimulus{tt});
            end
        end
    elseif filter_stimulus_flag==2
        % Apply triangular window to stimulus
        for tt=1:num_trials
            if ~isempty(find(stimulus{tt}>0)) % If stimulus
                stimulus{tt}(find(stimulus{tt}>0)) = triang(length(find(stimulus{tt}>0)));
            end
        end
    end

    % Plot smoothed stimulus over raw stimulus
    figure(h_stimulus)
    hold on
    set(gca,'ColorOrderIndex',1)
    for tt=1:num_trials
        plot(time{tt},stimulus{tt},'LineWidth',3)
    end
    hold off
end

%% Band Pass Filter ROI & Glomeruli + Z-score
if strcmp(filter_select,'Bandpass')
    [b_bp, a_bp] = butter(4,[fl fh]/(fs/2));
elseif strcmp(filter_select,'Lowpass')
    [b_bp, a_bp] = butter(4,fh/(fs/2));
end

% Initialize variables
fil_roi_data = cell(num_trials,1);
fil_glo_data = cell(num_trials,1);

% Loop through trials
for tt=1:num_trials

    % ROIs
    fil_roi_data{tt} = nan(size(raw_roi_data{tt}));
    for ii=1:size(raw_roi_data{tt},1)
        if roi_glo_filtering_flag
            fil_roi_data{tt}(ii,:) = ...
                filtfilt(b_bp, a_bp, dff_roi_data{tt}(ii,:));
        else
            fil_roi_data{tt}(ii,:) = dff_roi_data{tt}(ii,:);
        end
    end
    if roi_glo_zscore_flag
        fil_roi_data{tt} = zscore(fil_roi_data{tt}')'; % zscore
    end

    % GLOMERULI
    fil_glo_data{tt} = nan(size(raw_glo_data{tt}));
    for ii=1:size(raw_glo_data{tt},1)
        if roi_glo_filtering_flag
            fil_glo_data{tt}(ii,:) = ...
                filtfilt(b_bp, a_bp, dff_glo_data{tt}(ii,:));
        else
            fil_glo_data{tt}(ii,:) = dff_glo_data{tt}(ii,:);
        end
    end
    if roi_glo_zscore_flag
        fil_glo_data{tt} = zscore(fil_glo_data{tt}')'; % zscore
    end
end

% Plot filter
% Uncomment below to plot filter
% bp_fv = fvtool(b_bp, a_bp);
% set(bp_fv, 'MagnitudeDisplay', 'Magnitude squared')
% bp_fv.Fs = fs;
% set(gca,'FontSize',16)
% xlim([0 50])
% title('Filtering of ROIs/Glomeruli')
% objects = findall(gca);
% objects(2).LineWidth = 3;


%% Check out correlations between ROIs/Glomeruli

if run_individual_fit % only run if flag set

    c_limits_roi = cell(num_trials,1);
    c_limits_glo = cell(num_trials,1);

    % RAW ROIs
    for tt=1:num_trials

        c_limits_roi{tt} = nan([3 2]);

        RHO_raw_roi_data{tt} = corr(raw_roi_data{tt}');

        corr_fig_hdl{tt} = figure('NumberTitle', 'off', 'Name',...
                                  sprintf('Corr - Trial %i',trial_select(tt)));
        subplot(3,2,1)
        imagesc(RHO_raw_roi_data{tt})
        title(sprintf('RHO for raw roi data: Trial %i',trial_select(tt)))
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('roi')
        ylabel('roi')
        if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
            c_limits_roi{tt}(1,:) = caxis;
            c_max = max(max(RHO_raw_roi_data{tt}-diag(diag(RHO_raw_roi_data{tt}))));
            c_limits_roi{tt}(1,2) = c_max;
        else % Use c_axis if given, e.g. c_axis = [-1 1]
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image

    end


    % RAW GLOMERULI
    for tt=1:num_trials

        figure(corr_fig_hdl{tt})

        c_limits_glo{tt} = nan([3 2]);

        RHO_raw_glo_data{tt} = corr(raw_glo_data{tt}');

        %figure
        subplot(3,2,2)
        imagesc(RHO_raw_glo_data{tt})
        title(sprintf('RHO for raw glo data: Trial %i',trial_select(tt)))
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('glomerulus')
        ylabel('glomerulus')
        if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
            c_limits_glo{tt}(1,:) = caxis;
            c_max = max(max(RHO_raw_glo_data{tt}-diag(diag(RHO_raw_glo_data{tt}))));
            c_limits_glo{tt}(1,2) = c_max;
        else % Use c_axis if given, e.g. c_axis = [-1 1]
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image
    end

    % FIL ROIs
    for tt=1:num_trials

        figure(corr_fig_hdl{tt})

        RHO_fil_roi_data{tt} = corr(fil_roi_data{tt}');

        %figure
        subplot(3,2,3)
        imagesc(RHO_fil_roi_data{tt})
        title('RHO for filtered roi data')
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('roi')
        ylabel('roi')
        if isempty(c_axis)
            c_limits_roi{tt}(2,:) = caxis;
            c_max = max(max(RHO_fil_roi_data{tt}-diag(diag(RHO_fil_roi_data{tt}))));
            c_limits_roi{tt}(2,2) = c_max;
        else
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image
    end

    % FIL GLOMERULI
    for tt=1:num_trials

        figure(corr_fig_hdl{tt})

        RHO_fil_glo_data{tt} = corr(fil_glo_data{tt}');

        %figure
        subplot(3,2,4)
        imagesc(RHO_fil_glo_data{tt})
        title('RHO for filtered glomeruli data')
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('glomerulus')
        ylabel('glomerulus')
        if isempty(c_axis)
            c_limits_glo{tt}(2,:) = caxis;
            c_max = max(max(RHO_fil_glo_data{tt}-diag(diag(RHO_fil_glo_data{tt}))));
            c_limits_glo{tt}(2,2) = c_max;
        else
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image
    end

end

%% ALL TRIALS TOGETHER: Plot correlations between ROIs/Glomeruli across trials

% Compute correlations
RHO_all_raw_roi_data = corr(cell2mat(raw_roi_data')');
RHO_all_raw_glo_data = corr(cell2mat(raw_glo_data')');
RHO_all_fil_roi_data = corr(cell2mat(fil_roi_data')');
RHO_all_fil_glo_data = corr(cell2mat(fil_glo_data')');

% RAW ROIs: Plot correlations between all raw ROIs across all trials
all_corr_fig_hdl = figure('NumberTitle', 'off', 'Name','All Corr: Br Adj');
    subplot(3,2,1)
    imagesc(RHO_all_raw_roi_data)
    title('RHO for all raw roi data')
    h_bar = colorbar;
    ylabel(h_bar,'RHO')
    set(gca, 'FontSize', 24)
    xlabel('roi')
    ylabel('roi')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_roi_all_trial(1,:) = caxis;
        c_max = max(max(RHO_all_raw_roi_data-diag(diag(RHO_all_raw_roi_data))));
        c_limits_roi_all_trial(1,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image

% RAW GLOMERULI:
subplot(3,2,2)
    imagesc(RHO_all_raw_glo_data)
    title('RHO for all raw glo data')
    h_bar = colorbar;
    ylabel(h_bar,'RHO')
    set(gca, 'FontSize', 24)
    xlabel('glomerulus')
    ylabel('glomerulus')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_glo_all_trial(1,:) = caxis;
        c_max = max(max(RHO_all_raw_glo_data-diag(diag(RHO_all_raw_glo_data))));
        c_limits_glo_all_trial(1,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image

% FIL ROIs:
subplot(3,2,3)
    imagesc(RHO_all_fil_roi_data)
    title('RHO for all fil roi data')
    h_bar = colorbar;
    ylabel(h_bar,'RHO')
    set(gca, 'FontSize', 24)
    xlabel('roi')
    ylabel('roi')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_roi_all_trial(2,:) = caxis;
        c_max = max(max(RHO_all_fil_roi_data-diag(diag(RHO_all_fil_roi_data))));
        c_limits_roi_all_trial(2,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image

% FIL GLOMERULI:
subplot(3,2,4)
    imagesc(RHO_all_fil_glo_data)
    title('RHO for all fil glo data')
    h_bar = colorbar;
    ylabel(h_bar,'RHO')
    set(gca, 'FontSize', 24)
    xlabel('glomerulus')
    ylabel('glomerulus')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_glo_all_trial(2,:) = caxis;
        c_max = max(max(RHO_all_fil_glo_data-diag(diag(RHO_all_fil_glo_data))));
        c_limits_glo_all_trial(2,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image

%% Compare correlations between fil, br adj, stim adj, and br+stim adj
%  across all trials using filter learned across all trials
all_corr_fig_br_stim_hdl = figure('NumberTitle', 'off', 'Name','All Corr: Br+Stim Adj');
    subplot(2,2,1)
    imagesc(RHO_all_fil_glo_data)
    title('Before Adjustment')
    h_bar = colorbar;
    ylabel(h_bar,'\rho')
    set(gca, 'FontSize', 24)
    xlabel('Glomerulus')
    ylabel('Glomerulus')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_glo_all_trial_br_stim(1,:) = caxis;
        c_max = max(max(RHO_all_fil_glo_data-diag(diag(RHO_all_fil_glo_data))));
        c_limits_glo_all_trial_br_stim(1,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image



%% Data Parameters & Selection

% Loop through trials
for tt=1:num_trials
    m(tt) = min([length(breathing{tt})-n, ...
                 length(breathing{tt})-s]); % m - number of examples

    % Initialize for loop
    A{tt} = nan(m(tt),n); % Breathing
    S{tt} = nan(m(tt),s); % Stimulus
    window_end_index = b_starting_index-1; % End of window preceding b start

    % Loop through windows
    for ii = 1:m(tt)
        A{tt}(ii,:) = breathing{tt}((window_end_index-(n-1)):...
                            window_end_index); % Breathing signal - m x n
        S{tt}(ii,:) =  stimulus{tt}((window_end_index-(s-1)):...
                            window_end_index); % Stimulus signal  - m x s
        window_end_index = window_end_index+1; % Increment window_end_index
    end
end

%% Solve problem

% quad-smoothness
omega_breathing = sparse(toeplitz([2, -1, zeros(1, n-2)]));
omega_breathing_square = omega_breathing*omega_breathing;
omega_stimulus = sparse(toeplitz([2, -1, zeros(1, s-2)]));
omega_stimulus_square = omega_stimulus*omega_stimulus;

% create matrix for smoothness penalty & lambda learning
if (length(lambda_breathing)>1) || (length(lambda_stimulus)>1)
    delta = nan(n+s,n+s,length(lambda_breathing),length(lambda_stimulus));
    for ii=1:length(lambda_breathing)
        for jj=1:length(lambda_stimulus)
            delta(:,:,ii,jj) = [(lambda_breathing(ii)*omega_breathing_square)...
                                zeros(n,s);...
                                zeros(s,n)...
                                (lambda_stimulus(jj)*omega_stimulus_square)];
        end
    end
else
    delta = [(lambda_breathing*omega_breathing_square)  zeros(n,s);
             zeros(s,n)   (lambda_stimulus*omega_stimulus_square)];
end


if run_individual_fit % only run if flag set

    if (length(lambda_breathing)>1) || (length(lambda_stimulus)>1)
        % I didn't implement GCV for individual fit yet, so if a range of
        % parameters is given for lambdas then I will print an error
        error('Need to add GCV support for individual fit')
    end
    for tt=1:num_trials

        % ROIs
        roi_quad_beta{tt}  = nan(size(fil_roi_data{tt},1),n);
        roi_quad_alpha{tt} = nan(size(fil_roi_data{tt},1),s);
        roi_filtered_breathing_quad{tt} = nan(size(fil_roi_data{tt},1),size(A{tt},1));
        roi_filtered_stimulus_quad{tt}  = nan(size(fil_roi_data{tt},1),size(S{tt},1));

        X = [A{tt} S{tt}]; % Number of samples X total filter length (B+S)
        X_term = (((X')*X) + (m(tt)*delta))\(X');

        for ii=1:size(fil_roi_data{tt},1)   % ROIs
            b{tt} = fil_roi_data{tt}(ii,b_starting_index:end);
            x_quad{tt} = X_term*(b{tt}');
            roi_quad_beta{tt}(ii,:)  = x_quad{tt}(1:n); % Breathing filter
            roi_quad_alpha{tt}(ii,:) = x_quad{tt}((n+1):end); % Stimulus filter
            roi_filtered_breathing_quad{tt}(ii,:) = A{tt}*(roi_quad_beta{tt}(ii,:)');
            roi_filtered_stimulus_quad{tt}(ii,:)  = S{tt}*(roi_quad_alpha{tt}(ii,:)');
        end

        % Glomeruli
        glo_quad_beta{tt}  = nan(size(fil_glo_data{tt},1),n);
        glo_quad_alpha{tt} = nan(size(fil_glo_data{tt},1),s);
        glo_filtered_breathing_quad{tt} = nan(size(fil_glo_data{tt},1),size(A{tt},1));
        glo_filtered_stimulus_quad{tt}  = nan(size(fil_glo_data{tt},1),size(S{tt},1));

        for ii=1:size(fil_glo_data{tt},1)   % Glomeruli
            b{tt} = fil_glo_data{tt}(ii,b_starting_index:end);
            x_quad{tt} = X_term*(b{tt}');
            glo_quad_beta{tt}(ii,:)  = x_quad{tt}(1:n); % Breathing filter
            glo_quad_alpha{tt}(ii,:) = x_quad{tt}((n+1):end); % Stimulus filter
            glo_filtered_breathing_quad{tt}(ii,:) = A{tt}*(glo_quad_beta{tt}(ii,:)');
            glo_filtered_stimulus_quad{tt}(ii,:)  = S{tt}*(glo_quad_alpha{tt}(ii,:)');
        end

    end
end

%% Solve problem for ALL data across selected trials
A_all  = cell2mat(A');
S_all  = cell2mat(S');

X      = [A_all S_all];

% ROIs
if roi_solve_flag
    for ii=1:size(fil_roi_data{1},1)   % ROIs
        for tt=1:num_trials
            b_all_pre_cat{tt} = fil_roi_data{tt}(ii,b_starting_index:end);
        end
        b_all = cell2mat(b_all_pre_cat); % concatenate trials
        %X_term = (((X')*X) + (sum(m)*delta))\(X');
        [x_quad_all, roi_ii(ii), roi_jj(ii)] = gcv_solver( X, delta, b_all,...
                                                      lambda_breathing,...
                                                      lambda_stimulus);
        %x_quad_all = X_term*(b_all');
        roi_quad_beta_all(ii,:)  = x_quad_all(1:n); % Breathing filter
        roi_quad_alpha_all(ii,:) = x_quad_all((n+1):end); % Stimulus filter
        for tt=1:num_trials
            % Filter
            roi_filtered_breathing_quad_all{tt}(ii,:) = A{tt}*(roi_quad_beta_all(ii,:)');
            roi_filtered_stimulus_quad_all{tt}(ii,:)  = S{tt}*(roi_quad_alpha_all(ii,:)');
        end
    end
end

% Glomeruli
for ii=1:size(fil_glo_data{1},1)   % Glomeruli
    fprintf('Glomerulus %i\n',ii)
    for tt=1:num_trials
        b_all_pre_cat{tt} = fil_glo_data{tt}(ii,b_starting_index:end);
    end
    b_all = cell2mat(b_all_pre_cat); % concatenate trials
    [x_quad_all, glo_ii(ii), glo_jj(ii)] = gcv_solver( X, delta, b_all,...
                                                      lambda_breathing,...
                                                      lambda_stimulus );
    %x_quad_all = X_term*(b_all');
    glo_quad_beta_all(ii,:)  = x_quad_all(1:n); % Breathing filter
    glo_quad_alpha_all(ii,:) = x_quad_all((n+1):end); % Stimulus filter
    for tt=1:num_trials
        % Filter
        glo_filtered_breathing_quad_all{tt}(ii,:) = A{tt}*(glo_quad_beta_all(ii,:)');
        glo_filtered_stimulus_quad_all{tt}(ii,:)  = S{tt}*(glo_quad_alpha_all(ii,:)');
    end
end

%% SAVE INTERMEDIATE FILE
save([odorant_name '_all_INTERMEDIATE.mat'])

%% PLOT LEARNED FILTERS FOR TRIALS

if run_individual_fit % only run if flag set
    % ROIS
    % Breathing
    for tt=1:num_trials
        figure('NumberTitle', 'off', 'Name',...
                                  sprintf('ROI Filters - Trial %i',trial_select(tt)));
        subplot(2,1,1)
        plot(breathing_filter_time{tt}, zeros(size(roi_quad_beta{tt}(1,:))),'--')
        for ii=1:size(roi_quad_beta{tt},1)
            hold on
            plot(breathing_filter_time{tt},flip(roi_quad_beta{tt}(ii,:)), 'LineWidth', 3)
        end
        hold off
        xlim([min(breathing_filter_time{tt}) max(breathing_filter_time{tt})])
        title('Breathing Filters for ROIs')
        xlabel('Reverse time (sec)')
        set(gca, 'FontSize', 16)

        % Stimulus
        subplot(2,1,2)
        plot(stimulus_filter_time{tt}, zeros(size(roi_quad_alpha{tt}(1,:))),'--')
        for ii=1:size(roi_quad_alpha{tt},1)
            hold on
            plot(stimulus_filter_time{tt},flip(roi_quad_alpha{tt}(ii,:)), 'LineWidth', 3)
        end
        hold off
        xlim([min(stimulus_filter_time{tt}) max(stimulus_filter_time{tt})])
        title('Stimulus Filters for ROIs')
        xlabel('Reverse time (sec)')
        set(gca, 'FontSize', 16)

    end

    for tt=1:num_trials
        % GLOMERULI
        % Breathing
        figure('NumberTitle', 'off', 'Name',...
                                  sprintf('GLO Filters - Trial %i',trial_select(tt)));
        subplot(2,1,1)
        plot(breathing_filter_time{tt}, zeros(size(glo_quad_beta{tt}(1,:))),'--')
        for ii=1:size(glo_quad_beta{tt},1)
            hold on
            plot(breathing_filter_time{tt},flip(glo_quad_beta{tt}(ii,:)), 'LineWidth', 3)
        end
        hold off
        xlim([min(breathing_filter_time{tt}) max(breathing_filter_time{tt})])
        title('Breathing Filters for Glomeruli')
        xlabel('Reverse time (sec)')
        set(gca, 'FontSize', 16)

        % Stimulus
        subplot(2,1,2)
        plot(stimulus_filter_time{tt}, zeros(size(glo_quad_alpha{tt}(1,:))),'--')
        for ii=1:size(glo_quad_alpha{tt},1)
            hold on
            plot(stimulus_filter_time{tt},flip(glo_quad_alpha{tt}(ii,:)), 'LineWidth', 3)
        end
        hold off
        xlim([min(stimulus_filter_time{tt}) max(stimulus_filter_time{tt})])
        title('Stimulus Filters for Glomeruli')
        xlabel('Reverse time (sec)')
        set(gca, 'FontSize', 16)

    end
end

%% ALL TRIALS: Plot filters learned across all trials

if roi_solve_flag
    % ROIS
    % Breathing
    figure('NumberTitle', 'off', 'Name','All ROI Filters');
    subplot(2,1,1)
    plot(breathing_filter_time{1}, zeros(size(breathing_filter_time{1})),'--')
    for ii=1:size(roi_quad_beta_all,1) % Loop through ROIs
        hold on
        plot(breathing_filter_time{1},flip(roi_quad_beta_all(ii,:)), 'LineWidth', 3)
    end
    hold off
    xlim([min(breathing_filter_time{1}) max(breathing_filter_time{1})])
    title('Breathing Filters for ROIs')
    xlabel('Reverse time (sec)')
    set(gca, 'FontSize', 16)

    % Stimulus
    subplot(2,1,2)
    plot(stimulus_filter_time{1}, zeros(size(stimulus_filter_time{1})),'--')
    for ii=1:size(roi_quad_alpha_all,1) % Loop through ROIs
        hold on
        plot(stimulus_filter_time{1},flip(roi_quad_alpha_all(ii,:)), 'LineWidth', 3)
    end
    hold off
    xlim([min(stimulus_filter_time{1}) max(stimulus_filter_time{1})])
    title('Stimulus Filters for ROIs')
    xlabel('Reverse time (sec)')
    set(gca, 'FontSize', 16)
end


% Glomeruli
% Breathing
figure('NumberTitle', 'off', 'Name','All GLO Filters');
subplot(2,1,1)
%plot(breathing_filter_time{1}, zeros(size(breathing_filter_time{1})),'--')
% Create string for legend
% Plot synthetic filters also if simulation
leg_index = 1;
if contains(sub_data_folder, 'simulated')
    load([overall_data_folder '/parameters/' sub_data_folder '.mat'])
    plot(breathing_filter_time{1}, breathing_fil_weights',...
     'LineWidth',3)
    hold on
    set(gca,'ColorOrderIndex',1)
    for ii=1:size(glo_quad_beta_all,1)
        glom_leg_str{ii} = sprintf('G%i True',ii);
        leg_index = leg_index+1;
    end
end
% Plot filters
for ii=1:size(glo_quad_beta_all,1) % Loop through glomeruli
    %plot(breathing_filter_time{1},flip(glo_quad_beta_all(ii,:)),'--', ...
     %    'LineWidth', 3)
    if ismac
       plot(breathing_filter_time{1},flip(glo_quad_beta_all(ii,:)),'--', ...
            'LineWidth', 3)
    else
       plot(breathing_filter_time{1},flip(glo_quad_beta_all(ii,:)), ...
            'LineWidth', 3)
    end
    hold on
    glom_leg_str{leg_index} = sprintf('G%i Estimated',ii);
    leg_index = leg_index+1;
end
hold off
xlim([min(breathing_filter_time{1}) max(breathing_filter_time{1})])
title('Breathing Filters for Glomeruli')
xlabel('Reverse time (sec)')
set(gca, 'FontSize', 16)
legend(glom_leg_str)

% Stimulus
subplot(2,1,2)
%plot(stimulus_filter_time{1}, zeros(size(stimulus_filter_time{1})),'--')
% Plot synthetic filters also if simulation
if contains(sub_data_folder, 'simulated')
    plot(stimulus_filter_time{1}, stimulus_fil_weights',...
     'LineWidth',3)
    hold on
    set(gca,'ColorOrderIndex',1)
end
% Plot filters
for ii=1:size(glo_quad_alpha_all,1) % Loop through glomeruli
    %plot(stimulus_filter_time{1},flip(glo_quad_alpha_all(ii,:)),'--', ...
     %    'LineWidth', 3)
    plot(stimulus_filter_time{1},flip(glo_quad_alpha_all(ii,:)), ...
         'LineWidth', 3)
    hold on
end
hold off
xlim([min(stimulus_filter_time{1}) max(stimulus_filter_time{1})])
title('Stimulus Filters for Glomeruli')
xlabel('Reverse time (sec)')
set(gca, 'FontSize', 16)
legend(glom_leg_str)

%% PLOT ALL FILTERS LEARNED ACROSS ALL TRIALS AS SEPARATE LINES
figure('NumberTitle', 'off', 'Name','All GLO Filters (lines)');
for ii=1:size(glo_quad_alpha_all,1) % Loop through glomeruli
    plot(stimulus_filter_time{1},(flip(glo_quad_alpha_all(ii,:))/...
         max(glo_quad_alpha_all(ii,:)))+ii, ...
         'LineWidth', 3)
    hold on
end
hold off

% Plot settings
xlim([min(stimulus_filter_time{1}) max(stimulus_filter_time{1})])
title('Stimulus Filters for Glomeruli')
xlabel('Reverse time (sec)')
set(gca, 'FontSize', 16)
legend(glom_leg_str)
ylabel('Glomerulus')
box off
yticks(1:size(glo_quad_alpha_all,1))


%% BREATHING & STIMULUS ADJUSTED ROIs/GLOMERULI

if run_individual_fit % only run if flag set
    for tt=1:num_trials

        % ROIs
        roi_breathing_adjusted{tt} = nan(size(fil_roi_data{tt}(:,b_starting_index:end)));
        roi_stimulus_adjusted{tt}  = nan(size(fil_roi_data{tt}(:,b_starting_index:end)));
        roi_both_adjusted{tt}      = nan(size(fil_roi_data{tt}(:,b_starting_index:end)));
        for ii=1:size(fil_roi_data{tt},1)   % ROIs
            roi_breathing_adjusted{tt}(ii,:) = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_breathing_quad{tt}(ii,:);
            roi_stimulus_adjusted{tt}(ii,:)  = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_stimulus_quad{tt}(ii,:);
            roi_both_adjusted{tt}(ii,:)      = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_stimulus_quad{tt}(ii,:)...
                                           - roi_filtered_breathing_quad{tt}(ii,:);
        end

        % Glomeruli
        glo_breathing_adjusted{tt} = nan(size(fil_glo_data{tt}(:,b_starting_index:end)));
        glo_stimulus_adjusted{tt}  = nan(size(fil_glo_data{tt}(:,b_starting_index:end)));
        glo_both_adjusted{tt}      = nan(size(fil_glo_data{tt}(:,b_starting_index:end)));
        for ii=1:size(fil_glo_data{tt},1)
            glo_breathing_adjusted{tt}(ii,:) = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                           - glo_filtered_breathing_quad{tt}(ii,:);
            glo_stimulus_adjusted{tt}(ii,:)  = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                           - glo_filtered_stimulus_quad{tt}(ii,:);
            glo_both_adjusted{tt}(ii,:)      = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                           - glo_filtered_stimulus_quad{tt}(ii,:)...
                                           - glo_filtered_breathing_quad{tt}(ii,:);
        end

    end
end

%% ALL: BREATHING & STIMULUS ADJUSTED ROIs/GLOMERULI
for tt=1:num_trials

    if roi_solve_flag
        % ROIs
        for ii=1:size(fil_roi_data{tt},1)   % Loop through ROIs
            roi_breathing_adjusted_all{tt}(ii,:) = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_breathing_quad_all{tt}(ii,:);
            roi_stimulus_adjusted_all{tt}(ii,:)  = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_stimulus_quad_all{tt}(ii,:);
            roi_both_adjusted_all{tt}(ii,:)      = fil_roi_data{tt}(ii,b_starting_index:end) ...
                                           - roi_filtered_stimulus_quad_all{tt}(ii,:)...
                                           - roi_filtered_breathing_quad_all{tt}(ii,:);
        end
    end

    % Glomeruli
    for ii=1:size(fil_glo_data{tt},1)  % Loop through glomeruli
        glo_breathing_adjusted_all{tt}(ii,:) = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                       - glo_filtered_breathing_quad_all{tt}(ii,:);
        glo_stimulus_adjusted_all{tt}(ii,:)  = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                       - glo_filtered_stimulus_quad_all{tt}(ii,:);
        glo_both_adjusted_all{tt}(ii,:)      = fil_glo_data{tt}(ii,b_starting_index:end) ...
                                       - glo_filtered_stimulus_quad_all{tt}(ii,:)...
                                       - glo_filtered_breathing_quad_all{tt}(ii,:);
    end

end



%% COMPUTE FVU

for tt=1:num_trials

    if roi_solve_flag
        % ROIs
        % Total variance
        var_roi_total{tt}         = var(fil_roi_data{tt}(:,b_starting_index:end),[],2);
        if run_individual_fit % only run if flag set
            % Variances
            var_roi_breathing_adj{tt} = var(roi_breathing_adjusted{tt},[],2);
            var_roi_stimulus_adj{tt}  = var(roi_stimulus_adjusted{tt},[],2);
            var_roi_both_adj{tt}      = var(roi_both_adjusted{tt},[],2);

            % FVUs
            fvu_roi_breathing_adj{tt} = var_roi_breathing_adj{tt}./var_roi_total{tt};
            fvu_roi_stimulus_adj{tt}  = var_roi_stimulus_adj{tt}./var_roi_total{tt};
            fvu_roi_both_adj{tt}      = var_roi_both_adj{tt}./var_roi_total{tt};
        end
    end

    % Glomeruli
    % Total variance
    var_glo_total{tt}         = var(fil_glo_data{tt}(:,b_starting_index:end),[],2);
    if run_individual_fit % only run if flag set
        % Variances
        var_glo_breathing_adj{tt} = var(glo_breathing_adjusted{tt},[],2);
        var_glo_stimulus_adj{tt}  = var(glo_stimulus_adjusted{tt},[],2);
        var_glo_both_adj{tt}      = var(glo_both_adjusted{tt},[],2);

        % FVUs
        fvu_glo_breathing_adj{tt} = var_glo_breathing_adj{tt}./var_glo_total{tt};
        fvu_glo_stimulus_adj{tt}  = var_glo_stimulus_adj{tt}./var_glo_total{tt};
        fvu_glo_both_adj{tt}      = var_glo_both_adj{tt}./var_glo_total{tt};
    end

end


%% ALL: COMPUTE FVU
for tt=1:num_trials

    if roi_solve_flag
        % ROIs
        % Variances
        var_roi_breathing_adj_all{tt} = var(roi_breathing_adjusted_all{tt},[],2);
        var_roi_stimulus_adj_all{tt}  = var(roi_stimulus_adjusted_all{tt},[],2);
        var_roi_both_adj_all{tt}      = var(roi_both_adjusted_all{tt},[],2);

        % FVUs
        fvu_roi_breathing_adj_all{tt} = var_roi_breathing_adj_all{tt}./var_roi_total{tt};
        fvu_roi_stimulus_adj_all{tt}  = var_roi_stimulus_adj_all{tt}./var_roi_total{tt};
        fvu_roi_both_adj_all{tt}      = var_roi_both_adj_all{tt}./var_roi_total{tt};
    end

    % Glomeruli
    % Variances
    var_glo_breathing_adj_all{tt} = var(glo_breathing_adjusted_all{tt},[],2);
    var_glo_stimulus_adj_all{tt}  = var(glo_stimulus_adjusted_all{tt},[],2);
    var_glo_both_adj_all{tt}      = var(glo_both_adjusted_all{tt},[],2);

    % FVUs
    fvu_glo_breathing_adj_all{tt} = var_glo_breathing_adj_all{tt}./var_glo_total{tt};
    fvu_glo_stimulus_adj_all{tt}  = var_glo_stimulus_adj_all{tt}./var_glo_total{tt};
    fvu_glo_both_adj_all{tt}      = var_glo_both_adj_all{tt}./var_glo_total{tt};

end


%% Also compute across all trials to get a general FVU

if roi_solve_flag
    % ROIs
    % Extract relevant fil_roi_data values
    for tt=1:num_trials
        fil_roi_data_b_start{tt} = fil_roi_data{tt}(:,b_starting_index:end);
    end

    % Concatenate trials
    fil_roi_data_b_start_cat       = cell2mat(fil_roi_data_b_start);
    roi_breathing_adjusted_all_cat = cell2mat(roi_breathing_adjusted_all);
    roi_stimulus_adjusted_all_cat  = cell2mat(roi_stimulus_adjusted_all);
    roi_both_adjusted_all_cat      = cell2mat(roi_both_adjusted_all);

    % Variances
    var_roi_total_all_gen         = var(fil_roi_data_b_start_cat,[],2);
    var_roi_breathing_adj_all_gen = var(roi_breathing_adjusted_all_cat,[],2);
    var_roi_stimulus_adj_all_gen  = var(roi_stimulus_adjusted_all_cat,[],2);
    var_roi_both_adj_all_gen      = var(roi_both_adjusted_all_cat,[],2);

    % FVUs
    fvu_roi_breathing_adj_all_gen = var_roi_breathing_adj_all_gen./var_roi_total_all_gen;
    fvu_roi_stimulus_adj_all_gen  = var_roi_stimulus_adj_all_gen./var_roi_total_all_gen;
    fvu_roi_both_adj_all_gen      = var_roi_both_adj_all_gen./var_roi_total_all_gen;
end


% Glomeruli
% Extract relevant fil_glo_data values
for tt=1:num_trials
    fil_glo_data_b_start{tt} = fil_glo_data{tt}(:,b_starting_index:end);
end

% Concatenate trials
fil_glo_data_b_start_cat       = cell2mat(fil_glo_data_b_start);
glo_breathing_adjusted_all_cat = cell2mat(glo_breathing_adjusted_all);
glo_stimulus_adjusted_all_cat  = cell2mat(glo_stimulus_adjusted_all);
glo_both_adjusted_all_cat      = cell2mat(glo_both_adjusted_all);

% Variances
var_glo_total_all_gen         = var(fil_glo_data_b_start_cat,[],2);
var_glo_breathing_adj_all_gen = var(glo_breathing_adjusted_all_cat,[],2);
var_glo_stimulus_adj_all_gen  = var(glo_stimulus_adjusted_all_cat,[],2);
var_glo_both_adj_all_gen      = var(glo_both_adjusted_all_cat,[],2);

% FVUs
fvu_glo_breathing_adj_all_gen = var_glo_breathing_adj_all_gen./var_glo_total_all_gen;
fvu_glo_stimulus_adj_all_gen  = var_glo_stimulus_adj_all_gen./var_glo_total_all_gen;
fvu_glo_both_adj_all_gen      = var_glo_both_adj_all_gen./var_glo_total_all_gen;


%% PLOT FVUs

if run_individual_fit % Only run if flag set
    % Colors
    color_order = get(gca, 'ColorOrder');
    col_to_use = color_order(1,:);

    for tt=1:num_trials

        % Create figure
        figure('NumberTitle', 'off', 'Name', ...
               sprintf('FVU: Breathing - Trial %i',trial_select(tt)));

        if roi_solve_flag
            % ROIs
            subplot(2,1,1)
            bar(fvu_roi_breathing_adj{tt}, 'FaceColor', col_to_use)
            title(sprintf('FVUs for Breathing-Adjusted ROIs - Trial %i',...
                  trial_select(tt)))
            xlabel('ROI')
            ylabel('FVU')
            set(gca,'FontSize', 16)
            box off
            xlim([0.5 size(fil_roi_data{tt},1)+0.5])
        end

        % Glomeruli
        subplot(2,1,2)
        bar(fvu_glo_breathing_adj{tt}, 'FaceColor', col_to_use)
        title(sprintf('FVUs for Breathing-Adjusted Glomeruli - Trial %i',...
              trial_select(tt)))
        xlabel('Glomerulus')
        ylabel('FVU')
        set(gca,'FontSize', 16)
        box off
        xlim([0.5 size(fil_glo_data{tt},1)+0.5])

    end
end


%% ALL: PLOT FVUs


% Colors
color_order = get(gca, 'ColorOrder');
col_to_use = color_order(1,:);

for tt=1:num_trials

    % Create figure
    figure('NumberTitle', 'off', 'Name', ...
           sprintf('All FVU: Breathing - Trial %i',trial_select(tt)));

    if roi_solve_flag
        % ROIs
        subplot(2,1,1)
        bar(fvu_roi_breathing_adj_all{tt}, 'FaceColor', col_to_use)
        title(sprintf('FVUs for Breathing-Adjusted ROIs - Trial %i',...
              trial_select(tt)))
        xlabel('ROI')
        ylabel('FVU')
        set(gca,'FontSize', 16)
        box off
        xlim([0.5 size(fil_roi_data{tt},1)+0.5])
    end

    % Glomeruli
    subplot(2,1,2)
    bar(fvu_glo_breathing_adj_all{tt}, 'FaceColor', col_to_use)
    title(sprintf('FVUs for Breathing-Adjusted Glomeruli - Trial %i',...
          trial_select(tt)))
    xlabel('Glomerulus')
    ylabel('FVU')
    set(gca,'FontSize', 16)
    box off
    xlim([0.5 size(fil_glo_data{tt},1)+0.5])

end

%% ALL ACROSS TRIALS: PLOT FVUs ACROSS TRIALS FOR BREATHING

% Colors
color_order = get(gca, 'ColorOrder');
col_to_use = color_order(1,:);

% Create figure
figure('NumberTitle', 'off', 'Name','All ACROSS FVU: Breathing');

if roi_solve_flag
    % ROIs
    subplot(2,1,1)
    bar(fvu_roi_breathing_adj_all_gen, 'FaceColor', col_to_use)
    title('FVUs for Breathing-Adjusted ROIs (Across All Trials)')
    xlabel('ROI')
    ylabel('FVU')
    set(gca,'FontSize', 16)
    box off
    xlim([0.5 size(fil_roi_data{1},1)+0.5])
end

% Glomeruli
subplot(2,1,2)
bar(fvu_glo_breathing_adj_all_gen, 'FaceColor', col_to_use)
title('FVUs for Breathing-Adjusted Glomeruli (Across All Trials)')
xlabel('Glomerulus')
ylabel('FVU')
set(gca,'FontSize', 16)
box off
xlim([0.5 size(fil_glo_data{1},1)+0.5])

%% ALL ACROSS TRIALS: PLOT FVUs ACROSS TRIALS FOR STIMULUS

% Colors
color_order = get(gca, 'ColorOrder');
col_to_use = color_order(1,:);

% Create figure
figure('NumberTitle', 'off', 'Name','All ACROSS FVU: Stimulus');

if roi_solve_flag
    % ROIs
    subplot(2,1,1)
    bar(fvu_roi_stimulus_adj_all_gen, 'FaceColor', col_to_use)
    title('FVUs for Stimulus-Adjusted ROIs (Across All Trials)')
    xlabel('ROI')
    ylabel('FVU')
    set(gca,'FontSize', 16)
    box off
    xlim([0.5 size(fil_roi_data{1},1)+0.5])
end

% Glomeruli
subplot(2,1,2)
bar(fvu_glo_stimulus_adj_all_gen, 'FaceColor', col_to_use)
title('FVUs for Stimulus-Adjusted Glomeruli (Across All Trials)')
xlabel('Glomerulus')
ylabel('FVU')
set(gca,'FontSize', 16)
box off
xlim([0.5 size(fil_glo_data{1},1)+0.5])

%% PLOT FIL, FIL BREATHING/STIMULUS/BOTH, & ADJUSTED

line_width = 3;

if plot_roi_flag==1
    for tt=1:num_trials
        % ROIs
        % Breathing
        for ii=1:size(fil_roi_data{tt},1)
            figure('NumberTitle', 'off', 'Name', sprintf('Trial %i - R%i',...
                   trial_select(tt),ii));
            if ~contains(version, '2016')
                h_sup = suptitle(sprintf('Trial %i - ROI %i',trial_select(tt),ii));
            end

            ax1 = subplot(3,1,1);
            plot(time{tt}(b_starting_index:end), fil_roi_data{tt}(ii,b_starting_index:end),...
                 'DisplayName', sprintf('Filtered ROI %i',ii), 'LineWidth', line_width)
            hold on
            plot(time{tt}(b_starting_index:end), roi_filtered_breathing_quad{tt}(ii,:),...
                 'DisplayName', 'Filtered Breathing', 'LineWidth', line_width)
            hold off
            legend show
            xlim(time_range_to_plot)
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            grid on

            ax2 = subplot(3,1,2);
            plot(time{tt}(b_starting_index:end), fil_roi_data{tt}(ii,b_starting_index:end),...
                 'DisplayName', sprintf('Filtered ROI %i',ii), 'LineWidth', line_width)
            hold on
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(3,:);
            plot(time{tt}(b_starting_index:end), roi_breathing_adjusted{tt}(ii,:),...
                 'DisplayName', sprintf('Adjusted ROI %i',ii),...
                 'Color', col_to_use, 'LineWidth', line_width)
            hold off
            xlim(time_range_to_plot)
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            legend show
            title(sprintf('FVU: %0.2f', fvu_roi_breathing_adj{tt}(ii)))
            grid on

            ax3 = subplot(3,1,3);
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(3,:);
            plot(time{tt}(b_starting_index:end), roi_breathing_adjusted{tt}(ii,:),...
                 'DisplayName', sprintf('Adjusted ROI %i',ii),...
                 'Color', col_to_use, 'LineWidth', line_width)
            xlim(time_range_to_plot)
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            legend
            title(sprintf('FVU: %0.2f', fvu_roi_breathing_adj{tt}(ii)))
            grid on

            linkaxes([ax1 ax2 ax3],'xy')
        end
    end
end

% Glomeruli
% Breathing
if plot_glo_flag==1
    for tt=1:num_trials
        for ii=1:size(fil_glo_data{tt},1)
            figure('NumberTitle', 'off', 'Name', sprintf('Trial %i - G%i',...
                   trial_select(tt),ii));
            if ~contains(version, '2016')
                h_sup = suptitle(sprintf('Glomerulus %i',ii));
            end

            ax1 = subplot(3,1,1);
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(9,:);
            plot(time{tt}(b_starting_index:end), zeros(size(time{tt}(b_starting_index:end))),...
                 '--','Color',col_to_use,'LineWidth',2,'HandleVisibility','off')
            hold on
            set(gca,'ColorOrderIndex',1)
            plot(time{tt}(b_starting_index:end), fil_glo_data{tt}(ii,b_starting_index:end),...
                 'DisplayName', sprintf('Filtered Glomerulus %i',ii), ...
                 'LineWidth', line_width)
            plot(time{tt}(b_starting_index:end), glo_filtered_breathing_quad{tt}(ii,:),...
                 'DisplayName', 'Filtered Breathing', 'LineWidth', line_width)
            hold off
            legend show
            xlim(time_range_to_plot)
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            ylabel('\Delta F/F')
            %grid on

            text(mean(xlim),max(ylim)+(0.35*max(ylim)),...
                 sprintf('Glomerulus %i',ii),...
                 'FontSize', 30, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

            ax2 = subplot(3,1,2);
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(9,:);
            plot(time{tt}(b_starting_index:end), zeros(size(time{tt}(b_starting_index:end))),...
                 '--','Color',col_to_use,'LineWidth',2,'HandleVisibility','off')
            hold on
            set(gca,'ColorOrderIndex',1)
            plot(time{tt}(b_starting_index:end), fil_glo_data{tt}(ii,b_starting_index:end),...
                 'DisplayName', sprintf('Filtered Glomerulus %i',ii), ...
                 'LineWidth', line_width)
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(3,:);
            plot(time{tt}(b_starting_index:end), glo_breathing_adjusted{tt}(ii,:),...
                 'DisplayName', sprintf('Adjusted Glomerulus %i',ii),...
                 'Color', col_to_use, 'LineWidth', line_width)
            hold off
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            xlim(time_range_to_plot)
            title(sprintf('FVU: %0.2f', fvu_glo_breathing_adj{tt}(ii)))
            legend show
            ylabel('\Delta F/F')
            %grid on

            ax3 = subplot(3,1,3);
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(9,:);
            plot(time{tt}(b_starting_index:end), zeros(size(time{tt}(b_starting_index:end))),...
                 '--','Color',col_to_use,'LineWidth',2)
            hold on
            set(gca,'ColorOrderIndex',1)
            color_order = get(gca, 'ColorOrder');
            col_to_use = color_order(3,:);
            plot(time{tt}(b_starting_index:end), glo_breathing_adjusted{tt}(ii,:),...
                 'DisplayName', sprintf('Adjusted Glomerulus %i',ii),...
                 'Color', col_to_use, 'LineWidth', line_width)
            hold off
            xlabel('Time (s)')
            set(gca,'FontSize', 16)
            xlim(time_range_to_plot)
            title(sprintf('FVU: %0.2f', fvu_glo_breathing_adj{tt}(ii)))
            legend
            ylabel('\Delta F/F')
            %grid on

            linkaxes([ax1 ax2 ax3],'xy')
        end
    end
end



%% Check out correlations between ROIs/Glomeruli after breathing/stim removal

if run_individual_fit % only run if flag set
    for tt=1:num_trials
        figure(corr_fig_hdl{tt})

        % FIL ROIs with breathing adjusted
        RHO_roi_breathing_adjusted{tt} = corr(roi_breathing_adjusted{tt}');

        %figure
        subplot(3,2,5)
        imagesc(RHO_roi_breathing_adjusted{tt})
        title('RHO for filtered roi data with breathing adjusted')
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('roi')
        ylabel('roi')
        if isempty(c_axis)
            c_limits_roi{tt}(3,:) = caxis;
            c_max = max(max(RHO_roi_breathing_adjusted{tt}-...
                    diag(diag(RHO_roi_breathing_adjusted{tt}))));
            c_limits_roi{tt}(3,2) = c_max;
        else
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image


        % FIL GLOMERULI with breathing adjusted
        RHO_glo_breathing_adjusted{tt} = corr(glo_breathing_adjusted{tt}');

        %figure
        subplot(3,2,6)
        imagesc(RHO_glo_breathing_adjusted{tt})
        title('RHO for filtered glomeruli data with breathing adjusted')
        h_bar = colorbar;
        ylabel(h_bar,'RHO')
        set(gca, 'FontSize', 24)
        xlabel('glomerulus')
        ylabel('glomerulus')
        if isempty(c_axis)
            c_limits_glo{tt}(3,:) = caxis;
            c_max = max(max(RHO_glo_breathing_adjusted{tt}-...
                    diag(diag(RHO_glo_breathing_adjusted{tt}))));
            c_limits_glo{tt}(3,2) = c_max;
        else
            caxis(c_axis)
        end
        if isempty(colorbar_m)
            colormap(redblue())
        else
            colormap(redblue(colorbar_m))
        end
        axis image

        % Set caxis across all plots if caxis empty
        for ii =1:6
            subplot(3,2,ii)
            if rem(ii, 2) == 0 % even->glo
                lim_to_use = max([max(c_limits_glo{tt}(:)) abs(min(c_limits_glo{tt}(:)))]);
                caxis([-lim_to_use lim_to_use])
            else % odd-> roi
                lim_to_use = max([max(c_limits_roi{tt}(:)) abs(min(c_limits_roi{tt}(:)))]);
                caxis([-lim_to_use lim_to_use])
            end
        end
    end
end


%% ALL: Check out correlations between ROIs/Glomeruli after breathing removal

figure(all_corr_fig_hdl)

if roi_solve_flag
    % FIL ROIs with breathing adjusted
    RHO_all_roi_breathing_adjusted = corr(roi_breathing_adjusted_all_cat');

    subplot(3,2,5)
    imagesc(RHO_all_roi_breathing_adjusted)
    title('RHO for all filtered roi data with breathing adjusted')
    h_bar = colorbar;
    ylabel(h_bar,'RHO')
    set(gca, 'FontSize', 24)
    xlabel('roi')
    ylabel('roi')
    if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
        c_limits_roi_all_trial(3,:) = caxis;
        c_max = max(max(RHO_all_roi_breathing_adjusted-...
                diag(diag(RHO_all_roi_breathing_adjusted))));
        c_limits_roi_all_trial(3,2) = c_max;
    else % Use c_axis if given, e.g. c_axis = [-1 1]
        caxis(c_axis)
    end
    if isempty(colorbar_m)
        colormap(redblue())
    else
        colormap(redblue(colorbar_m))
    end
    axis image
end

% FIL Glomeruli with breathing adjusted
RHO_all_glo_breathing_adjusted = corr(glo_breathing_adjusted_all_cat');

subplot(3,2,6)
imagesc(RHO_all_glo_breathing_adjusted)
title('RHO for all filtered glomeruli data with breathing adjusted')
h_bar = colorbar;
ylabel(h_bar,'RHO')
set(gca, 'FontSize', 24)
xlabel('glomerulus')
ylabel('glomerulus')
if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
    c_limits_glo_all_trial(3,:) = caxis;
    c_max = max(max(RHO_all_glo_breathing_adjusted-...
            diag(diag(RHO_all_glo_breathing_adjusted))));
    c_limits_glo_all_trial(3,2) = c_max;
else % Use c_axis if given, e.g. c_axis = [-1 1]
    caxis(c_axis)
end
if isempty(colorbar_m)
    colormap(redblue())
else
    colormap(redblue(colorbar_m))
end
axis image

% Set caxis across all plots if caxis empty
for ii =1:6
    if (ignore_raw_for_c_normalization==1)&&(ii<3)

    else
        subplot(3,2,ii)
        if rem(ii, 2) == 0 % even->glo
            if ignore_raw_for_c_normalization==1
                c_limits_glo_all_trial(1,:) = 0;
            end
            lim_to_use = max([max(c_limits_glo_all_trial(:)) ...
                              abs(min(c_limits_glo_all_trial(:)))]);
            caxis([-lim_to_use lim_to_use])
        else % odd-> roi
            if ignore_raw_for_c_normalization==1
                c_limits_roi_all_trial(1,:) = 0;
            end
            lim_to_use = max([max(c_limits_roi_all_trial(:)) ...
                              abs(min(c_limits_roi_all_trial(:)))]);
            caxis([-lim_to_use lim_to_use])
        end
    end
end

%% ALL: Check out correlations between ROIs/Glomeruli after breathing+stim
%  removal
figure(all_corr_fig_br_stim_hdl)

% FIL Glomeruli with breathing adjusted
RHO_all_glo_breathing_adjusted = corr(glo_breathing_adjusted_all_cat');

subplot(2,2,2)
imagesc(RHO_all_glo_breathing_adjusted)
title('Adjusted for Breathing')
h_bar = colorbar;
ylabel(h_bar,'\rho')
set(gca, 'FontSize', 24)
xlabel('Glomerulus')
ylabel('Glomerulus')
if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
    c_limits_glo_all_trial_br_stim(2,:) = caxis;
    c_max = max(max(RHO_all_glo_breathing_adjusted-...
            diag(diag(RHO_all_glo_breathing_adjusted))));
    c_limits_glo_all_trial_br_stim(2,2) = c_max;
else % Use c_axis if given, e.g. c_axis = [-1 1]
    caxis(c_axis)
end
if isempty(colorbar_m)
    colormap(redblue())
else
    colormap(redblue(colorbar_m))
end
axis image



% FIL Glomeruli with stimulus adjusted
RHO_all_glo_stimulus_adjusted = corr(glo_stimulus_adjusted_all_cat');

subplot(2,2,3)
imagesc(RHO_all_glo_stimulus_adjusted)
title('Adjusted for Stimulus')
h_bar = colorbar;
ylabel(h_bar,'\rho')
set(gca, 'FontSize', 24)
xlabel('Glomerulus')
ylabel('Glomerulus')
if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
    c_limits_glo_all_trial_br_stim(3,:) = caxis;
    c_max = max(max(RHO_all_glo_stimulus_adjusted-...
            diag(diag(RHO_all_glo_stimulus_adjusted))));
    c_limits_glo_all_trial_br_stim(3,2) = c_max;
else % Use c_axis if given, e.g. c_axis = [-1 1]
    caxis(c_axis)
end
if isempty(colorbar_m)
    colormap(redblue())
else
    colormap(redblue(colorbar_m))
end
axis image




% FIL Glomeruli with breathing+stim adjusted
RHO_all_glo_breathing_stimulus_adjusted = corr(glo_both_adjusted_all_cat');

subplot(2,2,4)
imagesc(RHO_all_glo_breathing_stimulus_adjusted)
title('Adjusted for Br + Stim')
h_bar = colorbar;
ylabel(h_bar,'\rho')
set(gca, 'FontSize', 24)
xlabel('Glomerulus')
ylabel('Glomerulus')
if isempty(c_axis) % If c_axis unspecified/empty, normalize c_axis for trial
    c_limits_glo_all_trial_br_stim(4,:) = caxis;
    c_max = max(max(RHO_all_glo_breathing_stimulus_adjusted-...
            diag(diag(RHO_all_glo_breathing_stimulus_adjusted))));
    c_limits_glo_all_trial_br_stim(4,2) = c_max;
else % Use c_axis if given, e.g. c_axis = [-1 1]
    caxis(c_axis)
end
if isempty(colorbar_m)
    colormap(redblue())
else
    colormap(redblue(colorbar_m))
end
axis image




% Set caxis across all plots if caxis empty
for ii =1:4
    subplot(2,2,ii)
    lim_to_use = max([max(c_limits_glo_all_trial_br_stim(:)) ...
                      abs(min(c_limits_glo_all_trial_br_stim(:)))]);
    caxis([-lim_to_use lim_to_use])
end

%% PRODUCE CONCATENATED RAW & DFF WITH B START INDEXING
for tt=1:num_trials
        raw_glo_data_b_start{tt} = raw_glo_data{tt}(:,b_starting_index:end);
        dff_glo_data_b_start{tt} = dff_glo_data{tt}(:,b_starting_index:end);
end
raw_glo_data_cat = cell2mat(raw_glo_data_b_start);
dff_glo_data_cat = cell2mat(dff_glo_data_b_start);

%% PRODUCE CONCATENATED FILTERED STIM WITH B START INDEXING
glo_filtered_stimulus_quad_all_b_start_cat = cell2mat(glo_filtered_stimulus_quad_all);

%% SHOW GLOMERULAR TIME SERIES
h(1) = figure('NumberTitle', 'off', 'Name', 'Raw');
plot(time_generic(1:length(raw_glo_data_cat)),raw_glo_data_cat')
title('Raw')
set(gca,'FontSize',20)
ylim_save(1,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(2) = figure('NumberTitle', 'off', 'Name', 'DFF');
plot(time_generic(1:length(raw_glo_data_cat)),dff_glo_data_cat')
title('DFF')
set(gca,'FontSize',20)
ylim_save(2,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(3) = figure('NumberTitle', 'off', 'Name', 'Breathing Adjusted');
plot(time_generic(1:length(raw_glo_data_cat)),glo_breathing_adjusted_all_cat')
title('Breathing Adjusted')
set(gca,'FontSize',20)
ylim_save(3,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(4) = figure('NumberTitle', 'off', 'Name', 'Stimulus Adjusted');
plot(time_generic(1:length(raw_glo_data_cat)),glo_stimulus_adjusted_all_cat')
title('Stimulus Adjusted')
set(gca,'FontSize',20)
ylim_save(4,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(5) = figure('NumberTitle', 'off', 'Name', 'Both Adjusted');
plot(time_generic(1:length(raw_glo_data_cat)),glo_both_adjusted_all_cat')
title('Both Adjusted')
set(gca,'FontSize',20)
ylim_save(5,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(6) = figure('NumberTitle', 'off', 'Name', 'Convolved Stim');
plot(time_generic(1:length(raw_glo_data_cat)),...
        glo_filtered_stimulus_quad_all_b_start_cat', 'LineWidth',3)
title('Convolved Stimulus')
set(gca,'FontSize',20)
ylim_save(6,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

h(7) = figure('NumberTitle', 'off', 'Name', 'FIL');
plot(time_generic(1:length(raw_glo_data_cat)),fil_glo_data_b_start_cat')
title('Filtered Glomeruli')
set(gca,'FontSize',20)
ylim_save(7,:) = ylim;
xlabel('Time (s)')
box off
xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])

for ii=1:7
    figure(h(ii))
    ylim([min(ylim_save(:)) max(ylim_save(:))])
end

%% SHOW ABOVE PLOTS BUT WITH SEPARATE LINES FOR EACH GLOMERULUS
data_cell =         {raw_glo_data_cat;...
                     dff_glo_data_cat;...
                     fil_glo_data_b_start_cat;...
                     glo_breathing_adjusted_all_cat;...
                     glo_stimulus_adjusted_all_cat;...
                     glo_both_adjusted_all_cat;...
                     glo_filtered_stimulus_quad_all_b_start_cat};

data_cell_max =    max(cell2mat(cellfun(@(x) max(x'),data_cell,'UniformOutput',0)));


data_cell_titles = {'Raws';...
                    'DFFs';...
                    'FILs';...
                    'Breathing Adjusteds';...
                    'Stimulus Adjusteds';...
                    'Both Adjusteds';...
                    'Convolved Stims'};

for ii=1:length(data_cell)

    % Plot
    h(ii) = figure('NumberTitle', 'off', 'Name', data_cell_titles{ii});
    for jj=1:size(raw_glo_data_cat,1)
        if ii==7 % Linewidth=2
            plot(time_generic(1:length(raw_glo_data_cat)),...
             (data_cell{ii}(jj,:)/data_cell_max(jj))+jj,...
             'LineWidth',2)
        else     % Linewidth=1
            plot(time_generic(1:length(raw_glo_data_cat)),...
             (data_cell{ii}(jj,:)/data_cell_max(jj))+jj)
        end
        hold on
    end
    hold off

    % Plot settings
    title(data_cell_titles{ii})
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    box off
    xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])
    ylabel('Glomerulus (Individually Normalized to Max Across Processing')
    ylim([1 size(raw_glo_data_cat,1)+1])
    yticks(1:size(raw_glo_data_cat,1))
    ylim([min(ylim)-0.5 max(ylim)])
end


%% ANALYSIS OF RELATIONSHIPS BETWEEN glo_stimulus_adjusted_all_cat,
%   adjusted for part that stimulus is removed effectively

time_index_end = 5000;

chopped_indices = [];
for ii=1:length(trial_select)
    chopped_indices = [chopped_indices ((1:time_index_end)+...
        (size(data_cell{1},2)/length(trial_select)*(ii-1)))];
end

% Create chopped figures
for ii=1:length(data_cell)

    % Plot
    h(ii) = figure('NumberTitle', 'off', 'Name', ['Chopped: ' data_cell_titles{ii}]);
    for jj=1:size(raw_glo_data_cat,1)
        plot(time_generic(1:(time_index_end*length(trial_select))),...
             (data_cell{ii}(jj,chopped_indices)/data_cell_max(jj))+jj)
        hold on
    end
    hold off

    for jj=1:num_trials
        hold on
        fill(time_generic([1 1 time_index_end time_index_end]+((jj-1)*time_index_end)),...
             [11 11.25 11.25 11],color_order(jj,:))
        text(time_generic(100+((jj-1)*time_index_end)),11.1,...
            sprintf('Trial %i',trial_select(jj)))
    end
    hold off

    % Plot settings
    if exist('odorant_name','var')
        title([odorant_name ' - ' data_cell_titles{ii}])
    else
        title(data_cell_titles{ii})
    end
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    box off
    xlim([time_generic(1) time_generic(length(raw_glo_data_cat))])
    %ylabel('Glomerulus (Individually Normalized to Max Across Processing')
    ylabel('Glomerulus')
    ylim([1 size(raw_glo_data_cat,1)+1])
    yticks(1:size(raw_glo_data_cat,1))
    ylim([min(ylim)-0.5 max(ylim)+.25])

    % Format
    set(h(ii), 'Position', [1269         390         886         497])
    h(ii).Renderer='Painters';
    set(h(ii),'Units','Inches');
    pos = get(h(ii),'Position');
    set(h(ii),'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)])

    % Save fig
    saveas(h(ii), ['saved_figures/' get(get(gca,'Title'),'String') '.pdf'])
    savefig(h(ii), ['saved_figures/' get(get(gca,'Title'),'String') '.fig'])
end

%% Compute covariance & save
stim_adj_covered = cell2mat(cellfun(@(x) x(:,1:time_index_end), ...
    glo_stimulus_adjusted_all,'UniformOutput',0));
fil_glo_covered = cell2mat(cellfun(@(x) x(:,1:time_index_end), ...
    fil_glo_data_b_start,'UniformOutput',0));

stim_adj_covered_cov = cov(stim_adj_covered');
stim_adj_covered_corr = corr(stim_adj_covered');
stim_adj_covered_pcorr = partialcorr(stim_adj_covered');
save('var_stim_adj.mat','stim_adj_covered_cov')

fil_glo_covered_cov = cov(fil_glo_covered');
fil_glo_covered_corr = corr(fil_glo_covered');
fil_glo_covered_pcorr = partialcorr(fil_glo_covered');

% Plot pcorr
figure('NumberTitle', 'off', 'Name', 'Before Stimulus Adjusted: PCorr');
imagesc(fil_glo_covered_pcorr)
colorbar
caxis([min(caxis) max(max(fil_glo_covered_pcorr(~diag(ones(size(fil_glo_covered,1),1)))))])
colormap bluewhitered
set(gca,'FontSize',20)
title('Partial Correlations Between Glomeruli Before Adjustment')
xlabel('Glomerulus')
ylabel('Glomerulus')

figure('NumberTitle', 'off', 'Name', 'Stimulus Adjusted: PCorr');
imagesc(stim_adj_covered_pcorr)
colorbar
caxis([min(caxis) max(max(stim_adj_covered_pcorr(~diag(ones(size(stim_adj_covered,1),1)))))])
colormap bluewhitered
set(gca,'FontSize',20)
title('Partial Correlations Between Stimulus-Adjusted Glomeruli')
xlabel('Glomerulus')
ylabel('Glomerulus')

%% Plot graph of corr & pcorr
h_graph = figure('NumberTitle', 'off','Name','Graphs');
subplot(1,2,1)
h_bar = glom_plot(stim_adj_covered_corr);

title('Corr')
set(h_bar, 'NumberTitle', 'off', 'Name', 'Bar: Corr')

figure(h_graph)
subplot(1,2,2)
h_bar = glom_plot(stim_adj_covered_pcorr);

title('Pcorr')
set(h_bar, 'NumberTitle', 'off', 'Name', 'Bar: Pcorr')

figure(h_graph)
subplot(1,2,1)
title('Corr')
caxis_save(1,:) = caxis;
subplot(1,2,2)
title('Pcorr')
caxis_save(2,:) = caxis;
for ii=1:2
    subplot(1,2,ii)
    caxis([min(caxis_save(:)) max(caxis_save(:))])
end

%% FPCA
% Smoothing parameter for FPCA -> can be vector
alphavs = 400;
[p,N] = size(glo_both_adjusted_all{1}); % p neurons, N time points
Omegu = sparse(toeplitz([2, -1, zeros(1, p - 2)])); % neurons
Omegv = sparse(toeplitz([2, -1, zeros(1, N - 2)])); % time
Omegv = Omegv*Omegv; % Second order differences matrix

[U,V,d,optaus,optavs,Xhat,bicu,bicv] = fpca_nested_bic(glo_both_adjusted_all{1},2,0,...
                            alphavs,Omegu,Omegv,0,0,1000,5);


%% Save & email
save([odorant_name '_all_FINAL.mat'])
if ~ismac
   email_sending('olf_done','done')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
