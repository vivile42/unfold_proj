% unfold for cardiac phase separate
%cd needs to be set in the tsk folder (BBC/WP1/data/EEG/tsk)
%cd
%% preproc
% Load the dataset
G_N = ["g01", "g03", 'g08', 'g10', 'g11', 'g12', 'g15', 'g16', 'g17', 'g19', 'g23', 'g24','g25', 'g28', 'g32', 'g33', 'g34', 'g37', 'g38', 'g39', 'g40', 'g41', 'g42', 'g44', 'g45', 'g46', 'g47', 'g49', 'g51', 'g52']

%length(G_N)
for g=3:length(G_N)

    g_num=G_N(g);
    disp(g_num)
    EEG = pop_loadset(sprintf('./preproc/%s/%s_deconv/%s_n_tsk_eeg_data.set',g_num,g_num,g_num));
    % Load the eog log
    dir_log=sprintf('./preproc/%s/%s_epochs/%s_n_tsk_eog_log.txt',g_num,g_num,g_num);
    blink_cell=readcell(dir_log);
    blink_cell=blink_cell(2:end,1);
    cleaned = erase(blink_cell, '(');
    bad_blink_times = str2double(cleaned);
    % Extract events
    event_table=readtable(sprintf('./preproc/%s/%s_deconv/%s_n_tsk_event_table.csv',g_num,g_num, g_num));
    
    %add awareness to event strct in uf
    awareness_values= event_table.awareness;
    card_values=event_table.card_phase;
    rsp_values=event_table.rsp_phase;
    
    for i=1:length(EEG.event);
        EEG.event(i).awareness= awareness_values{i};
        EEG.event(i).card_phase= card_values{i};
        EEG.event(i).rsp_phase= rsp_values{i};
    end
    
    % Keep only the relevant columns: time and type
    event_table = event_table(:, {'time', 'event','awareness','card_phase','rsp_phase'});
    gh
    % Rename columns to match Unfold's expectations
    event_table.Properties.VariableNames = {'time', 'event','awareness','card_phase','rsp_phase'};
    
    % Display the event table
    %disp(event_table);
    
    %% light cleaning

    % Detect high-amplitude artifacts
    
    cfg = [];
    cfg.threshold = 100;  % Threshold in µV (e.g., ±200 µV)
    cfg.timewindow = [-0.2 0.8];  % Time window around events to check (in seconds)
    cfg.eventtypes = {'hep', 'vep'};  % Event types to consider
    sr=256
    % Run artifact detection
    EEG_exclud = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',200);
    
    
    % --- Convert times to sample indices ---
    bad_blink_tf = round(bad_blink_times * sr);
    
    pre_samples=250;
    post_samples=250;
    
    % --- Define blink exclusion intervals [start, end] ---
    start_tf = bad_blink_tf - pre_samples;
    end_tf   = bad_blink_tf + post_samples;
    bad_blink_intervals = [start_tf, end_tf];
    
    % --- Append to your existing exclusion matrix ---
    eeg_exclude = [EEG_exclud; bad_blink_intervals];
    
    % --- Sort everything by start time ---
    eeg_exclude = sortrows(eeg_exclude, 1);
    % Check the marked artifacts
    exc_seg=sum(EEG_exclud(:,2)-EEG_exclud(:,1));
    tot_seg=EEG.pnts;
    ratio_exc=round(exc_seg/tot_seg*100);
    %display ratio of rejected data
    log_rej=sprintf('percentage of data marked as artefact: %d percent',ratio_exc);
    log_rej_all=sprintf('\n %s percentage of data marked as artefact: %d percent',g_num, ratio_exc);
    disp(log_rej);
    %save log with the ratio
    
    fid=fopen(sprintf('./preproc/%s/%s_deconv/cardfit/%s_n_tsk_deconv_artrej_log_cardfit.txt',g_num,g_num,g_num),'w');
  
    fprintf(fid,log_rej);
    
    
    
    %% setup unfold design
    % Initialize Unfold
    cfg = [];
    
    % Define the time window for the basis function
    cfg.timelimits = [-0.2, 0.8];
    cfg.winrej = eeg_exclude;
    
    % Split event table by cardiac phase
    sys_mask = strcmp(event_table.card_phase, 'sys');
    event_table_sys = event_table(sys_mask, :);
    event_table_dia = event_table(~sys_mask, :);
    
    % Create separate EEG datasets by filtering EEG.event
    % Copy the EEG structure
    EEG_sys = EEG;
    EEG_dia = EEG;
    
    % Filter EEG.event for sys events
    sys_event_mask = strcmp({EEG.event.card_phase}, 'sys');
    EEG_sys.event = EEG.event(sys_event_mask);
    
    % Filter EEG.event for dia events
    dia_event_mask = ~sys_event_mask;
    EEG_dia.event = EEG.event(dia_event_mask);
    
    % Verify the number of events
    
    log_sys=sprintf('\n Number of sys events in EEG_sys: %d',length(EEG_sys.event))
    log_dia=sprintf('\n Number of dia events in EEG_dia: %d',length(EEG_dia.event))
    fprintf(fid,log_sys);
    fprintf(fid,log_dia);
    
    fclose(fid);
    % Initialize Unfold structures for sys and dia
    cfg_sys = cfg;
    cfg_sys.event = event_table_sys;
    cfg_sys.formula = {'y ~ 1', 'y ~ 1'};  % Intercept-only model for each event type
    cfg_sys.eventtypes = {'hep', 'vep'};
    uf_sys = uf_designmat(EEG_sys, cfg_sys);
    
    cfg_dia = cfg;
    cfg_dia.event = event_table_dia;
    cfg_dia.formula = {'y ~ 1', 'y ~ 1'};
    cfg_dia.eventtypes = {'hep', 'vep'};
    uf_dia = uf_designmat(EEG_dia, cfg_dia);
    
    % Plot design matrices to verify (optional)
    mkdir(sprintf('./preproc/%s/%s_deconv/figures',g_num,g_num))
    uf_plotDesignmat(uf_sys, 'sort', 1);
    title('design mat sys');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_design_mat_sys.png',  g_num,g_num,g_num), 'png');
    
    uf_plotDesignmat(uf_dia, 'sort', 1);
    title('design mat dia');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_design_mat_dia.png',  g_num,g_num,g_num), 'png');
    
    %% fit models
    cfgTimeexpand = [];
    cfgTimeexpand.timelimits = [-0.2, 0.8];
    
    % Fit sys model
    uf_sys = uf_timeexpandDesignmat(uf_sys, cfgTimeexpand);
    uf_sys = uf_continuousArtifactExclude(uf_sys);
    uf_sys = uf_glmfit(uf_sys);
    EEG_epoch_sys = uf_epoch(uf_sys, 'timelimits', cfgTimeexpand.timelimits);
    EEG_epoch_sys = uf_glmfit_nodc(EEG_epoch_sys);
    ufresult_sys = uf_condense(EEG_epoch_sys);
    ufresult_sys = uf_predictContinuous(ufresult_sys);
    
    % Fit dia model
    uf_dia = uf_timeexpandDesignmat(uf_dia, cfgTimeexpand);
    uf_dia = uf_continuousArtifactExclude(uf_dia);
    uf_dia = uf_glmfit(uf_dia);
    EEG_epoch_dia = uf_epoch(uf_dia, 'timelimits', cfgTimeexpand.timelimits);
    EEG_epoch_dia = uf_glmfit_nodc(EEG_epoch_dia);
    ufresult_dia = uf_condense(EEG_epoch_dia);
    ufresult_dia = uf_predictContinuous(ufresult_dia);
    
    
    
    % Save results for both conditions
    save_beta_results(ufresult_sys, g_num, 'sys');
    save_beta_results(ufresult_dia, g_num, 'dia');
    
    disp('Export complete');
    
    %% plot result
    % Plot parameters for specific channels with condition in title
    
    figure;
    uf_plotParam(ufresult_sys, 'channel', 'A11', 'plotSeparate', 'event');
    title('sys: Channel A11');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_sys_A11.png', g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParam(ufresult_sys, 'channel', 'B14', 'plotSeparate', 'event');
    title('sys: Channel B14');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_sys_B14.png', g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParam(ufresult_dia, 'channel', 'A11', 'plotSeparate', 'event');
    title('dia: Channel A11');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_dia_A11.png', g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParam(ufresult_dia, 'channel', 'B14', 'plotSeparate', 'event');
    title('dia: Channel B14');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_dia_B14.png',  g_num,g_num,g_num), 'png');
    
    % Plot topographic maps with condition in title
    figure;
    uf_plotParamTopo(ufresult_sys, 'betaSetName', 'beta');
    title('sys: Topographic Beta');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_sys_topo_beta.png',  g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParamTopo(ufresult_sys, 'betaSetName', 'beta_nodc');
    title('sys: Topographic Beta NoDC');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_sys_topo_beta_nodc.png',  g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParamTopo(ufresult_dia, 'betaSetName', 'beta');
    title('dia: Topographic Beta');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_dia_topo_beta.png',  g_num,g_num,g_num), 'png');
    
    figure;
    uf_plotParamTopo(ufresult_dia, 'betaSetName', 'beta_nodc');
    title('dia: Topographic Beta NoDC');
    saveas(gcf, sprintf('./preproc/%s/%s_deconv/figures/%s_plot_dia_topo_beta_nodc.png', g_num,g_num,g_num), 'png');
    
    close all
end

%% save Beta to CSV
% Function to save beta results for a given condition
function save_beta_results(ufresult, g_num, condition)
    betas = ufresult.beta;              % [num_predictors x num_timepoints]
    times = ufresult.times;             % 1 x num_timepoints
    labels = ufresult.unfold.colnames;  % one label per beta
    eventtypes = ufresult.unfold.eventtypes;  % one per eventtype
    chan_labels = {ufresult.chanlocs.labels};

    % Debug: Inspect labels and eventtypes
    disp('Debug: labels content:');
    disp(labels);
    disp('Debug: eventtypes content:');
    disp(eventtypes);

    n_evtypes = length(eventtypes);
    n_betas_per_type = length(labels) / n_evtypes;
    final_labels = cell(1, length(labels));

    for i = 1:n_evtypes
        for j = 1:n_betas_per_type
            idx = (i-1)*n_betas_per_type + j;
            % Convert labels{idx} to string, handling cell or non-string cases
            label_str = string(labels{idx});
            if iscell(label_str)
                label_str = label_str{1}; % Extract first element if still a cell
                
            end
            disp(label_str)
            disp(eventtypes{i})
            final_labels{idx} = sprintf('%s_%s_%s', condition, string(eventtypes{i}), label_str);
        end
    end

    %ufresult.deconv = 1;
    table = uf_unfold2csv(ufresult,'deconv',1);
    mkdir(sprintf('./preproc/%s/%s_deconv/cardfit/beta',g_num,g_num))
    dir_log=sprintf('./preproc/%s/%s_deconv/cardfit/beta/%s_beta_dc_%s.csv',g_num,g_num,g_num,condition);
    writetable(table, dir_log);
    
    %ufresult.deconv = 0;
    table = uf_unfold2csv(ufresult,'deconv',0);
    mkdir(sprintf('./preproc/%s/%s_deconv/cardfit/beta',g_num,g_num))
    dir_log=sprintf('./preproc/%s/%s_deconv/cardfit/beta/%s_beta_nodc_%s.csv',g_num,g_num,g_num,condition);
    writetable(table, dir_log);
   
    % Create output folder
    output_dir ='./ana/deconvolution/ep_betas/cardfit/dc';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    % Export one .ep file per predictor
    for p = 1:length(final_labels)
        beta_ep = squeeze(ufresult.beta(:,:,p));  % [time x chan]
        beta_ep = beta_ep';                       % → transpose to [chan x time]
        % Make filename safe
        fname = matlab.lang.makeValidName(final_labels{p})
        fname = fullfile(output_dir, sprintf('%s_%s_dc.ep',g_num, fname));
        % Save as ASCII .ep (space-delimited, like BrainVision)
        dlmwrite(fname, beta_ep, 'delimiter', '\t', 'precision', 6);
    end
    % Create output folder
    output_dir ='./ana/deconvolution/ep_betas/cardfit/nodc';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

        % Export one .ep file per predictor
    for p = 1:length(final_labels)
        beta_ep = squeeze(ufresult.beta_nodc(:,:,p));  % [time x chan]
        beta_ep = beta_ep';                       % → transpose to [chan x time]
        % Make filename safe
        fname = matlab.lang.makeValidName(final_labels{p})
        fname = fullfile(output_dir, sprintf('%s_%s_nodc.ep',g_num, fname));
        % Save as ASCII .ep (space-delimited, like BrainVision)
        dlmwrite(fname, beta_ep, 'delimiter', '\t', 'precision', 6);
    end
end