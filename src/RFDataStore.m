function RFDataStore(RData, varargin)
global TOGGLE_RFDataStore;
global VERSACE_PARAMS;

if TOGGLE_RFDataStore
    save_dir = VERSACE_PARAMS.save_dir;
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non exists
    
    % get save filename
    fileName = string(datetime('now'),'yyMMdd_HHmmss_SSSS') + ".mat";
    fnm = fullfile(save_dir, fileName); % filename

    % get current settings
    vs = update_vstruct();
    vs.RData = RData;
    vs.VERSACE_PARAMS = VERSACE_PARAMS;
    
    % save raw data
    save(fnm, '-v7.3', '-nocompression', '-struct', 'vs'); % TODO: maybe save as .dat file?
    
    disp("RF data saved to " + fnm);
    TOGGLE_RFDataStore = false;
end
end
