function RFDataStore(RData, varargin)
global TOGGLE_RFDataStore;
global VSXOOD_SAVE_DIR;

if TOGGLE_RFDataStore
    save_dir = VSXOOD_SAVE_DIR;
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non exists
    
    % get save filename
    fileName = string(datetime('now'),'yyMMdd_HHmmss_SSSS') + ".mat";
    fnm = fullfile(save_dir, fileName); % filename
    
    % save raw data
    save(fnm, '-v7.3', '-nocompression', 'RData'); % TODO: maybe save as .dat file?
    
    disp("RF data saved to " + fnm);
    TOGGLE_RFDataStore = false;
end
end
