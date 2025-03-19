function RFDataStore(RData, varargin)
global TOGGLE_RFDataStore;
global VERSACE_PARAMS;

if TOGGLE_RFDataStore
    % Defaults
    if ~isfield(VERSACE_PARAMS, 'save_dir'), VERSACE_PARAMS.save_dir = '.'; end
    if ~isfield(VERSACE_PARAMS, 'verbose'),  VERSACE_PARAMS.verbose = false; end
    

    save_dir = VERSACE_PARAMS.save_dir;
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non exists
    
    % get save filename
    fileName = string(datetime('now'),'yyMMdd_HHmmss_SSSS') + ".mat";
    fnm = fullfile(save_dir, fileName); % filename

    % get current settings
    vs = update_vstruct();
    if isfield(vs, 'UI'), vs = rmfield(vs, 'UI'); end % exclude UI
    [vs.RData, vs.VERSACE_PARAMS] = deal(RData, VERSACE_PARAMS); % append
    
    % save raw data
    save(fnm, '-v7.3', '-nocompression', '-struct', 'vs'); % TODO: maybe save as .dat file?
    
    if ~isfield(VERSACE_PARAMS, 'verbose'), VERSACE_PARAMS.verbose = false; end
    if VERSACE_PARAMS.verbose, disp("RF data saved to " + fnm); end
    TOGGLE_RFDataStore = false;
end
end
