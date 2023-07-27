function RFDataStore(RData)
    global TOGGLE_RFDataStore;
    save_dir  = fullfile(pwd, "RF_DATA_L7_4"); % make a path relative to the current location
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non exists
    if TOGGLE_RFDataStore
        % get save filename
        disp('Saving RF data ...')
        fileName = datestr(now,'yyyymmdd_HHMMSS') + ".mat";
        fnm = fullfile(save_dir, fileName); % filename
        
        % save raw data
        save(fnm, '-v7.3', '-nocompression', 'RData'); % TODO: maybe save as .dat file?

        % append configuration data
        vsx = load("MatFiles/qups-vsx.mat");
        cnf = load("MatFiles/qups-conf.mat");
        for f = string(fieldnames(cnf))', vsx.(f) = cnf.(f); end
        save(fnm, '-append', '-struct', 'vsx');
        
        disp("RF data saved to " + fnm);
        TOGGLE_RFDataStore = false; % unset variable saving
    end
    return
end
