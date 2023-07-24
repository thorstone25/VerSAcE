function RFDataStore(RData)
    global TOGGLE_RFDataStore;
    save_dir  = fullfile(pwd, "RF_DATA_L7_4"); % make a path relative to the current location
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non existss
    if TOGGLE_RFDataStore
        disp('Saving RF data ...')
        fnm = fullfile(save_dir, [datestr(now,'yyyymmdd_HHMMSS') '.mat']); % filename
        save(fnm, '-v7.3', '-nocompression', 'RData'); % TODO: maybe save as .dat file?
        disp("RF data saved to " + fnm);
        toggle = false; % unset variable saving
    end
    return
end
