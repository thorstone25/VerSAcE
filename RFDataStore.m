function RFDataStore(RData)
    global toggle;
    save_dir  = fullfile(pwd, 'RF_DATA_L7_4'); % make a path relative to the current location
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end % make a directory if non existss
    if toggle
        % full_path = 'Y:/coop_03_Stanford/verasonics/Ameya'; % where to save?
        disp('Saving RF data...')
        save(fullfile(save_dir, [datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
            '-v7.3', '-nocompression', 'RData'); % TODO: maybe save as .dat file?
        disp('RF data saved!');
        toggle = false; % unset variable saving
    end
    return
end