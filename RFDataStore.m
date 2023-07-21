function RFDataStore(Resource)
    global toggle;
    if toggle
        disp('RF')
        RcvData{1} = Resource.RcvBuffer;%rcvbuf;
        full_path = 'Y:/coop_03_Stanford/verasonics/Ameya'; % where to save?
        disp('Saving RF data...')
        save('-v7.3', '-nocompression', [full_path, '/RF_DATA_L7_4' datestr(now,'yyyymmdd_HHMMSS') '.mat'],...
         'RcvData'); % maybe save as .dat file?
        disp('RF data saved!')

        % clear rcvbuf
    end
    return
end