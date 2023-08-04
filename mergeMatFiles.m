function mergeMatFiles(fileNames, saveLocation)
    
    mergedData = struct();
    
    
    for i = 1:numel(fileNames)
        data = load(fileNames{i});
        if isempty(fieldnames(mergedData))
            mergedData = data;
        else
            f = fieldnames(data);
            for j = 1:length(f)
                mergedData.(f{j}) = data.(f{j});
            end
        end
    end
    
    save(fullfile(saveLocation, 'mergedData.mat'), 'mergedData');
    disp('Merged data saved to mergedData.mat');
    
end