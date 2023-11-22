file1 = './MatFiles/L12-3vFlashAngles';
file2 = './MatFiles/L12-3vFlashAngles_test';
file3 = './MatFiles/L12-3vFlashAngles_demo';
file4 = 'Matfiles/qups_config.mat';
file5 = 'MatFiles/L11-5vFlashAngles';
% file3 = './MatFiles/L12-3v_192RyLns_Hadamard.mat';

og = load(file5);
vsx = load(file4);
vsx_demo = load(file2);
% hadamard = load(file3);



% % if isequaln(og, vsx)
% %     disp('exact same values!');
% % end

og_vars = fieldnames(og);
vsx_vars = fieldnames(vsx);
% % 
% % % for i = 1:numel(vsx.UI(2).Callback)
% % %     if ~isequal(vsx.UI(2).Callback(i), vsx.UI(2).Callback(i))
% % %         disp(i);
% % %     end
% % % end
% % 
% % 
for i = 1: numel(vsx_vars)
    var = vsx_vars{i};
    if ~isequal(og.(var), vsx.(var))
        fprintf('variable "%s" is different.\n', var);
    end
end

% % 
% % % if isequal(vsx, vsx)
% % %     disp('they are same');
% % % else
% % %     disp('they are still different');
% % % end
% % 
% % 
