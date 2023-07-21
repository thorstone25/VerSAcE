% classdef VSXSystem < matlab.mixin.Copyable
%     properties
%         LEsysExpected (1,1) double = 0
%         AcqSlotsExpected (1,:) double = [1 1 1 1 1 1 1 1]
%         UTAexpected double
%         Product
%         FrequencyExpected (1,1) string = "SF"
%         activeClampExpected (1,1) double = 0
%         Frequency (1,1) string
%         UTA (1,1) string
%     end
%     methods
%         function obj = VSXSystem()
%             if isempty(obj.Frequency)
%                 obj.Frequency = string(VSXTrans.frequency);
%             end
%         end
%     end
% end
% 
