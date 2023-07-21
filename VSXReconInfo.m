classdef VSXReconInfo < matlab.mixin.Copyable
    properties
        mode (1,:) char 
        txnum (1,1) double 
        rcvnum (1,1) double 
        regionnum (1,1)
        
    end
    methods
        function obj = VSXReconInfo(kwargs)
            arguments
                kwargs.?VSXReconInfo
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end