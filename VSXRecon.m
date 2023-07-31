classdef VSXRecon < matlab.mixin.Copyable
    properties
        senscutoff (1,1) double = 0.6
        pdatanum (1,1) VSXPData
        rcvBufFrame (1,1) double {mustBeInteger} = -1 % -1 to process last frame acquired
        IntBufDest VSXInterBuffer {mustBeScalarOrEmpty} = VSXInterBuffer.empty
        IntBufDestFrm (1,1) double {mustBeInteger} = -1 % -1 to process last frame acquired
        ImgBufDest VSXImageBuffer {mustBeScalarOrEmpty} = VSXImageBuffer.empty
        ImgBufDestFrm (1,1) double {mustBeInteger} = -1 % -1 to process last frame acquired
        RINums (1,:) VSXReconInfo = VSXReconInfo.empty
        newFrameTimeout (1,1) double = 1000 % timeout in ms
        rcvLUT (:,:) uint16 = []
    end
    methods
        function obj = VSXRecon(kwargs)
            arguments
                kwargs.?VSXRecon
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end
