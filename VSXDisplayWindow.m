classdef VSXDisplayWindow < matlab.mixin.Copyable
    properties
        Type (1,:) char = 'Verasonics' %{mustBeMember(Type,["Verasonics", "Matlab"])} = "Verasonics"
        Title (1,:) char = 'Window title'
% %         mode (1,:) char = '2d'
% %         Orientation (1,:) char = 'xz' %{mustBeMember(Orientation,["xz","yz","xy"])}="xz"
        AxesUnits (1,:) char = 'wavelengths'  %{mustBeMember(AxesUnits,["wavelengths","mm"])}="wavelengths"
        Position (1,:) double = [1 1 512 512]
        ReferencePt (1,3) double
        pdelta (1,1) double = 0.25
        Colormap (256,3) double = "grey"
        numFrames (1,1) double = 1
% %         firstFrame (1,1) double 
% %         lastFrame (1,1) double
% %         clrWindow (1,1) double
    end
    methods
        function obj = VSXDisplayWindow(kwargs)
            arguments
                kwargs.?VSXDisplayWindow
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
end