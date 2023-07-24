classdef VSXDisplayWindow < matlab.mixin.Copyable
    properties
        Type (1,1) string {mustBeMember(Type,["Verasonics", "Matlab"])} = "Verasonics"
        Title (1,1) string = 'Window title'
% %         mode (1,:) char = '2d'
% %         Orientation (1,:) char = 'xz' %{mustBeMember(Orientation,["xz","yz","xy"])}="xz"
        AxesUnits (1,1) string {mustBeMember(AxesUnits,["wavelengths","mm"])}="wavelengths"
        Position (1,:) double = [1 1 512 512]
        ReferencePt (1,3) double
        pdelta (1,1) double = 0.25
        Colormap (256,3) double = gray(256)
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
    methods(Static)
        function vDisplayWindow = QUPS(scan,kwargs)
            arguments
                scan (1,1) ScanCartesian
		kwargs.?VSXDisplayWindow
            end
            args = namedargs2cell(kwargs);
            vDisplayWindow = VSXDisplayWindow( ...
                'ReferencePt', [scan.xb(1), scan.yb(1), scan.zb(1)], ...
                'pdelta', scan.dx, ... 0.35
                'Position', [1, 1, scan.nx, scan.nz], ...
                args{:} ...
            ); 
        end
    end
end
