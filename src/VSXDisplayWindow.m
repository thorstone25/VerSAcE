classdef VSXDisplayWindow < matlab.mixin.Copyable
    properties
        Type (1,1) string {mustBeMember(Type,["Verasonics", "Matlab"])} = "Verasonics"
        Title (1,1) string = 'Window title'
        mode (1,1) string {mustBeMember(mode,["2d"])} = '2d'
        Orientation (1,1) string {mustBeMember(Orientation,["xz","yz","xy"])} = "xz"
        AxesUnits (1,1) string {mustBeMember(AxesUnits,["wavelengths","mm"])}="wavelengths"
        Position (1,:) double {mustBeInteger} = [1 1 512 512]
        ReferencePt (1,3) double
        pdelta (1,1) double = 0.25
        Colormap (256,3) double = gray(256)
        numFrames (1,1) double = 1
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
            ps = get(0, 'screensize'); % (primary) monitor size
            dp = min(scan.dx, scan.dz); % set spacing to smallest difference
            sz = round(range([scan.xb; scan.zb],2) ./ dp)'; % set size based on image range
            scl = 0.5 * min(ps(3:4) ./ sz);  % monitor scaling
            [dp, sz] = deal(dp ./ scl, sz * scl); % scale to monitor size
            vDisplayWindow = VSXDisplayWindow( ...
                'ReferencePt', [scan.xb(1), scan.yb(1), scan.zb(1)], ... upper left corner (2D, cartesian)
                'pdelta', dp, ... 0.35
                'Position', [1 1 round(sz)], ...
                args{:} ...
                );
        end
    end
end
