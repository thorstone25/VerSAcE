classdef VSXPData < matlab.mixin.Copyable
    properties
        Coord (1,1) string {mustBeMember(Coord, ["rectangular",...
                                                 "polar",...
                                                 "spherical"])} = "rectangular"
        PDelta (1,:) double {mustBeFinite}
        Size (1,:) double
        Origin (1,:) double
        Region (1,:) struct = struct('Shape',struct('Name','PData'),'numPixels',[],'PixelsLA',[]);
    end
    methods
        function obj = VSXPData(kwargs)
            arguments
                kwargs.?VSXPData
            end
            for f = string(fieldnames(kwargs))'
                obj.(f) = kwargs.(f); 
            end
        end
    end
    methods(Static)
        function vPData = QUPS(scan)
    	    arguments
    	        scan (1,1) ScanCartesian
    	    end
    	    vPData = VSXPData(...
        	    'PDelta', [scan.dx, 0, scan.dz],...
        	    'Size'  , [scan.nz, scan.nx, scan.ny],...
        	    'Origin', [scan.xb(1), scan.yb(1), scan.zb(1)] ...
        	    );
        end
    end
end
