function [xdc, Trans] = TransducerVerasonics(name, units)
% TransducerVerasonics - Construct a Verasonics Transducer by name
% 
% xdc = TransducerVerasonics(name) calls `computeTrans` to construct a QUPS
% Transducer `xdc` from the given `name`.
%
% [xdc, Trans] = TransducerVerasonics(name) additionally returns the
% Vantage struct `Trans`.
%
% [...] =  TransducerVerasonics(name, units) additionally sets the units to
% either "mm" or "wavelengths". The default is "mm".
%
% Example:
% [xdc, Trans] = TransducerVerasonics("L12-3v"),
% 
% See also Transducer.Verasonics computeTrans
arguments
    name (1,1) string {mustBeMember(name, ["L7-4","L10-5","L10-4v","L11-4v","L11-5","L11-5v","L11-5gH","L12-3v","L12-5 38mm","L12-5 50mm","L22-8v","L22-14v","L22-14v LF","L22-14vX","L22-14vX LF","L35-16vX","L38-22v","L39-21gD","GE9LD","4DL7","GEC1-6D","GE4CD","GEIC5-9D","GEL3-12D","GEM5ScD","CL10-5","CL15-7","C4-2","C5-2","C5-2v","C5-2gH","C7-4","C8-4V","C8-5","C9-5ICT","P3-2","P4-1","P4-2","P4-2v","P4-2gH","P5-3","P6-3","P7-4","RC6gV","P5-64vN","P5-128vN","IP-104","IP-105","H-101","H-104","H-106","H-301","H-302","H-313","H-313A","Matrix1024-3","Matrix1024-8","UTA 1024-MUX","UTA 256-Direct","Adapter Embedded S.T.E","260 ZIF Backshell Kit","408 ZIF Backshell Kit","260 ZIF S.T.E. Fixture","408 ZIF S.T.E. Fixture","260 ZIF Calibration STE"])}
    units (1,1) {mustBeMember(units, ["mm", "wavelengths"])} = "mm"
end

Trans = computeTrans(struct("name", char(name), 'units', char(units)));
xdc = Transducer.Verasonics(Trans);




