file1 = './MatFiles/L12-3vFlashAngles';
file2 = './MatFiles/L12-3vFlashAngles_test';
% file3 = './MatFiles/L12-3v_192RyLns_Hadamard.mat';

og = load(file1);
vsx = load(file2);
% hadamard = load(file3);

mydata = og;

%---------------------

AxesUnit = vsx.AxesUnit;

DwHeight = 543;%  
DwWidth = 637;%

Event = vsx.Event;
Media = vsx.Media;

MinMaxVal = [11.0382 51.7415 33.1146];%

P = vsx.P;
PData = vsx.PData;
Process = vsx.Process;
Receive = vsx.Receive;
Recon = vsx.Recon;
ReconInfo = vsx.ReconInfo;
Resource = vsx.Resource;

ScrnSize = [1,1,1920,1080]; %

SeqControl = vsx.SeqControl;
TGC = vsx.TGC;
TW = vsx.TW;
TX = vsx.TX;
Trans = vsx.Trans;
UI = vsx.UI;


angle = 0.3142;%
dtheta = 0.1047;
frameRateFactor = 2;
i = 40;
j = 13;
k = 546;
maxAcqLength = 294;
n = 601;
na = 7;
nsc = 45;
pers = 20;
startAngle = -0.3142;


filename = 'L12-3vFlashAngles_demo';

for f = string(fieldnames(mydata))'
    save(filename, (f), '-append');
end

%--------------------------

vsx__compare_debugger;

filename = 'MatFiles/L12-3vFlashAngles_demo';
VSX;
