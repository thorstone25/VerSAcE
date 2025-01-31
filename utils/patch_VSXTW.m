% patch_VSXTW
%
% patch_VSXTW enables the `VSXTW.computeTWWaveform()` method. A copy of
% Vantage must be on the path and activated.
% 
% See also: computeTWWaveform VSXTW

% only run the 1st time
ismth = ismember("computeTWWaveform", methods("VSXTW"));
if ismth, warning("VSXTW.computeTWWaveform() is available: continuing."); return; end 

w = whitespacePattern(0,Inf); % alias: whitespace maybe 
str = readlines(which("computeTWWaveform")); % grab Vantage files

% fix the nonlin subfunction (add 'Resource')
pat = "nonlin"+whitespacePattern(0,Inf)+"("+wildcardPattern+")"; % pattern
ins = ", Resource"; % insert
for i = find(contains(str, pat))' % matching lines
    e = strfind(str(i), ")"); % ending paren
    e(e >= min([strfind(str(i), "%"), Inf])) = []; % ignore comments
    str(i) = insertBefore(str(i), e(end), ins); % insertion
end

% % % workspace variables
% assert that TPC will exist
pat = "evalin('base', 'exist(''TPC'', ''var'')')";
str = replace(str, pat, "~isempty(TPC)");

% remove the "x = evalin('base', 'x')" 
var = ("Resource" | "Trans" | "TPC"); % variables
pat = optionalPattern(var +w+"="+w) + ("evalin" | "assignin")+whitespacePattern(0,Inf)...
    +"("+wildcardPattern+"base"+wildcardPattern+var+wildcardPattern+")"+optionalPattern(";"); % pattern
str = replace(str, pat, ""); % delete these lines

% set returnTW always true
pat = "returnTW"+w+"="+w+"0";
str = replace(str, pat, "returnTW = true");

% % % struct -> class conversions
% sync TWout/TWin
pat = "TWout"+w+"="+w+"TWin";
i = find(contains(str, pat));
str(i) = replace(str(i), "TWin", "copy(TWin)"); % explicit copy for classes
str = replace(str, "TW"+("out"|"in"), "TW"); % TWout / TWin -> TW

% replace "TW(i).TriLvlWvfm_Sim" with "TW(i).TriLvlWvfm" 
mch = ".TriLvlWvfm";
pat = "TW"+optionalPattern("("+alphanumericsPattern(0,Inf)+")")+mch+"_Sim";
i = find(contains(str, pat));
str(i) = replace(str(i), mch+"_Sim", mch);

% isprop <- isfield for classes
var = ("Resource" | "TPC" | "TW"+optionalPattern("in"|"out")); % variables
mch = "isfield("+var;
i = find(contains(str, mch));
str(i) = replace(str(i), "isfield", "isprop");

% % % entry
% change the signature
pat = "function"+wildcardPattern+"computeTWWaveform("+wildcardPattern+")";
rep = "function [TW, Trans, TPC] = computeTWWaveform(TW, Trans, Resource, TPC)";
str = replace(str, pat, rep);

% add the argument validation
j = find(contains(str, "%% Initialization"), 1, 'first');
i = "    "; % indent
ins = [ ...
"arguments", ...
i + "TW VSXTW", ...
i + "Trans (1,1) struct", ...
i + "Resource (1,1) VSXResource", ...
i + "TPC VSXTPC {mustBeScalarOrEmpty} = VSXTPC.empty", ...
"end    ", ...
]';
str = cat(1,str(1:j), ins, str(j+1:end));

% % % MATLAB ... not C ...
str = replace(str, "fprintf(2,"+w, "error(");

% % % write as a method for VSXTW
% fld = fullfile(fileparts(which('VSXTW.m')), "+VSXTW"); % folder
% mkdir(fld); writelines(str, fullfile(fld, "computeTWWaveform2.m")); % write

% % % insert into VSXTW
i = find(contains(str, "%% Heuristic Nonlinearity"));
fmth = str(1:i-1); % method
fsub = str(i:end); % subfunction

flnm = which("VSXTW.m");
str_tw = readlines(flnm); % read the classdef file
j = find(contains(str_tw, "methods"), 1, 'last'); % start of a methods block
i = extract(str_tw(j), lineBoundary('start')+w); % get the indentation whitespace
while(startsWith(str_tw(j-1),w+"%")), j = j - 1; end % avoid comment blocks

ins = i+cat(1, "","methods","    "+fmth,"end",""); % insertion text (method block)
str_tw = cat(1, str_tw(1:j-1), ins, str_tw(j:end)); % insert
str_tw = cat(1, str_tw, fsub); % append the subfunction
writelines(str_tw, flnm); % write back out to file