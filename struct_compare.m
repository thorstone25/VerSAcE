import matlab.unittest.TestCase
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.constraints.StructComparator
import matlab.unittest.constraints.NumericComparator
import matlab.unittest.constraints.LogicalComparator
import matlab.unittest.constraints.StringComparator

%% preprocess
% activate;
% addpath ~/verasonics/Ameya/VSXOOD/
% addpath ~/verasonics/Ameya/vsxTest/
% addpath(genpath('~/verasonics/Ameya/qups/'));
% addpath(genpath('~/verasonics/Vantage-4.8.4-2305101400'));
sequence_types_test;
vsx__compare_debugger;

%%
% sub field
fsub = struct('a', 1, 'b', "2", 'c', 'three');

% cases to test for equality
% exp = struct("f1",zeros(1,10),"f2",'a',"f3",       {'b','c'} , "fs", fsub);
st1 = struct("f1",zeros(1,10),"f2","a","f3",       {'b','c'} , "fs", fsub); % f2 is string
st2 = struct("f1",zeros(1,10),"f2",'a',"f3",string({'b','c'}), "fs", fsub); % f3 is string
% cpy = exp; % just a copy

for f = string(fieldnames(vsx))'
   exp = og.(f);
   cpy = vsx.(f);



% simple tests
testCase = TestCase.forInteractiveUse;
testCase.verifyThat(cpy,IsEqualTo(exp));
% testCase.verifyThat(cpy,IsEqualTo(exp)); % ,"Using",StructComparator(NumericComparator)));
% testCase.verifyThat(st1,IsEqualTo(exp)); % ,"Using",StructComparator(NumericComparator)));
% testCase.verifyThat(st2,IsEqualTo(exp)); % ,"Using",StructComparator(NumericComparator)));

% more complicated test
% [cpy.f1] = deal(zeros(0,10)); % wrong size, but ignore
% [cpy.f2] = st1.f2; % f2 string instead of char
testCase.verifyThat(cpy,IsEqualTo(exp, ...
    "Using", StructComparator( ...
    [NumericComparator, LogicalComparator, StringComparator], ...
    "Recursively", true, ...
    "IgnoringFields", ["Parameters"] ...
    )...
    ));
end