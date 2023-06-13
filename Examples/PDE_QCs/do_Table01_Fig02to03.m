function do_Table01_Fig02to03()
% this lists table 2 and plots Fig 2 and Fig 3 
clear all;
close all;
% change RBFs: For Gaussian: type:'g' with par=any values
%              FOR MQ:  type='mq' with par= 0.5
%              FOR IMQ: type='mq' with par= -0.5
%              For W2:  type='w' with par= 2 **see below note.
%              For Matern with \nu=5:  type='ms' with par= 5
type='g'; 
par=0.5;

n=[500;1000;2000]; %list of number of collocation points
nt=900; % number of test points

% F1 is the function with singularity near the boundary
% F2-F6 are six Frankes' benchmark functions
%only F1 and F2 are presented in the paper
% "CASE=1" indicates F1.
CASE=2;  

% main code
error_list(type, par, n, nt, CASE)

%** note on Wendland RBF: When using Wendland, change "scaleList" in 
%"Ex3_errList" function to "scaleList = RBFscalestart:0.1:20;"
%since the Wendland RBF needs higher scale values to converge.
% other RBFs are fine with default "scaleList = RBFscalestart:0.01:1;"

%%% main code and instruction can be found on 
% "https://num.math.uni-goettingen.de/schaback/research/papers/MPfKBM.pdf"

