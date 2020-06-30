%% EGB342 Assignment 2A - Data Generator
%  Use this file to generate the Data and signals that you will
%  be using in Assignment 2A.
%  A single .mat file will be generated for you which includes
%  all of the signals that you will be using in this assignment.

%  - Fill in your student number on Line 21.
%  - Run the code.

% Preparing MATLAB workspace
clear, clc, close all;

%% %%%% This is an individual assignment %%%% %%


% Enter Your student ID. 
% Discard the leading 'n'.
% If your number begins with '0', also discard the '0'.
% E.g. n08857181 becomes 8857171, 

n = 9839950; % Enter your student number here

%% STOP %% Do not change anyting below this line until Part 2.7
%
EbN0 = 7; % (in dB) Later in Part 2.9 You will change the this value. %
%
%
%% Do not change
%
% Generating your signals.
[rx, t_signal, preamble, t_preamble] = DataGen(n,EbN0);
save A2AData.mat rx t_signal preamble t_preamble;
fprintf('Your .mat file has been saved to:\n\r');
disp(pwd)
fprintf('In the file "A2AData.mat".\n\r');
fprintf('Use the "load A2AData.mat" Command to load the saved variables into the workspace.\n\r');
% You only need to run this file ONCE to generate your signals.