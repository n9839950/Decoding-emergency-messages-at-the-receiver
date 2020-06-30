%% EGB342 Assignment 2A Template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Run the GenerateAssignment2AData.m script if you have not already done
%  so.
clear all, close all, clc % clearing and preparing a clean workspace.
%
% Please take advantage of the naming convention laid out in the brief
%
%==========================================================================
% Enter your code below this line:
%==========================================================================
% Student Name: Manan Patel     
% Student Number: n9839950 

%% 1.4 Create the signal waveform S1 = S1(t) and S0 = S0(t) with 100 data points

t0 =500e-9; % Time duration of the symbol
NumPts = 100; % Number of data points per symbol 
A = 2.5;
x=linspace(-t0,t0,NumPts+1);
x(end)= []; % Time vector for symbol 

S1_pos =  A*(1-(-x(1:(end/2))/t0)); 
S1_neg = A*(1-(x((end/2 +1):end)/t0)); 
S1 = [S1_pos S1_neg];            
S0 = -S1; % Antipodal Signal 

% plot of S0 and S1
figure(1)
subplot(2,1,1),plot(x,S0),title('S_0(t)')
ylabel('Voltage(V)'),xlabel('time(s)')
subplot(2,1,2)  
plot(x,S1),title('S_1(t)')
ylabel('Voltage(V)'),xlabel('time(s)')


%% 1.5 Find the average energy of the modulation scheme using the s1 and s0
t_sym = (x(2)-x(1)) ; % sampling time
E_avg =(sum((S1).^2)*t_sym); %average Energy 


%% 1.6 Write a MATLAB code to modulate your name and ID number
nameID_str = double('MANAN B PATEL N9893950');
nameID_bits = de2bi(nameID_str ,7 , 'left-msb');
t_msg = reshape(nameID_bits',1,154); 


%% 1.7 Plot the spectrum of the baseband modulation signal generated in part  and find the null to null bandwidth

N_sym = length(t_msg);    % number of data symbols
T_sym = 2*t0;            % Symbol Duration
t = linspace(0,N_sym*T_sym, N_sym*NumPts+1); 
t(end)=[] ; % Time vector of the message signal
fs= 1/(t(2)-t(1));

info_t_msg = t_msg;
info_t_msg(t_msg == 0) = -1; 
b_msg = zeros(1,length(t)); % create an upsampled message
b_msg(1:NumPts:end) = info_t_msg; %

% Replace impulses with corresponding symbol
msg_S = conv(b_msg, S1) ; % modulate signal
msg_S = msg_S(1:length(b_msg));

figure(2)
plot(t,msg_S); 
title('Time domiain')
xlabel('time(S)')
ylabel('Amplitude(V)')

nfft = 2^15; %dft length

f = linspace(-fs/2,fs/2,nfft+1);
f=f(1:end-1); %frequency vector
Pxx1 = abs(fft(msg_S,nfft)/nfft);
figure(3)
plot(f,10*log10(fftshift(Pxx1)));
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density')
title('Transmit spectrum');


%% 1.8 Calculate Signal-to-Noise ratio (SNR) at the destinations
%null_BW = 2.1e6;  
%E_b = E_avg;
L_A = 15*4;
L_B = 25*4;
N_0 = 1e-13;

% transmitted power
tx_p = 10*log10((E_avg)/T_sym);  

% calculating power received 
pr_xa_db = (tx_p - L_A);
pr_xa = 10^(pr_xa_db/10);
pr_xb_db = (tx_p - L_B); 
pr_xb = 10^(pr_xb_db/10);

%calculating SNR in dB 
SNR_A = 10*log10((pr_xa*T_sym)/N_0) + 10*log10(2) - 10*log10(NumPts);
SNR_B = 10*log10((pr_xb*T_sym)/N_0) + 10*log10(2) - 10*log10(NumPts);
 

%% 1.9 Generate and plot these noisy signals received at locations A and B

%adding AWGN to the signal 
noise_a = awgn(msg_S, SNR_A,'measured');
noise_b = awgn(msg_S, SNR_B,'measured'); 

% plot 
figure(4)
subplot(2,1,1),plot(t,noise_a)
title('Received Signal at A')
xlabel('time(S)'),ylabel('Amplitude(V)')
subplot(2,1,2), plot(t,noise_b)
title('Received Signal at B')
xlabel('time(S)'),ylabel('Amplitude(V)')

%% 1.11 Create and plot the impulse response hopt(t) of the optimum receiver filter matched to detect symbol S1.
%Impulse response of matched filter
h_opt_s1 = conj(fliplr(S1)); 
figure(5)
plot(x,h_opt_s1)
xlabel ('Time(s)')
ylabel ('Amplitude(V)')
title('Impulse Response')
h_opt_s0 = conj(fliplr(S0)); 

%% 1.12 Filter the clean and noisy messages generated in Part 1.4 using the filter() command

Rx = filter(h_opt_s1, 1, msg_S)/NumPts;
Rx_noise_a = filter(h_opt_s1, 1, noise_a)/NumPts;
Rx_noise_b = filter(h_opt_s1, 1, noise_b)/NumPts;

%% 1.13 Plot and compare the three filtered messages

figure(6) 
subplot(3,1,1) 
plot(t, Rx)
xlabel('Time(s)'), ylabel('Amplitude(V)'), title('Filtered message')

subplot(3,1,2) 
plot(t, Rx_noise_a)
xlabel('Time(s)'), ylabel('Amplitude(V)'), title('Filtered message at A')

subplot(3,1,3) 
plot(t, Rx_noise_b)
xlabel('Time(s)'), ylabel('Amplitude(V)'), title('Filtered message at B')


%% 1.14 Decode the output and display the received message

zm_a = h_opt_s0(NumPts:NumPts:end);
decode_zm = sign(zm_a);
decode_zm(decode_zm == -1)= 0;

zm_x = Rx_noise_a (NumPts: NumPts :end); % amplitude required for filtered message  
decode_zm_x = sign(zm_x); % detecting symbol 
decode_zm_x (decode_zm_x == -1) = 0 ;

decode_msgA = reshape(decode_zm_x,[7,22])' ;
decode_msgA = bi2de(decode_msgA,'left-msb');
decode_msgA = reshape(decode_msgA,[1,22]);
decode_msgA = char(decode_msgA);  

zm_x_b = Rx_noise_b (NumPts: NumPts :end); % amplitude required for filtered message  
decode_zm_x = sign(zm_x_b); % detecting symbol 
decode_zm_x  (decode_zm_x == -1) = 0 ; 

decode_msgB = reshape(decode_zm_x,[7,22])' ;
decode_msgB = bi2de(decode_msgB,'left-msb');
decode_msgB = reshape(decode_msgB,[1,22]);
decode_msgB = char(decode_msgB);  

% number of bit error  
bit_err_a = abs(decode_zm - decode_zm_x)/2;
bit_erra = sum(bit_err_a); % for signal at loc A 
bit_err_b = abs (decode_zm - decode_zm_x)/2;
bit_errb = abs(bit_err_b); % for signal at loc B 


%% PART 2 Decoding Emergency message at the receiver

%% 2.1 Find the impulse response h_pre(t) of the optimum receiver

load A2AData.mat % load the data
h_preamble = conj(fliplr(preamble)) ;

figure(7)
plot(t_preamble, preamble)
xlabel('Time(s)'),ylabel('Amplitude'),title('Preamble')
plot(t_preamble, h_preamble)
xlabel('Time(s)'),ylabel('Amplitude'),title('Optimum Receiver')


%% 2.2 Plot the signal respresenting the noisy message sequence 

msg_rec = filter(h_preamble,1,rx)/NumPts;
plot(t_signal,msg_rec);

[Peak_value , EndLoc] = max(msg_rec);
msg_str = EndLoc +1;

figure(8)
rx_msg_t = t_signal(msg_str: end);
rx_msg_emg = rx(msg_str: end );
plot(rx_msg_t, rx_msg_emg);
xlabel('Time(s)'),ylabel('Amplitude'),title('Signal')


%% 2.3 Design the impulse response(s) needed for matched filter
rx_filtered_msg = filter(h_opt_s1,1,rx_msg_emg)/NumPts ;
figure(9) 
plot(rx_msg_t,rx_filtered_msg )

dcd_sec = find(t_signal == 0.0007975);
L_sqt = rx_msg_t (1:dcd_sec);
L_sq = rx_msg_emg(1:dcd_sec);
figure(10)
plot(L_sqt,L_sq) 
xlabel('Time(s)'),ylabel('Amplitude'),title('Filtered Signal')


%% 2.4 Implement a function (MF_receiver) for the matched filter receiver. The function header is to conform to the following syntax.

decode_emg_msg = MF_receiver(L_sq, NumPts,h_opt_s1);





