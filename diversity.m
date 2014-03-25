%% Wireless Comms mini Matlab 2
%Neema Aggarwal
%Shivam Mevawala
%Nicolas Castro

close all;
clear all;
SNR = -4:2:20; %list of SNR values to run algorithm
%intialize vecs
BER_reg=zeros(length(SNR));
BER_mrrc2=zeros(length(SNR));
BER_mrrc4=zeros(length(SNR));

n=3072; %number of samples
m=2; %BPSK is 2-QAM

% delayVector = [0 1 2 3 4] * 1e-5; % Discrete delays of four-path channel (s)
% gainVector = [0 -4 -6 -9 -14];

rchan_flat=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.
rchan_flat2=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat2.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat2.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

rchan_flat3=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat3.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat3.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.
rchan_flat4=rayleighchan(1e-5,1e4); %Set up the channel fading object
rchan_flat4.StoreHistory = 1; %If this value is 1, channel state information
%needed by the channel visualization tool is stored. 
rchan_flat4.StorePathGains = 1;%If set to 1, the complex path gain vector is
%stored.

%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));


%loop over SNR values
for k=1:length(SNR)
    %generate a random vector of 4 symbols
    X=randi([0 m-1],1,n);
    %modulate
    Y=qammod(X,m);
    %add noise and rayleigh fading
    A=filter(rchan_flat,Y);
    A2=filter(rchan_flat2,Y);
    A3=filter(rchan_flat3,Y);
    A4=filter(rchan_flat4,Y);    
    Ag = awgn(A, SNR(k),'measured');
    Ag2 = awgn(A2, SNR(k),'measured'); 
    Ag3 = awgn(A3, SNR(k),'measured'); 
    Ag4 = awgn(A4, SNR(k),'measured'); 
    Ae=Ag.*conj(rchan_flat.PathGains.');
    Ae2=Ag2.*conj(rchan_flat2.PathGains.');
    Ae3=Ag3.*conj(rchan_flat3.PathGains.');
    Ae4=Ag4.*conj(rchan_flat4.PathGains.');
    %demodulate
    Z=qamdemod(Ae,m);
    Z_mrrc2=qamdemod(Ae+Ae2,m);
    Z_mrrc4=qamdemod(Ae+Ae2+Ae3+Ae4,m);
    %calculate bit error rate
    BER_reg(k)=biterr(Z,X)/(2*n);
    BER_mrrc2(k)=biterr(Z_mrrc2,X)/(2*n);
    BER_mrrc4(k)=biterr(Z_mrrc4,X)/(2*n);
end

%plots
figure 
semilogy(EbNo, BER_reg,'kx');
hold on;
semilogy(EbNo, BER_mrrc2,'go');
semilogy(EbNo, BER_mrrc4,'rx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots- Flat Fading')
legend('no diversity', 'MRRC 2 RX', 'MRRC 4 RX')

