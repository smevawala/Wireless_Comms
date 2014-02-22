%% Wireless Comms mini Matlab 2
%Neema Aggarwal
%Shivam Mevawala
%nicobitch

close all;
SNR = -4:1:20; %list of SNR values to run algorithm
%intialize vecs
BER_flat=zeros(length(SNR));
BER_freq=zeros(length(SNR));

n=3072; %number of samples
m=4; %QPSK is 4-QAM

rchan_flat=rayleighchan(1e-5,1e4);
rchan_flat.StoreHistory = 1;
rchan_flat.StorePathGains = 1;
%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));

%loop over SNR values
for k=1:length(SNR)
    %generate a random vector of 4 symbols
    X=randi([0 m-1],1,n);
    %modulate
    Y=qammod(X,m);
    %add noise
%     Y =awgn(Y, SNR(k),'measured');
    A=filter(rchan_flat,Y);
    A = awgn(A, SNR(k),'measured');
    %demodulate
    Z=qamdemod(A./rchan_flat.PathGains.',m);
    %calculkate bit error rate
    ber=biterr(Z,X)/(2*n);
%     ber = sum(Z ~= X)/length(X);
    BER_flat(k)=ber;

    
end

%plots



figure
% SPECT = distspec(trellis,7);
semilogy(EbNo,berfading(EbNo,'qam',4,1),'m-');
hold on;
semilogy(EbNo, BER_flat,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots')
legend('theoretical', 'actual')


