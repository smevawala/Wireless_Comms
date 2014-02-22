%% Wireless Comms mini Matlab 2
%Neema Aggarwal
%Shivam Mevawala
%nicobitch

close all;
SNR = -4:1:20; %list of SNR values to run algorithm
%intialize vecs
BERc=zeros(length(SNR));

n=3072; %number of samples
m=4; %QPSK is 4-QAM

rchan=rayleighchan(1e-5,1e4);
rchan.StoreHistory = 1;
rchan.StorePathGains = 1;
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
    A=filter(rchan,Y);
    A = awgn(A, SNR(k),'measured');
    %demodulate
    Z=qamdemod(A./rchan.PathGains.',m);
    %calculkate bit error rate
    ber=biterr(Z,X)/(2*n);
%     ber = sum(Z ~= X)/length(X);
    BERc(k)=ber;

end

%plots



figure
% SPECT = distspec(trellis,7);
semilogy(EbNo,berfading(EbNo,'qam',4,1),'m-');
hold on;
semilogy(EbNo, BERc,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots (Convolutional Coder)')
legend('theoretical', 'actual')


