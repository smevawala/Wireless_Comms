%% Wireless Comms Mini Matlab OFDM
%Neema Aggarwal
%Shivam Mevawala
%Nico Castro

close all;

SNR = -4:1:8; %list of SNR values to run algorithm
%intialize vecs
BERc=zeros(length(SNR));
tblen =16; %will handle delay for convolution coder

n=3072; %msg length, must be mult of 6,4, and 128
m=2; %QPSK is 2-QAM whick is BPSK
nsubc=64; % number of subcarriers

%use the SNR to calculate EbNo
EbNo_c = SNR -10*log10(log2(m));


delayVector = [0 1 2 3 4] * 1e-5; % Discrete delays of five-path channel (s)
gainVector = linspace(0, 8, length(delayVector));
taps=5;
rchan_sel=rayleighchan(1e-5,0, delayVector, gainVector);%Set up the
%channel fading object with delay and gain vecs
% rchan_sel.StoreHistory = 1;
rchan_sel.StorePathGains = 1;

%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));

%loop over SNR values
for k=1:length(SNR)
X_bin=randi([0 m-1],n,1);

%modulate
Yc=qammod(X_bin,m);
ny=length(Yc)/nsubc;

% parallelize
Pc = reshape(Yc, [], nsubc);

% IFFT
ifft_sig=ifft(Pc,nsubc,2);

% Adding Cyclic Extension to enable circular conv
cext_data=[ifft_sig(:,(nsubc-15:nsubc)) ifft_sig];

%filter through chan
% Fc = zeros(size(cext_data,1),nsubc+16);
% PG = zeros(size(cext_data,1), taps);
for kkk=1:size(cext_data,1)
    Fc(kkk, :)=filter(rchan_sel,cext_data(kkk,:));
    PG(kkk,:)=rchan_sel.PathGains;
end

%add noise
Ac=awgn(Fc, SNR(k),'measured');

%Removing Cyclic Extension
rxed_sig=Ac(:,17:end);

% FFT
ff_sig=fft(rxed_sig,nsubc,2);

% ZF EQ
PG_fft = fft(PG, nsubc, 2);
Ye = ff_sig./PG_fft;

% serialize
Sc=reshape(Ye,[],1); 

%demodulate
Zc=qamdemod(Sc,m);

%calculate bit error rate
ber=biterr(Zc,X_bin)/(n);
BERc(k)=ber;
end

% plot BER
figure

semilogy(EbNo_c, BERc,'kx-');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plot')
legend('ZF', 'MMSE')


