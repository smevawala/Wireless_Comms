%% Wireless Comms mini Matlab 2
%Neema Aggarwal
%Shivam Mevawala
%nicobitch

close all;
SNR = -4:2:20; %list of SNR values to run algorithm
%intialize vecs
BER_flat=zeros(length(SNR));
BER_sel=zeros(length(SNR));

n=10000; %number of samples
m=4; %QPSK is 4-QAM
% delayVector = 1.0e-004 * [0 0.0400 0.0800 0.1200];  % Discrete delays of
                                                    % four-path channel (s)
% gainVector = [0 -3 -6 -9];

delayVector = [0 1 2 3 4] * 1e-5;
gainVector = [0 -4 -6 -9 -14];

rchan_flat=rayleighchan(1e-5,1e4);
rchan_flat.StoreHistory = 1;
rchan_flat.StorePathGains = 1;

rchan_sel=rayleighchan(1e-5,1, delayVector, gainVector);
rchan_sel.StoreHistory = 1;
rchan_sel.StorePathGains = 1;

% chan = [.986; .845; .237; .123+.31i]
%use the SNR to calculate EbNo for the normal sytem and the convolutional
%coder
EbNo = SNR -10*log10(log2(m));
% eq = lineareq(8, lms(0.001));
eq = dfe(5, 5, rls(.99));
eq.SigConst=qammod(0:3,4);
% eqlms.RefTap = 4;
bers=zeros(1,10);
%loop over SNR values
for k=1:length(SNR)
    k
    %generate a random vector of 4 symbols
    X=randi([0 m-1],1,n);
    %modulate
    Y=qammod(X,m);
    %add noise
%     Y =awgn(Y, SNR(k),'measured');
    A=filter(rchan_flat,Y);

%     As = filter(chan, 1, Y);
    A = awgn(A, SNR(k),'measured');
    Ae=A./rchan_flat.PathGains.';
    %demodulate
    Z=qamdemod(Ae,m);

    %calculate bit error rate
%     ber = sum(Z ~= X)/length(X);
    BER_flat(k)=biterr(Z,X)/(2*n);
    for kk=1:10
        kk
        As=filter(rchan_sel,Y);
        As = awgn(As, SNR(k),'measured');
        Ase=equalize(eq,As,Y(1:1000));
        Zs=qamdemod(Ase,m);
        bers(kk)=biterr(Zs,X)/(2*n);
    end
    BER_sel(k)=mean(bers);
end

%plots



figure
% SPECT = distspec(trellis,7);
semilogy(EbNo,berfading(EbNo,'qam',4,1),'m-');
hold on;
semilogy(EbNo, BER_flat,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots- Frequency Fading')
legend('theoretical', 'actual')


figure
% semilogy(EbNo,berfading(EbNo,'qam',4,1),'m-');
% hold on;
semilogy(EbNo, BER_sel,'kx');

xlabel('EbNo (dB)')
ylabel('BER')
title('Waterfall Plots- Frequency Selective')
legend('theoretical', 'actual')


