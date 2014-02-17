hConvEnc = comm.ConvolutionalEncoder(poly2trellis(7, [171 133]));
hConvEnc.PuncturePatternSource = 'Property';
hConvEnc.PuncturePattern = [1;1;0;1;1;0];
hMod = comm.QPSKModulator;
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)',...
  'SignalPower', 1, 'SamplesPerSymbol', 1);

hDMod = comm.QPSKDemodulator;
hVitDec = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
  'InputFormat', 'Unquantized');

hVitDec.PuncturePatternSource =  'Property';
hVitDec.PuncturePattern = hConvEnc.PuncturePattern;
hVitDec.TracebackDepth = 96;
hErrorCalc = comm.ErrorRate('ReceiveDelay', hVitDec.TracebackDepth);

EbNoEncoderInput = 2:.5:5; % in dB

EbNoEncoderOutput = EbNoEncoderInput + 10*log10(log2(m)*3/4);

frameLength = 3000; % this value must be an integer multiple of 3
targetErrors = 300;
maxNumTransmissions = 5e6;

BERVec = zeros(3,length(EbNoEncoderOutput)); % Allocate memory to store results
for n=1:length(EbNoEncoderOutput)
  reset(hErrorCalc)
  reset(hConvEnc)
  reset(hVitDec)
  hChan.EbNo = EbNoEncoderOutput(n); % Set the channel EbNo value for simulation
  while (BERVec(2,n) < targetErrors) && (BERVec(3,n) < maxNumTransmissions)
    % Generate binary frames of size specified by the frameLength variable
    data = randi([0 1], frameLength, 1);
    % Convolutionally encode the data
    encData = step(hConvEnc, data);
    % Modulate the encoded data
    modData = step(hMod, encData);
    % Pass the modulated signal through an AWGN channel
    channelOutput = step(hChan, modData);
    % Pass the real part of the channel complex outputs as the unquantized
    demData = step(hDMod,channelOutput); 
    % input to the Viterbi decoder.
    decData = step(hVitDec,demData);
    % Compute and accumulate errors
    BERVec(:,n) = step(hErrorCalc, data, decData);
  end
end


dist = 5:11;
nerr = [42 201 1492 10469 62935 379644 2253373];
codeRate = 3/4;
bound = nerr*(1/6)*erfc(sqrt(codeRate*(10.0.^((2:.02:5)/10))'*dist))';

berfit(EbNoEncoderInput,BERVec(1,:)); % Curve-fitted simulation results
hold on;
semilogy((2:.02:5),bound,'g'); % Theoretical results
legend('Empirical BER','Fit for simulated BER', 'Theoretical bound on BER')
axis([1 6 10^-6 10^-1])


