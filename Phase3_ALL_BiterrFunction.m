% Phase 3: Basic error control coding to improve the performance

clear all; close all; clc; 

% Var Declarations
fc = 10000;                                     % carrier freq
fs = 16*fc;                                     % sample freq
bbDataRate = 1000;                              % baseband DataRate 1kbps
nBits = 1024;                                   % data bits
SNRdb = -15:1:15;                               % SNR in db
SNR = (10.^(SNRdb/10));                         % convert SNR from db to dec
samplePeriod = fs/bbDataRate;
encoded_nBits = nBits/4*7;                      % Data bits for Cyclic/Hamming/Linear Code (7,4)

sigLen = fs*nBits/bbDataRate +1;                % signal length
sigLen_encoded = fs*encoded_nBits/bbDataRate +1; % Encoded signal length
amp = 2;                                        % amplitude
t = 0: 1/fs : nBits/bbDataRate;                 % time intervals --not encoded
t_encoded = 0: 1/fs : encoded_nBits/bbDataRate;  % time intervals --encoded

carrier = amp .* cos(2*pi*fc*t);
carrierLow_encoded = amp .* cos(2*pi*fc*t_encoded);             % for BFSK 1st lower freq carrier
carrierHigh_encoded = amp .* cos(2*pi*5*fc*t_encoded);          % for BFSK 2nd higher freq carrier
carrier_encoded = amp .* cos(2*pi*fc*t_encoded);      %For Encoded Carrier

[b,a] = butter(6,0.2);                          % 6th Order low pass filter coefficients


% create arrays to store Bit Error Rates BER
OOKerrorArr = zeros(1,length(SNR));
BPSKerrorArr = zeros(1,length(SNR));
BFSKerrorArr = zeros(1,length(SNR));
linear_OOKerrorArr = zeros(1,length(SNR));
linear_BPSKerrorArr = zeros(1,length(SNR));
linear_BFSKerrorArr = zeros(1,length(SNR));
cyclic_OOKerrorArr = zeros(1,length(SNR));
cyclic_BPSKerrorArr = zeros(1,length(SNR));
cyclic_BFSKerrorArr = zeros(1,length(SNR));
hamming_OOKerrorArr = zeros(1,length(SNR));
hamming_BPSKerrorArr = zeros(1,length(SNR));
hamming_BFSKerrorArr = zeros(1,length(SNR));

runCycles = 20;

plotSNRdb = -15;                                          % change in multiple of 5 to match; for plotting

for i = 1 : length(SNR)                                 % loop for diff SNR values
        aveOOKerror = 0;
        linear_aveOOKerror = 0;
        linear_aveBFSKerror = 0;
        linear_aveBPSKerror = 0;
        cyclic_aveOOKerror = 0;
        cyclic_aveBFSKerror = 0;
        cyclic_aveBPSKerror = 0;
        hamming_aveOOKerror = 0;
        hamming_aveBFSKerror = 0;
        hamming_aveBPSKerror = 0;

    for j = 1 : runCycles                               % loop to cal ave error rate over runCycles
        data = round(rand(1,nBits));                    % data Gen
        genpoly = cyclpoly(7,4);
        parmat1 = cyclgen(7,genpoly);
        trt = syndtable(parmat1);
        cyclic_sig = encode (data,7,4,'cyclic/binary',genpoly); %cyclic data
        hamming_sig = encode (data,7,4,'hamming/binary'); %hamming data
        pol = cyclpoly(7,4);
        parmat2 = cyclgen(7,pol);
        genmat = gen2par(parmat2);
        linear_sig = encode (data,7,4,'linear/binary', genmat); %linear data
        % signal Gen
        sig = zeros(1,sigLen);
        sig_linear = zeros(1,sigLen_encoded);
        sig_hamming = zeros(1,sigLen_encoded);
        sig_cyclic = zeros(1,sigLen_encoded);
        for k = 1 : sigLen - 1
            sig(k) = data(ceil(k*bbDataRate/fs));
        end
        for k = 1 : sigLen_encoded - 1
            sig_linear(k) = cyclic_sig(ceil(k*bbDataRate/fs));
        end
        for k = 1 : sigLen_encoded - 1
            sig_hamming(k) = linear_sig(ceil(k*bbDataRate/fs));
        end
        for k = 1 : sigLen_encoded - 1
            sig_cyclic(k) = hamming_sig(ceil(k*bbDataRate/fs));
        end
        sig(sigLen) = sig(sigLen-1);
        sig_linear(sigLen_encoded) = sig_linear(sigLen_encoded - 1);
        sig_hamming(sigLen_encoded) = sig_hamming(sigLen_encoded - 1);
        sig_cyclic(sigLen_encoded) = sig_cyclic(sigLen_encoded - 1);
        %***** On-Off Keying (OOK) *****
           % modulation
        linear_OOKsig = sig_linear .* carrier_encoded;                        % if 1 = carrier; if 0 = 0;
        hamming_OOKsig = sig_hamming .* carrier_encoded;
        cyclic_OOKsig = sig_cyclic .* carrier_encoded;
        OOKsig = sig .* carrier;

            % noise Gen
        OOKsigPwr = bandpower(OOKsig);                  % sum of absolute square/length
        linear_OOKsigPwr = bandpower(linear_OOKsig);
        hamming_OOKsigPwr = bandpower(hamming_OOKsig);
        cyclic_OOKsigPwr = bandpower(cyclic_OOKsig);

        OOKnoisePwr = OOKsigPwr ./ SNR(i);
        OOKnoise = sqrt(OOKnoisePwr/2) .* randn(1,sigLen); %%%% align with Phase1
        
        linear_OOKnoisePwr = linear_OOKsigPwr ./ SNR(i);
        linear_OOKnoise = sqrt(linear_OOKnoisePwr/2) .* randn(1,sigLen_encoded);

        hamming_OOKnoisePwr = hamming_OOKsigPwr ./ SNR(i);
        hamming_OOKnoise = sqrt(hamming_OOKnoisePwr/2) .* randn(1,sigLen_encoded);

        cyclic_OOKnoisePwr = cyclic_OOKsigPwr ./ SNR(i);
        cyclic_OOKnoise = sqrt(cyclic_OOKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        OOKreceived = OOKsig + OOKnoise;
        linear_OOKreceived = linear_OOKsig + linear_OOKnoise;

        hamming_OOKreceived = hamming_OOKsig + hamming_OOKnoise;

        cyclic_OOKreceived = cyclic_OOKsig + cyclic_OOKnoise;

            % demodulation
        OOKdemod = OOKreceived .* (2.* carrier);        %coherent
        OOKdemod = filtfilt(b, a, OOKdemod);
        linear_OOKdemod = linear_OOKreceived .* (2.* carrier_encoded);        %coherent
        linear_OOKdemod = filtfilt(b, a, linear_OOKdemod);

        hamming_OOKdemod = hamming_OOKreceived .* (2.* carrier_encoded);        %coherent
        hamming_OOKdemod = filtfilt(b, a, hamming_OOKdemod);

        cyclic_OOKdemod = cyclic_OOKreceived .* (2.* carrier_encoded);        %coherent
        cyclic_OOKdemod = filtfilt(b, a, cyclic_OOKdemod);

            % sampling and decision logic
        OOK = samplingANDdecision(OOKdemod, samplePeriod, nBits, amp/2);

        linear_OOK = samplingANDdecision(linear_OOKdemod, samplePeriod, encoded_nBits, amp/2);

        hamming_OOK = samplingANDdecision(hamming_OOKdemod, samplePeriod, encoded_nBits, amp/2);

        cyclic_OOK = samplingANDdecision(cyclic_OOKdemod, samplePeriod, encoded_nBits, amp/2);

            %decoding Code
        OOK_linear = decode(linear_OOK,7,4,'linear/binary',genmat);

        OOK_hamming = decode(hamming_OOK,7,4,'hamming/binary');

        OOK_cyclic = decode(cyclic_OOK,7,4,'cyclic/binary',genpoly,trt);


        %***** Binary Phase Shift Keying (BPSK) *****
            % modulation
        linear_BPSKsig = sig_linear .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        linear_BPSKsig = linear_BPSKsig .* carrier_encoded;

        hamming_BPSKsig = sig_hamming .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        hamming_BPSKsig = hamming_BPSKsig .* carrier_encoded;

        cyclic_BPSKsig = sig_cyclic .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        cyclic_BPSKsig = cyclic_BPSKsig .* carrier_encoded;

            % noise Gen
        linear_BPSKsigPwr = bandpower(linear_BPSKsig);
        linear_BPSKnoisePwr = linear_BPSKsigPwr ./ SNR(i);
        linear_BPSKnoise = sqrt(linear_BPSKnoisePwr/2) .* randn(1,sigLen_encoded);

        hamming_BPSKsigPwr = bandpower(hamming_BPSKsig);
        hamming_BPSKnoisePwr = hamming_BPSKsigPwr ./ SNR(i);
        hamming_BPSKnoise = sqrt(hamming_BPSKnoisePwr/2) .* randn(1,sigLen_encoded);

        cyclic_BPSKsigPwr = bandpower(cyclic_BPSKsig);
        cyclic_BPSKnoisePwr = cyclic_BPSKsigPwr ./ SNR(i);
        cyclic_BPSKnoise = sqrt(cyclic_BPSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        linear_BPSKreceived = linear_BPSKsig + linear_BPSKnoise;

        hamming_BPSKreceived = hamming_BPSKsig + hamming_BPSKnoise;

        cyclic_BPSKreceived = cyclic_BPSKsig + cyclic_BPSKnoise;

            % demodulation
        linear_BPSKdemod = linear_BPSKreceived .* (2.* carrier_encoded);
        linear_BPSKdemod = filtfilt(b, a, linear_BPSKdemod);

        hamming_BPSKdemod = hamming_BPSKreceived .* (2.* carrier_encoded);
        hamming_BPSKdemod = filtfilt(b, a, hamming_BPSKdemod);

        cyclic_BPSKdemod = cyclic_BPSKreceived .* (2.* carrier_encoded);
        cyclic_BPSKdemod = filtfilt(b, a, cyclic_BPSKdemod);

            % sampling and decision logic
        linear_BPSK = samplingANDdecision(linear_BPSKdemod, samplePeriod, encoded_nBits, 0);

        hamming_BPSK = samplingANDdecision(hamming_BPSKdemod, samplePeriod, encoded_nBits, 0);

        cyclic_BPSK = samplingANDdecision(cyclic_BPSKdemod, samplePeriod, encoded_nBits, 0);

            %decoding Code
        BPSK_linear = decode(linear_BPSK,7,4,'linear/binary',genmat);

        BPSK_hamming = decode(hamming_BPSK,7,4,'hamming/binary');

        BPSK_cyclic = decode(cyclic_BPSK,7,4,'cyclic/binary',genpoly,trt);

        %***** Binary Frequency Shift Keying *****
            % modulation
        linear_BFSKsigLow = carrierLow_encoded .* (sig_linear == 0);
        linear_BFSKsigHigh = carrierHigh_encoded .* (sig_linear == 1);
        linear_BFSKsig = linear_BFSKsigLow + linear_BFSKsigHigh;

        hamming_BFSKsigLow = carrierLow_encoded .* (sig_hamming == 0);
        hamming_BFSKsigHigh = carrierHigh_encoded .* (sig_hamming == 1);
        hamming_BFSKsig = hamming_BFSKsigLow + hamming_BFSKsigHigh;

        cyclic_BFSKsigLow = carrierLow_encoded .* (sig_cyclic == 0);
        cyclic_BFSKsigHigh = carrierHigh_encoded .* (sig_cyclic == 1);
        cyclic_BFSKsig = cyclic_BFSKsigLow + cyclic_BFSKsigHigh;

            % noise Gen
        linear_BFSKsigPwr = bandpower(linear_BFSKsig);
        linear_BFSKnoisePwr = linear_BFSKsigPwr ./ SNR(i);
        linear_BFSKnoise = sqrt(linear_BFSKnoisePwr/2) .* randn(1,sigLen_encoded);

        hamming_BFSKsigPwr = bandpower(hamming_BFSKsig);
        hamming_BFSKnoisePwr = hamming_BFSKsigPwr ./ SNR(i);
        hamming_BFSKnoise = sqrt(hamming_BFSKnoisePwr/2) .* randn(1,sigLen_encoded);

        cyclic_BFSKsigPwr = bandpower(cyclic_BFSKsig);
        cyclic_BFSKnoisePwr = cyclic_BFSKsigPwr ./ SNR(i);
        cyclic_BFSKnoise = sqrt(cyclic_BFSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        linear_BFSKreceived = linear_BFSKsig + linear_BFSKnoise;

        hamming_BFSKreceived = hamming_BFSKsig + hamming_BFSKnoise;

        cyclic_BFSKreceived = cyclic_BFSKsig + cyclic_BFSKnoise;

            % demodulation
        linear_BFSKdemodLow = linear_BFSKreceived .* (2.* carrierLow_encoded);
        linear_BFSKdemodLow = filtfilt(b, a, linear_BFSKdemodLow);
        linear_BFSKdemodHigh = linear_BFSKreceived .* (2.* carrierHigh_encoded);
        linear_BFSKdemodHigh = filtfilt(b, a, linear_BFSKdemodHigh);
        linear_BFSKdemod = linear_BFSKdemodHigh - linear_BFSKdemodLow;

        hamming_BFSKdemodLow = hamming_BFSKreceived .* (2.* carrierLow_encoded);
        hamming_BFSKdemodLow = filtfilt(b, a, hamming_BFSKdemodLow);
        hamming_BFSKdemodHigh = hamming_BFSKreceived .* (2.* carrierHigh_encoded);
        hamming_BFSKdemodHigh = filtfilt(b, a, hamming_BFSKdemodHigh);
        hamming_BFSKdemod = hamming_BFSKdemodHigh - hamming_BFSKdemodLow;

        cyclic_BFSKdemodLow = cyclic_BFSKreceived .* (2.* carrierLow_encoded);
        cyclic_BFSKdemodLow = filtfilt(b, a, cyclic_BFSKdemodLow);
        cyclic_BFSKdemodHigh = cyclic_BFSKreceived .* (2.* carrierHigh_encoded);
        cyclic_BFSKdemodHigh = filtfilt(b, a, cyclic_BFSKdemodHigh);
        cyclic_BFSKdemod = cyclic_BFSKdemodHigh - cyclic_BFSKdemodLow;

            % sampling and decision logic
        linear_BFSK = samplingANDdecision(linear_BFSKdemod, samplePeriod, encoded_nBits, 0);

        hamming_BFSK = samplingANDdecision(hamming_BFSKdemod, samplePeriod, encoded_nBits, 0);

        cyclic_BFSK = samplingANDdecision(cyclic_BFSKdemod, samplePeriod, encoded_nBits, 0);

             %decoding Code
        BFSK_linear = decode(linear_BFSK,7,4,'linear/binary',genmat);

        BFSK_hamming = decode(hamming_BFSK,7,4,'hamming/binary');

        BFSK_cyclic = decode(cyclic_BFSK,7,4,'cyclic/binary',genpoly,trt);

        %Bit error
        linear_OOKerror =  biterr(OOK_linear, data) ./encoded_nBits;
        OOKerror = biterr(OOK, data) ./nBits;
        linear_BPSKerror = biterr(BPSK_linear, data) ./encoded_nBits;
        linear_BFSKerror = biterr(BFSK_linear, data) ./encoded_nBits;

        linear_aveOOKerror = linear_OOKerror + linear_aveOOKerror;
        aveOOKerror = OOKerror + aveOOKerror ;
        linear_aveBPSKerror = linear_BPSKerror + linear_aveBPSKerror;
        linear_aveBFSKerror = linear_BFSKerror + linear_aveBFSKerror;

        hamming_OOKerror =  biterr(OOK_hamming, data) ./encoded_nBits;
        hamming_BPSKerror = biterr(BPSK_hamming, data) ./encoded_nBits;
        hamming_BFSKerror = biterr(BFSK_hamming, data) ./encoded_nBits;

        hamming_aveOOKerror = hamming_OOKerror + hamming_aveOOKerror;
        hamming_aveBPSKerror = hamming_BPSKerror + hamming_aveBPSKerror;
        hamming_aveBFSKerror = hamming_BFSKerror + hamming_aveBFSKerror;

        cyclic_OOKerror =  biterr(OOK_cyclic, data) ./encoded_nBits;
        cyclic_BPSKerror = biterr(BPSK_cyclic, data) ./encoded_nBits;
        cyclic_BFSKerror = biterr(BFSK_cyclic, data) ./encoded_nBits;

        cyclic_aveOOKerror = cyclic_OOKerror + cyclic_aveOOKerror;
        cyclic_aveBPSKerror = cyclic_BPSKerror + cyclic_aveBPSKerror;
        cyclic_aveBFSKerror = cyclic_BFSKerror + cyclic_aveBFSKerror;

    end     % end for runCycles loop

    OOKerrorArr(i) = aveOOKerror/runCycles;
    linear_OOKerrorArr(i) = linear_aveOOKerror/runCycles;
    linear_BPSKerrorArr(i) = linear_aveBPSKerror/runCycles;
    linear_BFSKerrorArr(i) = linear_aveBFSKerror/runCycles;

    hamming_OOKerrorArr(i) = hamming_aveOOKerror/runCycles;
    hamming_BPSKerrorArr(i) = hamming_aveBPSKerror/runCycles;
    hamming_BFSKerrorArr(i) = hamming_aveBFSKerror/runCycles;

    cyclic_OOKerrorArr(i) = cyclic_aveOOKerror/runCycles;
    cyclic_BPSKerrorArr(i) = cyclic_aveBPSKerror/runCycles;
    cyclic_BFSKerrorArr(i) = cyclic_aveBFSKerror/runCycles;

    % store variable for specified SNR value for plotting
    if (SNRdb(i) == plotSNRdb)
        dataPlot = data;
        cyclicData = cyclic_sig;
        sigPlot = sig;
        linearsigPlot = sig_linear;
        linearmodOOKplot = linear_OOKsig;
        linearrecOOKplot = linear_OOKreceived;
        lineardemodOOKplot = linear_OOK;
        lineardecodeOOKplot = OOK_linear;

        linearmodBPSKplot = linear_BPSKsig;
        linearrecBPSKplot = linear_BPSKreceived;
        lineardemodBPSKplot = linear_BPSK;
        lineardecodeBPSKplot = BPSK_linear;

        linearmodBFSKplot = linear_BFSKsig;
        linearrecBFSKplot = linear_BFSKreceived;
        lineardemodBFSKplot = linear_BFSK;
        lineardecodeBFSKplot = BFSK_linear;

        hammingsigPlot = sig_hamming;
        hammingmodOOKplot = hamming_OOKsig;
        hammingrecOOKplot = hamming_OOKreceived;
        hammingdemodOOKplot = hamming_OOK;
        hammingdecodeOOKplot = OOK_hamming;

        hammingmodBPSKplot = hamming_BPSKsig;
        hammingrecBPSKplot = hamming_BPSKreceived;
        hammingdemodBPSKplot = hamming_BPSK;
        hammingdecodeBPSKplot = BPSK_hamming;

        hammingmodBFSKplot = hamming_BFSKsig;
        hammingrecBFSKplot = hamming_BFSKreceived;
        hammingdemodBFSKplot = hamming_BFSK;
        hammingdecodeBFSKplot = BFSK_hamming;

        cyclicsigPlot = sig_cyclic;
        cyclicmodOOKplot = cyclic_OOKsig;
        cyclicrecOOKplot = cyclic_OOKreceived;
        cyclicdemodOOKplot = cyclic_OOK;
        cyclicdecodeOOKplot = OOK_cyclic;

        cyclicmodBPSKplot = cyclic_BPSKsig;
        cyclicrecBPSKplot = cyclic_BPSKreceived;
        cyclicdemodBPSKplot = cyclic_BPSK;
        cyclicdecodeBPSKplot = BPSK_cyclic;

        cyclicmodBFSKplot = cyclic_BFSKsig;
        cyclicrecBFSKplot = cyclic_BFSKreceived;
        cyclicdemodBFSKplot = cyclic_BFSK;
        cyclicdecodeBFSKplot = BFSK_cyclic;
    end
    
end     % end for SNR loop


% error semilogy plot
figure(1);
p1 = semilogy (SNRdb,linear_OOKerrorArr, 'r-*'); hold on;
p2 = semilogy (SNRdb,OOKerrorArr, 'k-*'); hold;hold;
p3 = semilogy (SNRdb,linear_BPSKerrorArr, 'g-*'); hold;hold;
p4 = semilogy (SNRdb,linear_BFSKerrorArr, 'b-*'); hold;hold
p5 = semilogy (SNRdb,hamming_OOKerrorArr, 'g-*'); hold;hold
p6 = semilogy (SNRdb,hamming_BPSKerrorArr, 'b-*'); hold;hold
p7 = semilogy (SNRdb,hamming_BFSKerrorArr, 'r-*'); hold;hold
p8 = semilogy (SNRdb,cyclic_OOKerrorArr, 'b-*'); hold;hold
p9 = semilogy (SNRdb,cyclic_BPSKerrorArr, 'r-*'); hold;hold;
p10 = semilogy (SNRdb,cyclic_BFSKerrorArr, 'g-*'); hold;
title('Bit Error Rate for different SNR');
legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1) p8(1) p9(1) p10(1)], {'Linear/OOK', 'Unencoded/OOK', 'Linear/BPSK', 'Linear/BFSK','Hamming/OOK', 'Hamming/BPSK', 'Hamming/BFSK', 'Cyclic/OOK', 'Cyclic/BPSK', 'Cyclic/BFSK'}); ylabel('BER'); xlabel('SNR(dB)');
xlim([0 50]);

function [binary] = samplingANDdecision(x, samplePeriod, nBits, threshold)
    sampled = zeros(1, nBits);
    binary = zeros(1, nBits);
    for n = 1 : nBits
        sampled(n) = x((2 * n - 1) * samplePeriod / 2);
        if (sampled(n) > threshold)
            binary(n) = 1;
        else
            binary(n) = 0;
        end
    end
end