       %***** On-Off Keying (OOK) *****
           % modulation
        encoded_OOKsig = sig_encoded .* carrier_encoded;                        % if 1 = carrier; if 0 = 0;
        OOKsig = sig .* carrier;

            % noise Gen
        OOKsigPwr = bandpower(OOKsig);                  % sum of absolute square/length
        encoded_OOKsigPwr = bandpower(encoded_OOKsig);

        OOKnoisePwr = OOKsigPwr ./ SNR(i);
        OOKnoise = sqrt(OOKnoisePwr/2) .* randn(1,sigLen); %%%% align with Phase1
        encoded_OOKnoisePwr = encoded_OOKsigPwr ./ SNR(i);
        encoded_OOKnoise = sqrt(encoded_OOKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        OOKreceived = OOKsig + OOKnoise;
        encoded_OOKreceived = encoded_OOKsig + encoded_OOKnoise;

            % demodulation
        OOKdemod = OOKreceived .* (2.* carrier);        %coherent
        OOKdemod = filtfilt(b, a, OOKdemod);
        encoded_OOKdemod = encoded_OOKreceived .* (2.* carrier_encoded);        %coherent
        encoded_OOKdemod = filtfilt(b, a, encoded_OOKdemod);

            % sampling and decision logic
        OOK = samplingANDdecision(OOKdemod, samplePeriod, nBits, amp/2);
        encoded_OOK = samplingANDdecision(encoded_OOKdemod, samplePeriod, encoded_nBits, amp/2);

            %decoding Cyclic Code
        OOK_decoded = decode(encoded_OOK,7,4,'cyclic/binary');



        %***** Binary Phase Shift Keying (BPSK) *****
            % modulation
        encoded_BPSKsig = sig_encoded .* 2 - 1;                         % if 1 -> 1*2-1 = 1; if 0 -> 0*2-1 = -1
        encoded_BPSKsig = encoded_BPSKsig .* carrier_encoded;

            % noise Gen
        encoded_BPSKsigPwr = bandpower(encoded_BPSKsig);
        encoded_BPSKnoisePwr = encoded_BPSKsigPwr ./ SNR(i);
        encoded_BPSKnoise = sqrt(encoded_BPSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        encoded_BPSKreceived = encoded_BPSKsig + encoded_BPSKnoise;

            % demodulation
        encoded_BPSKdemod = encoded_BPSKreceived .* (2.* carrier_encoded);
        encoded_BPSKdemod = filtfilt(b, a, encoded_BPSKdemod);

            % sampling and decision logic
        encoded_BPSK = samplingANDdecision(encoded_BPSKdemod, samplePeriod, encoded_nBits, 0);

            %decoding Cyclic Code
        BPSK_decoded = decode(encoded_BPSK,7,4,'cyclic/binary');

        %***** Binary Frequency Shift Keying *****
            % modulation
        encoded_BFSKsigLow = carrierLow_encoded .* (sig_encoded == 0);
        encoded_BFSKsigHigh = carrierHigh_encoded .* (sig_encoded == 1);
        encoded_BFSKsig = encoded_BFSKsigLow + encoded_BFSKsigHigh;

            % noise Gen
        encoded_BFSKsigPwr = bandpower(encoded_BFSKsig);
        encoded_BFSKnoisePwr = encoded_BFSKsigPwr ./ SNR(i);
        encoded_BFSKnoise = sqrt(encoded_BFSKnoisePwr/2) .* randn(1,sigLen_encoded);

            % transmission
        encoded_BFSKreceived = encoded_BFSKsig + encoded_BFSKnoise;

            % demodulation
        encoded_BFSKdemodLow = encoded_BFSKreceived .* (2.* carrierLow_encoded);
        encoded_BFSKdemodLow = filtfilt(b, a, encoded_BFSKdemodLow);
        encoded_BFSKdemodHigh = encoded_BFSKreceived .* (2.* carrierHigh_encoded);
        encoded_BFSKdemodHigh = filtfilt(b, a, encoded_BFSKdemodHigh);
        encoded_BFSKdemod = encoded_BFSKdemodHigh - encoded_BFSKdemodLow;

            % sampling and decision logic
        encoded_BFSK = samplingANDdecision(encoded_BFSKdemod, samplePeriod, encoded_nBits, 0);

             %decoding Cyclic Code
        BFSK_decoded = decode(encoded_BFSK,7,4,'cyclic/binary');

        %Bit error
        encoded_OOKerror =  biterr(OOK_decoded, data) ./encoded_nBits;
        OOKerror = biterr(OOK, data) ./nBits;
        encoded_BPSKerror = biterr(BPSK_decoded, data) ./encoded_nBits;
        encoded_BFSKerror = biterr(BFSK_decoded, data) ./encoded_nBits;

        encoded_aveOOKerror = encoded_OOKerror + encoded_aveOOKerror;
        aveOOKerror = OOKerror + aveOOKerror ;
        encoded_aveBPSKerror = encoded_BPSKerror + encoded_aveBPSKerror;
        encoded_aveBFSKerror = encoded_BFSKerror + encoded_aveBFSKerror;

    end     % end for runCycles loop

    OOKerrorArr(i) = aveOOKerror/runCycles;
    encoded_OOKerrorArr(i) = encoded_aveOOKerror/runCycles;
    encoded_BPSKerrorArr(i) = encoded_aveBPSKerror/runCycles;
    encoded_BFSKerrorArr(i) = encoded_aveBFSKerror/runCycles;

    % store variable for specified SNR value for plotting
    if (SNRdb(i) == plotSNRdb)
        dataPlot = data;
        cyclicData = cyclic_sig;
        sigPlot = sig;
        cyclicsigPlot = sig_encoded;
        modOOKplot = encoded_OOKsig;
        recOOKplot = encoded_OOKreceived;
        demodOOKplot = encoded_OOK;
        decodeOOKplot = OOK_decoded;

        modBPSKplot = encoded_BPSKsig;
        recBPSKplot = encoded_BPSKreceived;
        demodBPSKplot = encoded_BPSK;
        decodeBPSKplot = BPSK_decoded;

        modBFSKplot = encoded_BFSKsig;
        recBFSKplot = encoded_BFSKreceived;
        demodBFSKplot = encoded_BFSK;
        decodeBFSKplot = BFSK_decoded;
    end