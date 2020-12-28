% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

% Configuration Values
%<<<<<<< Updated upstream
clear variables

conf.audiosystem = 'bypass';% Values: 'matlab','native','bypass'
%=======
%conf.audiosystem = 'matlab';%'matlab'; % Values: 'matlab','native','bypass'
%>>>>>>> Stashed changes

conf.f_s     = 48000;   % sampling rate  
conf.f_data   = 1000;     % data rate (bps)
conf.nframes = 1;       % number of frames to transmit; this parameter is overridden below
conf.nbits   = 2000;    % number of bits 
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 20000;
conf.N = 250;
conf.f_sep = 8;         %Sub carrier frequency spearation
conf.ncp = 0.5;         %Cyclic prefix length (relative to the symbol length)


conf.npreamble  = 100;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;    % Downconversion frequency offset in ppm, this parameter is overriden below

% Init Section
% all calculations that you only have to do once
conf.os_factor_ofdm  = conf.f_s/(conf.N * conf.f_sep);
% Single carrier os_factor used for the preamble, preamble is BPSK so f_sym
% = f_data
conf.os_factor_sc = conf.f_s/(conf.f_data);
if mod(conf.os_factor_sc,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end
assert(mod(conf.os_factor_ofdm,1) == 0, 'OFDM Oversampling Factor is not an integer')


conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.

offsets = 2:2:10;
f_symbs = [100 200 500 1000 1500 2000];
%%
ber = zeros(length(offsets),length(f_symbs),conf.nframes);
per = zeros(length(offsets),length(f_symbs),1);
% Results
 %for j = 1:length(offsets)
  %   conf.offset = offsets(j);
   %  for i = 1:length(f_symbs)
    %     conf.f_sym = f_symbs(i);
     %    conf.os_factor_sc  = conf.f_s/conf.f_data;
%         
        % if mod(conf.os_factor,1) ~= 0
         %    disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
         %end
        for k=1:conf.nframes

            % Generate random data
            txbits = randi([0 1],conf.nbits,1);

            % TODO: Implement tx() Transmit Function
            [txsignal, conf] = tx_ofdm(txbits,conf,k);

            % % % % % % % % % % % %
            % Begin
            % Audio Transmission
            %

            % normalize values
            peakvalue       = max(abs(txsignal));
            normtxsignal    = txsignal / (peakvalue + 0.3);
            
            
            % create vector for transmission
            rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
            rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
            txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal

        %     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
            audiowrite('out.wav',rawtxsignal,conf.f_s)  

            % Platform native audio mode 
            if strcmp(conf.audiosystem,'native')

                % Windows WAV mode 
                if ispc()
                    disp('Windows WAV');
                    wavplay(rawtxsignal,conf.f_s,'async');
                    disp('Recording in Progress');
                    rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
                    disp('Recording complete')
                    rxsignal = rawrxsignal(1:end,1);

                % ALSA WAV mode 
                elseif isunix()
                    disp('Linux ALSA');
                    cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
                    system(cmd); 
                    disp('Recording in Progress');
                    system('aplay  out.wav')
                    pause(2);
                    disp('Recording complete')
                    rawrxsignal = wavread('in.wav');
                    rxsignal    = rawrxsignal(1:end,1);
                end

            % MATLAB audio mode
            elseif strcmp(conf.audiosystem,'matlab')
                disp('MATLAB generic');
                playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
                recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
                record(recobj);
                disp('Recording in Progress');
                playblocking(playobj)
                pause(0.5);
                stop(recobj);
                disp('Recording complete')
                rawrxsignal  = getaudiodata(recobj,'int16');
                rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;

            elseif strcmp(conf.audiosystem,'bypass')
                rawrxsignal = rawtxsignal(:,1);
                rxsignal    = rawrxsignal;
            end

            % Plot received signal for debgging
%             figure;
%             plot(rxsignal);
%             title('Received Signal')

            %
            % End
            % Audio Transmission   
            % % % % % % % % % % % %

            % TODO: Implement rx() Receive Function
            [rxbits conf]       = rx_abed(rxsignal,conf);

            res.rxnbits(k)      = length(rxbits);  
            res.biterrors(k)    = sum(rxbits ~= txbits);

            ber(1,1,k)=(res.biterrors(k))/(res.rxnbits(k));
            
        %    ber(j,i,k) = (res.biterrors(k))/(res.rxnbits(k));
       %  end
      %   per(j,i) = sum(ber(j,i,:) > 0)/conf.nframes;
     %end
 end


%%
figure
for i = 1:length(offsets);
    hold on
    plot(f_symbs,mean(ber(i,:,:),3),'DisplayName', ['offset = ',num2str(offsets(i)), ' parts per million']);
end
hold off
xlabel('Symbol Rate (Bd)')
ylabel('BER')
legend('show')
title('Bit Error Rate as a function of Symbol Frequency')
