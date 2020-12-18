function [rx_signal,conf] = rx_ofdm(rxsignal,conf,k)


%Down conversion
t = 0:1/conf.f_s:((length(rxsignal)-1)/conf.f_s);
down_converted = rxsignal.*exp(-conf.f_c*(1 - conf.offset*1e-6)*2*pi*1i*t');

%Low pass
epsilon_factor=1.1 %increases the cut off freq slightly above BW_bb
cut_of_freq=ceil(((conf.N+1)/2)*conf.f_sep)*(epsilon_factor)
filtered_rx=ofdmlowpass(down_converted,conf,cut_of_freq)


% Frame sync
preamble = preamble_generate(conf.npreamble);
% Map preamble to BPSK
preamble = 2*preamble - 1;
[data_idx, peak_phase] = frame_sync(filtered_rx,preamble,conf.os_factor_sc);
disp(['FrameSyncIdx ', num2str(data_idx)])


%OFDM size: os_factor_ofdm echantillons
%serial : os_factor_ofdm*N 
%filtered_rx[]
%DFT Back to frequency domain

freq_signal = reshape(rx_signal, [], conf.N);
for i = 1:conf.N
   freq_signal(:,i) = osfft(mapped(:,i),conf.os_factor_ofdm); 
end

freq_signal=reshape(freq_signal,[],1);


%freq_signal=osfft(,conf.os_factor_ofdm)


%Phase estimation

rx_signal[








%demapping
BPSK_map = [-1 1];
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];

switch conf.modulation_order
    case BPSK %BPSK
        disp('BPSK')
        [~,ind] = min(abs(ones(data_length,2)*diag(BPSK_map) - diag(data)*ones(data_length,2)),[],2);
        rxbits = de2bi(ind-1);
        % Unfold into a single column stream
        
        rxbits = rxbits(1:conf.nbits);
        
    case QPSK %QPSK
        disp('QPSK')
        [~,ind] = min(abs(ones(data_length,4)*diag(QPSK_map) - diag(data)*ones(data_length,4)),[],2);
        ind = ind+2;
        rxbits = de2bi(ind-1);
        % Unfold into a single column stream
        
        rxbits = rxbits(1:conf.nbits)';
    otherwise
        disp('WTF?')
        rxbits = zeros(conf.nbits,1);



rxbits=rxbits(1:conf.N); %discard the padded 0's

end











