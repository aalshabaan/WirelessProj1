function [rxbits,conf] = rx_ofdm(rxsignal,conf,k)


%Down conversion
t = 0:1/conf.f_s:((length(rxsignal)-1)/conf.f_s);
down_converted = rxsignal.*exp(-conf.f_c*(1 - conf.offset*1e-6)*2*pi*1i*t');

%Low pass
epsilon_factor=1.1 ;%increases the cut off freq slightly above BW_bb
cut_of_freq=ceil(((conf.N+1)/2)*conf.f_sep)*(epsilon_factor);
filtered_rx=ofdmlowpass(down_converted,conf,cut_of_freq);

%Copy of filtered_rx to, avoids damaging ofdm symbols 
copy_filtered_rx= filtered_rx;
% Matched-filter
rolloff = 0.22;
mf_length=20;
pulse = rrc(conf.os_factor_sc, rolloff, mf_length);
copy_filtered_rx = conv(pulse,copy_filtered_rx);


% Frame sync
preamble = preamble_generate(conf.npreamble);
% Map preamble to BPSK
preamble =  1-2*preamble;
[data_idx, peak_phase] = frame_sync(copy_filtered_rx,preamble,conf.os_factor_sc);
disp(['FrameSyncIdx ', num2str(data_idx)])


%OFDM size: os_factor_ofdm echantillons
%serial : os_factor_ofdm*N 

%We select the filtered signal starting from the begining index given by
%the frame sync
filtered_rx=filtered_rx(data_idx:end);
% post_frame_sync_filtered_rx=filtered_rx;
%DFT Back to frequency domain
 
T=1/conf.f_sep;
len_input_fft=T*conf.f_s;

%trash_len= length(filtered_rx)-mod 

trash_len=len_input_fft-mod(length(filtered_rx),len_input_fft);
trash = zeros(trash_len,1);

filtered_rx = [filtered_rx; trash];

filtered_rx = reshape(filtered_rx,  len_input_fft, []);

num_channels=size(filtered_rx);

for i = 1:num_channels(2)
   freq_signal(:,i) = osfft(filtered_rx(:,i),conf.os_factor_ofdm); 
end

freq_signal=reshape(freq_signal,[],1);


%Equalizer

%Phase estimation


%phase_error=abs(freq_signal(1)-pi)  


%rx_signal[






%demapping
data_length=floor(length(freq_signal)); 
data = freq_signal;

BPSK=1;
QPSK=2;

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


rxbits=rxbits(1:conf.nbits); %discard the padded 0's

end

