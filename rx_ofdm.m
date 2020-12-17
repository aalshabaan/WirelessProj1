function [rx_signal,conf] = rx_ofdm(rx_bits,conf,k)


%Low pass
%f=
%filtered_rx=ofdmlowpass(rx_bits,conf,f)


% Frame sync
preamble = preamble_generate(conf.npreamble);
% Map preamble to BPSK
preamble = 2*preamble - 1;
[data_idx, peak_phase] = frame_sync(filtered_rx,preamble,conf.os_factor);
disp(['FrameSyncIdx ', num2str(data_idx)])


%remove cp

conf.os_factor_ofdm



% Training



%Data


%Back to frequency domain

freq_signal=osfft(,conf.os_factor_ofdm)



end