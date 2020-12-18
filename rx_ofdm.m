function [rx_signal,conf] = rx_ofdm(rx_bits,conf,k)

%Low pass
%cut_of_freq=
%filtered_rx=ofdmlowpass(rx_bits,conf,cut_of_freq)


% Frame sync
preamble = preamble_generate(conf.npreamble);
% Map preamble to BPSK
preamble = 2*preamble - 1;
[data_idx, peak_phase] = frame_sync(filtered_rx,preamble,conf.os_factor);
disp(['FrameSyncIdx ', num2str(data_idx)])


%remove cp


%cp_length=floor(conf.ncp*conf.os_factor_ofdm)
%remove_cp(filtered_rx, conf.os_factor_ofdm, cp_length)



% Training



%Data


%Back to frequency domain

%freq_signal=osfft(,conf.os_factor_ofdm)

%Equalizer

%a=inv(osfft(h))*freq_signal


%demapping

end



function output=remove_cp(input, len_symbol, cp_length)
%remove_cp: removes cyclic prefix from recieved signal
%input: input vector to remove the cp from
%len_symbol: length of a symbol (without cp)
%cp_length: length of cp

%output: cleaned input vector from any cp

    num_symbols=floor(size(input,1)/(cp_length+len_symbol));
    output=[];
    
        for i=1:num_symbols
            output(:,i)=input((i-1)*(len_symbol+cp_length)+cp_length+1:(i-1)*(len_symbol+cp_length)+cp_length+len_symbol);

        end
        
    output=reshape(output,[],1);
end









