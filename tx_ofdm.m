function [tx_signal,conf] = tx_ofdm(tx_bits,conf,k)
%tx_ofdm Transforms transmitted bitstream into an OFDM signal
%   tx_bits: (nbits, 1) The transmitted bitstream, %   conf: The global configuration object
%   k: the current frame index

% mapped is of size ((tx_bits/modulation_order), 1)

mapped = map(tx_bits, conf.modulation_order);

trash_len = mod(conf.N - mod(size(mapped,1),conf.N),conf.N);
trash = zeros(trash_len,1);
mapped = [mapped; trash];
mapped = reshape(mapped, [], conf.N);

% Add training symbol (BPSK OFDM, all -1)
training_sym = -ones(1,conf.N);
mapped = [training_sym ;mapped];
for i = 1:size(mapped,1)
   time_signal(:,i) = osifft(mapped(i,:),conf.os_factor_ofdm); 
end

%Add the cyclic prefix
padding_start_index = floor(size(time_signal,1)*conf.ncp);
cyclic_prefix = time_signal(padding_start_index:end,:);
padded_signal = [cyclic_prefix; time_signal];

padded_signal = reshape(padded_signal,[],1);

% time_signal=reshape(time_signal,[],1);
% 
% %Add cyclic prefix to symbols
% data_length = size(mapped,1);
% %relative index from the start of the OFDM symbol
% padding_start_idx = floor(conf.ncp*data_length*conf.os_factor_ofdm);
% padded_signal = [];
% 
% for i=1:data_length
%     idx_start = 1 + (i-1)*data_length*conf.os_factor_ofdm;
%     idx_range = idx_start:idx_start + data_length*conf.os_factor_ofdm-1;
%     ofdm_symbol = time_signal(idx_range,:);
%     cyclic_prefix = ofdm_symbol(padding_start_idx:end,:);
%     
%     %concatenate the signal with the cyclic prefix and the new symbol
%     padded_signal = [padded_signal;cyclic_prefix;ofdm_symbol];
% end

%Generate a single-carrier preamble, oversample and then pulse shape it
preamble = upsample(map(preamble_generate(conf.npreamble),1), conf.os_factor_sc);
analog_preamble = conv(preamble,rrc(conf.os_factor_sc,0.22,20));

%Normalize the signal's energy
preamble_energy = mean(abs(analog_preamble).^2);
signal_energy = mean(abs(padded_signal).^2);

normalized_signal = padded_signal /(2e4*signal_energy);

%Append the preamble to the padded signal
baseband_signal = [analog_preamble;normalized_signal];

%Upmixing to the carrier frequency
t = 0:1/conf.f_s:(length(baseband_signal)-1)/conf.f_s;
tx_signal = real(baseband_signal.*exp(2*pi*1i*conf.f_c*t'));


end



function [mapped_signal] = map(tx_bits, mapping)

switch mapping
    case 1 %BPSK
        mapped_signal = 1-2*tx_bits;
    case 2 %QPSK
        %if the number of lines (time samples) is odd, pad the end with
        %zeros
        if (mod(size(tx_bits,1),2) ~= 0)
            tx_bits = [tx_bits;0];
        end
        map = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
            sub_signal = reshape(tx_bits,[],2);
            mapped_signal = map(bi2de(sub_signal)+1).';
   
        
    otherwise
        disp('incorrect mapping type');
end

end