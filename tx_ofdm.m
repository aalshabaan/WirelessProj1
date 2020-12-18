function [tx_signal,conf] = tx_ofdm(tx_bits,conf,k)
%tx_ofdm Transforms transmitted bitstream into an OFDM signal
%   tx_bits: (nbits, 1) The transmitted bitstream, %   conf: The global configuration object
%   k: the current frame index

% mapped is of size ((tx_bits/modulation_order)+1, 1)
%The +1 is the training OFDM symbol used for phase estimation;
mapped = map(tx_bits, conf.modulation_order);


trash_len = conf.N - mod(size(mapped,1),conf.N);
trash = zeros(trash_len,1);
mapped = [mapped; trash];
mapped = reshape(mapped, [], conf.N);
for i = 1:conf.N
   time_signal(:,i) = osifft(mapped(:,i),conf.os_factor_ofdm); 
end

time_signal=reshape(time_signal,[],1);

%Add cyclic prefix to symbols
data_length = size(mapped,1);
%relative index from the start of the OFDM symbol
padding_start_idx = floor(conf.ncp*conf.os_factor_ofdm);
padded_signal = [];

for i=1:data_length
    idx_start = 1 + (i-1)*conf.os_factor_ofdm;
    idx_range = idx_start:idx_start + conf.os_factor_ofdm-1;
    ofdm_symbol = time_signal(idx_range,:);
    cyclic_prefix = ofdm_symbol(padding_start_idx:end,:);
    
    %concatenate the signal with the cyclic prefix and the new symbol
    padded_signal = [padded_signal;cyclic_prefix;ofdm_symbol];
end

%Generate a single-carrier preamble, oversample and then pulse shape it
preamble = upsample(map(preamble_generate(conf.npreamble),1), conf.os_factor_sc);
analog_preamble = conv(preamble,rrc(conf.os_factor_sc,0.22,20));

%Assign the preamble to the first carrier
ofdm_preamble = zeros(size(analog_preamble,1),size(padded_signal,2));
ofdm_preamble(:,1) = analog_preamble;

%Append the preamble to the padded signal
baseband_signal = [ofdm_preamble;padded_signal];

%Upmixing to the carrier frequency
t = 0:1/conf.f_s:(length(baseband_signal)-1)/conf.f_s;
tx_signal = real(baseband_signal.*exp(2*pi*1i*conf.f_c*t'));


end



function [mapped_signal] = map(tx_bits, mapping)
% A BPSK-encoded 1
training_symbol = -1;

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

mapped_signal = [training_symbol; mapped_signal];


end