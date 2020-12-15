function [tx_signal,conf] = tx_ofdm(tx_bits,conf)
%tx_ofdm Transforms transmitted bitstream into an OFDM signal
%   tx_bits: (nbits, N) The transmitted bitstream, one subcarrier per
%   column
%   conf: The global configuration object

% mapped is of shape ((tx_bits/modulation_order)+1, N)
%The +1 is the training OFDM symbol used for phase estimation;
mapped = map(tx_bits, conf.modulation_order);

% time_signal is of shape(os_factor*((tx_bits/modulation_order)+1),N)
time_signal = osifft(mapped, conf.os_factor);

%Add cyclic prefix to symbols
data_length = shape(tx_bits,1);
padding_start_idx = floor(conf.ncp*conf.os_factor);
padded_signal = [];

for i=1:data_length
    idx_start = 1 + (i-1)*conf.os_factor;
    idx_range = idx_start:idx_start + conf.os_factor;
    ofdm_symbol = time_signal(idx_range);
    cyclic_prefix = ofdm_symbol(padding_start_idx:end);
    padded_signal = [padded_signal;cyclic_prefix;ofdm_symbol];
end

%Generate a single-carrier preamble, oversample and then pulse shape it
preamble = upsample(map(preamble_generate(conf.npreamble),1), conf.os_factor);
analog_preamble = conv(preamble,rrc(conf.os_factor,0.22,20));

%Assign the preamble to the first carrier
ofdm_preamble = zeros(shape(analog_preamble,1),conf.N);
ofdm_preamble(:,1) = analog_preamble;

%Append the preamble to the padded signal
baseband_signal = [ofdm_preamble;padded_signa];

%Upmixing to the carrier frequency
t = 0:1/conf.f_s:(length(baseband_signal)-1)/conf.f_s;
tx_signal = baseband_signal.*exp(2*pi*1i*conf.f_c*t');


end



function [mapped_signal] = map(tx_bits, mapping)
tx_bits = [ones(1,shape(tx_bits,2)); tx_bits];
switch mapping
    case 1 %BPSK
        mapped_signal = 1-2*tx_bits;
    case 2 %QPSK
        mapped_signal = zeros(ceil(shape(tx_bits,1)/2),shape(tx_bits,2));
        map = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
        for i = 1:shape(tx_bits,2)
            sub_signal = reshape(tx_bits(:,i),[],2);
            mapped_signal(:,i) = map(bi2de(sub_signal)+1).';
        end
        
    otherwise
        disp('incorrect mapping type');
end




end