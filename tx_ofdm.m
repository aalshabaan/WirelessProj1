function [tx_signal,conf] = tx_ofdm(tx_bits,conf)
%tx_ofdm Transforms transmitted bitstream into an OFDM signal
%   tx_bits: (nbits, N) The transmitted bitstream, one subcarrier per
%   column
%   conf: The global configuration object



end



function [mapped_signal] = map(tx_bits, mapping)

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