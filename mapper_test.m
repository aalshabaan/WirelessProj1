clc
close all;

data_length = 8;
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
data = QPSK_map;
tx_bits = [0;1;1;0;0;0;0;0;1;1;0;0];
if (mod(size(tx_bits,1),2) ~= 0)
            tx_bits = [tx_bits;0];
        end

bits = 2 * (tx_bits - 0.5);
            bits2 = reshape(bits,2,[]);

            imag = ((bits2(1,:) > 0)-0.5)*sqrt(2);
            real = ((bits2(2,:) > 0)-0.5)*sqrt(2);

            symbol = (real + 1i*imag).';