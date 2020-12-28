disp('QPSK')

data_length = 8;
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
data = QPSK_map;
tx_bits = [0;0;0;1;1;0;1;1;0;1;0;1;1;0;1;0];
if (mod(size(tx_bits,1),2) ~= 0)
            tx_bits = [tx_bits;0];
        end

% bits = 2 * (tx_bits - 0.5);
% bits2 = reshape(bits, [], 2);
% 
% real = ((bits2(:,2) > 0)-0.5)*sqrt(2);
% imag = ((bits2(:,1) > 0)-0.5)*sqrt(2);
% 
% symbol = (real + 1i*imag);

bits = 2 * (tx_bits - 0.5);
            bits2 = reshape(bits, [], 2);

            real = ((bits2(:,2) > 0)-0.5)*sqrt(2);
            imag = ((bits2(:,1) > 0)-0.5)*sqrt(2);

            symbol = (real + 1i*imag);


disp('QPSK')
[~,ind] = min(abs(ones(data_length,4)*diag(QPSK_map) - diag(symbol)*ones(data_length,4)),[],2);
%         ind = ind+2;
    foo= de2bi(ind-1);
        % Unfold into a single column stream
        
    
bar = reshape(foo,[],1);

ber = mean(bar ~= tx_bits)


