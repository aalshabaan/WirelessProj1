function [rxbits, conf] = rx_abed(rxsignal,conf)

%Downconversion

t = 0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s;
baseband_rx = rxsignal .* exp(-2*pi*1i*conf.f_c*t');
% baseband_rx = rxsignal;

%Lowpass filtering

freq_margin = 1.1;
bandwidth = ceil((conf.N+1)/2)*conf.f_sep;
filtered_rx = ofdmlowpass(baseband_rx,conf,freq_margin*bandwidth);

% filtered_rx = baseband_rx;
%Frame synchronisation

preamble = 1 - 2*preamble_generate(conf.npreamble); %Mapped to BPSK
[data_idx, ~] = frame_sync(filtered_rx,preamble,conf.os_factor_sc);
% len_symbol = floor(conf.N*conf.os_factor_ofdm*(1+conf.ncp));
len_symbol = floor(conf.N*conf.os_factor_ofdm + conf.ncp);
relevant_rx = filtered_rx(data_idx:data_idx + conf.nsymbols*len_symbol-1);


%Parallelize the symbols
padded_ofdm_symbols = reshape(relevant_rx, len_symbol, []);
% len_cp = len_symbol * (conf.ncp/(1+conf.ncp));
%Remove the cyclic prefix
% ofdm_time_symbols = padded_ofdm_symbols(len_cp + 1:end,:);
ofdm_time_symbols = padded_ofdm_symbols(conf.ncp + 1:end,:);

% %%%DEBUG%%%
% figure
% plot(abs(ofdm_time_symbols(:,2)))
% hold on
% plot(abs(conf.debug))
% hold off
% %%%DEBUG%%%


%Fourier transfrom
for idx = 1:conf.nsymbols
   ofdm_symbols(idx,:) = osfft(ofdm_time_symbols(:,idx), conf.os_factor_ofdm); 
end

training_sym = ofdm_symbols(1,:);
data_symbs = ofdm_symbols(2:end,:);

%Equalization
H = training_sym/-1;

equalized_data = data_symbs./H;

conf.H = H;

%Serialize data symbols
data = reshape(equalized_data,[],1);



%%

%Demapping
BPSK=1;
QPSK=2;
BPSK_map = [1 -1];
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
data_length = length(data);
switch conf.modulation_order
    case BPSK %BPSK
        disp('BPSK')
        [~,ind] = min(abs(ones(data_length,2)*diag(BPSK_map) - diag(data)*ones(data_length,2)),[],2);
        rxbits = de2bi(ind-1);
        
        
    case QPSK %QPSK
%         disp('QPSK')
%         [~,ind] = min(abs(ones(data_length,4)*diag(QPSK_map) - diag(data)*ones(data_length,4)),[],2);
%         rxbits = de2bi(ind-1)';
%         % Unfold into a single column stream

        LSB = real(data) > 0;
        MSB = imag(data) > 0;
        
        rxbits = [MSB'; LSB'];
         
        rxbits = reshape(rxbits,[],1);
        rxbits = rxbits(1:conf.nbits);
        
        
    otherwise
        disp('WTF?')
        rxbits = zeros(conf.nbits,1);

end


%debug_ber=sum(conf.debug~=rxbits)/length(rxbits)




%%%DEBUG%%%
% figure
% plot(real(data(1:1000)), real(conf.debug_2))
% hold on
% plot(imag(data(1:1000)), imag(conf.debug_2))
% legend('real', 'imag')
% hold off



REAL_ERRS = sum(real(equalized_data) .* real(conf.debug_2) < 0, 'all')
IMAG_ERRS = sum(imag(equalized_data) .* imag(conf.debug_2) < 0, 'all')
%%%DEBUG%%%

end