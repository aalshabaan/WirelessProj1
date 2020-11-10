function [rxbits, conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

% dummy 

% Downconversion
down_converted = rxsignal*e^(conf.f_c*2*pi*1i);

%Lowpass filtering:

down_converted = lowpass(down_converted, conf);
% Matched-filter
rolloff = 0.22;
pulse = rrc(conf.os_factor, rolloff, tx_filterlen);
filtered_rx = conv(pulse,down_converted);

% Frame sync
preamble = preamble_generate(100);
[data_idx, peak_phase] = frame_sync(preamble,conf.os_factor);

% Time and phase estimation and Interpolation

data_length = length(filtered(data_idx:end));
cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);
data = zeros(1,data_length);
theta_hat = zeros(1,data_length);
theta_hat(1) = peak_phase;


for i=1:data_length
    
     idx_start  = data_idx+(i-1)*os_factor;
     
     idx_range  = idx_start:idx_start+os_factor-1;
     segment    = filtered_rx_signal(idx_range);
    
     % Time error estimatuon
     pwr         = abs(segment).^2;
     diff_err(i) = [1 -1j -1 1j]*pwr; % Calculate the power spectrum
     cum_err     = cum_err + diff_err(i);
     epsilon(i)  = -1/(2*pi)*angle(cum_err);
     
     sample_diff   = floor(epsilon(ii)*os_factor); % integer
     int_diff      = mod(epsilon(ii)*os_factor,1); % interval [0 1)
    
     % Phase Estimation
     deltaTheta = 1/4*angle(-payload_data(i)^4) + pi/2*(-1:4);
     [~, ind] = min(abs(deltaTheta - theta_hat(i)));
     theta = deltaTheta(ind);
     theta_hat(i+1) = mod(0.01*theta + 0.99*theta_hat(i), 2*pi);
     
     % linear interpolateion
     y     = filtered_rx(idx_start+sample_diff:idx_start+sample_diff+1);
     y_hat = y(1)+int_diff*(y(2)-y(1));
     data(i) = y_hat;
     
     data(i) = data(i) * exp(-1i * theta_hat(i+1));
     
          
end


% Downsample
downsampled_data = data(1:conf.os_factor:end);
L = length(downsampled_data);

% Demapping
BPSK_map = [-1 1];
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];

switch conf.os_factor
    case 1 %BPSK
        [~,ind] = min((ones(L,2)*diag(BPSK_map) - diag(downsampled_data)*ones(L,2)),[],2);
        rx_bits = de2bi(ind-1);
        
    case 2 %QPSK
        [~,ind] = min((ones(L,4)*diag(QPSK_map) - diag(downsampled_data)*ones(L,4)),[],2);
        rx_bits = de2bi(ind-1);
    otherwise
        disp('WTF?');
        rx_bits = zeros(conf.nbits,1);
end
