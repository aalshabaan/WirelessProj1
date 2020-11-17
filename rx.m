function [rxbits, conf, raw_bits] = rx(rxsignal,conf,k)
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
t = 0:1/conf.f_c:((length(rxsignal)-1)/conf.f_c);
down_converted = rxsignal.*exp(-conf.f_c*2*pi*1i*t');

%Lowpass filtering:

down_converted = lowpass(down_converted, conf);
% Matched-filter
rolloff = 0.22;
pulse = rrc(conf.os_factor, rolloff, conf.mf_length);
filtered_rx = conv(pulse,down_converted);

% Frame sync
preamble = preamble_generate(conf.npreamble);
% Map preamble to BPSK
preamble = 2*preamble - 1;
[data_idx, peak_phase] = frame_sync(filtered_rx,preamble,conf.os_factor);
disp(['FrameSyncIdx ', num2str(data_idx)])

% Time and phase estimation and Interpolation

data_length = floor(length(filtered_rx(data_idx:end))/conf.os_factor);
disp(['nb_symbs: ', num2str(data_length)])
% data_length = conf.nbits;
cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);
data = zeros(1,data_length);
theta_hat = zeros(1,data_length);
theta_hat(1) = peak_phase;
y_array=zeros(data_length,2);

for i=1:data_length
    
     idx_start  = data_idx+(i-1)*conf.os_factor;
     
     idx_range  = idx_start:idx_start+conf.os_factor-1;
     segment    = filtered_rx(idx_range);
    
     % Time error estimatuon
     pwr         = abs(segment).^2;
     fourier_exponents = exp(2*pi*1i*(0:conf.os_factor-1)/conf.os_factor);
     diff_err(i) = fourier_exponents*pwr; % Calculate the power spectrum
     cum_err     = cum_err + diff_err(i);
     epsilon(i)  = -1/(2*pi)*angle(cum_err);
     
     sample_diff   = floor(epsilon(i)*conf.os_factor); % integer
     int_diff      = mod(epsilon(i)*conf.os_factor,1); % interval [0 1)
    
     % linear interpolateion
     y     = filtered_rx(idx_start+sample_diff:idx_start+sample_diff+1);
     y_array(i,:)=y; %debugging purposes
     y_hat = y(1)+int_diff*(y(2)-y(1));
     data(i) = y_hat;
     
     % Phase Estimation
     deltaTheta = 1/4*angle(-data(i)^4) + pi/2*(-1:4);
     [~, ind] = min(abs(deltaTheta - theta_hat(i)));
     theta = deltaTheta(ind);
     theta_hat(i+1) = mod(0.01*theta + 0.99*theta_hat(i), 2*pi);
     
     
     data(i) = data(i) * exp(-1i * theta_hat(i+1));
     
          
end


% Demapping
BPSK_map = [-1 1];
QPSK_map =  1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];

switch conf.modulation_order
    case 1 %BPSK
        disp('BPSK')
        [~,ind] = min((ones(data_length,2)*diag(BPSK_map) - diag(data)*ones(data_length,2)),[],2);
        rxbits = de2bi(ind-1);
        % Unfold into a single column stream
        raw_bits = rxbits;
        rxbits = rxbits(1:conf.nbits);
        
    case 2 %QPSK
        disp('QPSK')
        [~,ind] = min((ones(data_length,4)*diag(QPSK_map) - diag(data)*ones(data_length,4)),[],2);
        rxbits = de2bi(ind-1);
        % Unfold into a single column stream
        raw_bits = rxbits;
        rxbits = rxbits(1:conf.nbits);
    otherwise
        disp('WTF?')
        rxbits = zeros(conf.nbits,1);
end
