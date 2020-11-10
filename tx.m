function [txsignal conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

% dummy 400Hz sinus generation
%time = 1:1/conf.f_s:4;
%txsignal = 0.3*sin(2*pi*400 * time.');


%Begining of processing


%Preamble



signal=[preamble_generate(conf.npreamble);txbits];


%Mapping
%BPSK FOR PREAMBLE
BPSK=1;
signal(1:conf.npreamble)= mapping(signal(1:conf.npreamble),BPSK);

signal(conf.npreamble+1:length(signal))= mapping(signal(conf.npreamble+1:length(signal)),conf.modulation_order);

%Upsampling

%signal=[zeros(1,conf.os_factor);signal];

signal= upsample(signal,conf.os_factor);

%Matched filter RRC
MF_LENGTH=20;
baseband_signal=matched_filter(signal, conf.os_factor, MF_LENGTH);

%Up conversion
txsignal=up_conversion(baseband_signal,conf.f_c);

end


%Here are the prototypes of all the functions used for the filter

%Do not forget to normalize the energy!



function out_signal=mapping(input_signal,mapping_type)

BPSK=1;
QPSK=2;

    if(mapping_type==BPSK)
        out_signal=2*input_signal-1;
   
    elseif(mapping_type==QPSK)
        
        NonGrayMap = 1/sqrt(2) * [( 1-1j) ( 1+1j) (-1+1j) (-1-1j)];
        out_signal = NonGrayMap(bi2de(input_signal)+1).' ;
        
    else
        disp('Incorrect mapping type')
    end
            
end




function filtered_signal = matched_filter(signal, os_factor, mf_length)

    rolloff_factor = 0.22;

    h = rrc(os_factor, rolloff_factor, mf_length);
    filtered_signal = conv(h, signal);

end

function passband_signal=up_conversion(baseband_signal, fc)

    passband_signal= real(baseband_signal.*exp(i*2*pi*fc));

end

function h = rrc(os_factor, rolloff_factor, filterlength)
% Returns a the FIR coefficients of a Root Raised Cosine filter.
% os_factor is the oversampling factor, typically set to 4.
% rolloff_factor is a in [0; 1]. UMTS for instance uses a rolloff factor of 0.22.
% filterlength is the _onesided_ filterlength, i.e. the total number of taps is 2*filterlength+1.
%
% Note that the current implementation of this function does not handle the case when the denominator becomes zero,
% which means that the rolloff_factor must be chosen so that os_factor/(4*rolloff_factor) is not be an integer.
% In this case, choose a slightly different rolloff factor, or fix this implementation...

n = (-filterlength : filterlength)' / os_factor;

if rolloff_factor == 0, % special case, just return a sinc
	h = sinc(n);
else
	if floor(os_factor/(4*rolloff_factor)) == os_factor/(4*rolloff_factor), error('os_factor/(4*rolloff_factor) must not be an integer.'); end
	h = (4*rolloff_factor/pi * cos((1+rolloff_factor)*pi*n) + (1-rolloff_factor)*sinc((1-rolloff_factor)*n)) ./ (1 - (4*rolloff_factor*n).^2);
end

% Normalize to a total power of 1.
h = h / sqrt(sum(abs(h).^2));
end

