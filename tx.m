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
conf.mf_length=20;
baseband_signal=matched_filter(signal, conf.os_factor, conf.mf_length);

%Up conversion
txsignal=up_conversion(baseband_signal,conf.f_c);

end


%Here are the prototypes of all the functions used for the filter


function out_signal=mapping(input_signal,mapping_type)

BPSK=1;
QPSK=2;

    if(mapping_type==BPSK)
        out_signal=2*input_signal-1;
   
    elseif(mapping_type==QPSK)
        
        reshape(input_signal, [], 2);
        GrayMap = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
        out_signal = GrayMap(bi2de(input_signal)+1).' ; % This doesn't work cuz signal is a single column here
        
%         size_a=length(temp_out_signal)/2;
%         out_signal=reshape(temp_out_signal,[size_a,2]);
%         size(out_signal);
        
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
    
    t = 0:1/fc:((length(baseband_signal)- 1)/fc);
    passband_signal= real(baseband_signal.*exp(1i*2*pi*fc*t'));

end


