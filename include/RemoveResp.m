function [ corrected_data ] = RemoveResp(dataS,z2skip,wl)

% [ corrected_data ] = RemoveResp(dataS,z2skip,wl)
% 
% calculates complex response from input poles and zeros
% and removes response from data vector
% dataS is data structure from irisFetch using the includePZ option
% z2skip = 0 for displacement (default),
%          1 for velocity,
%          2 for acceleration,
% wl is waterlevel to use (default is .001)
debug=0;

nnfft=2^nextpow2(length(dataS.data)*4);
nnfft=ceil(nnfft/2)*2;
 freqs=dataS.sampleRate.*[0:nnfft/2 -nnfft/2+1:-1]'/(nnfft*1); 
 
omegas=2*pi*freqs;

if nargin < 2,
    z2skip = 0;
    wl=.001;
end
if nargin < 3,
    wl=.001;
end

hsts2s=ones(size(omegas));
for n = (1+z2skip):length(dataS.sacpz.zeros),
    hsts2s=hsts2s.*(1i.*omegas - dataS.sacpz.zeros(n));
end
for n = 1:length(dataS.sacpz.poles),
    hsts2s=hsts2s./(1i.*omegas - dataS.sacpz.poles(n));
end
complex_response = hsts2s.*dataS.sacpz.constant;
[yy,ii]=min(abs(freqs-dataS.sensitivityFrequency));

minwl=(abs(complex_response(ii)))*wl;
ii=find(abs(complex_response));
minwl=max([ min(abs(complex_response(ii))) minwl]);

complex_response2 = (abs(complex_response) + minwl) .*exp(1i*unwrap(angle(complex_response)));

if debug==1,
    figure(9); clf
    subplot(2,1,1)
    loglog(freqs,abs(complex_response2))
    subplot(2,1,2)
    semilogx(freqs,(angle(complex_response2)))
    pause
end

datafft1=fft(detrend(dataS.data),nnfft);

 datafft2=datafft1.*(conj(complex_response2))./(conj(complex_response2).*complex_response2);
 datafilt1=real(ifft(datafft2));
 corrected_data=datafilt1(1:length(dataS.data));

return;

