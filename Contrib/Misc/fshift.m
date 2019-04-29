function y = fshift(x,s)
% FSHIFT Fractional circular shift
%   Syntax:
%
%       >> y = fshift(x,s)
%
%   FSHIFT circularly shifts the elements of vector x by a (possibly
%   non-integer) number of elements s. FSHIFT works by applying a linear
%   phase in the spectrum domain and is equivalent to CIRCSHIFT for integer
%   values of argument s (to machine precision).

% Author:   Fran�ois Bouffard
%           fbouffard@gmail.com
% Copyright (c) 2015, Francois Bouffard
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
needtr = 0; 
if size(x,1) == 1; 
    x = x(:); 
    needtr = 1; 
end;
N = size(x,1); 
r = floor(N/2)+1; 
f = ((1:N)-r)/(N/2); 
p = exp(-1j*s*pi*f).';
if ~mod(N,2)
    % N is even. This becomes important for complex signals.
    % Thanks to Ahmed Fasih for pointing out the bug.
    % For even signals, f(1) = -1 and phase is sampled at -pi. 
    % The correct value for p(1) should be the average of the f = -1 and
    % f = 1 cases. Since f has antisymmetric phase and unitary (and thus
    % symmetric) magnitude, the average is the real part of p.
    p(1) = real(p(1));
end
y = ifft(fft(x).*ifftshift(p)); 
if isreal(x);
    y = real(y); 
end;
if needtr; 
    y = y.'; 
end;
