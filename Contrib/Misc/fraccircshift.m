function B = fraccircshift(A,shiftsize)
%fraccircshift expands circshift to fractional shifts values, using linear
%interpolation. In contrast to other approaches to non-integer shifts 
%of matrices on the base of fft2 or interp2, the number of dimensions of A 
%is not limited, and the syntax of circshift applies. For integer elements 
%of shiftsize, fracircshift and circshift give the same results. 
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
int = floor(shiftsize);     %integer portions of shiftsize
fra = shiftsize - int;      %fractional portions of shiftsize
dim = numel(shiftsize);
B = A;
for n = 1:numel(shiftsize)  %The dimensions are treated one after another.
    intn = int(n);
    fran = fra(n);
    shift1 = zeros(dim,1);
    shift1(n) = intn;
    shift2 = zeros(dim,1);
    shift2(n) = intn+1;
    %Linear intepolation:
    B = (1-fran)*circshift(B,shift1) + fran*circshift(B,shift2);
end
