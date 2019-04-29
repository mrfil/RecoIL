function AA = int_tim_seghanning(tt,L)
%function AA = int_tim_seg(tt,L)
%  to be used with function fast_mr.m in @fast_mr
%  tt is time vector
%  L is number of segments, i.e. number of interpolation points
%  This version just uses a Hanning window
%    For use of fewer time segments choose the optimal time segmented interpolator in
%     int_tim_seg.m == From Fastmr paper, Sutton, 2003.


mint = min(tt(:));
maxt = max(tt(:));
rangt = maxt-mint;
ndat = length(tt);

tau = (rangt+eps)/(L);

AA = zeros(L+1,ndat);

if L==0
  AA = ones(1,ndat); 
  return
end

tt = tt(:)-mint;


for ll = 1:L+1
   argument_inner = tt - (ll-1)*tau;
   AA(ll,:) = (0.5 + 0.5.*cos(pi*argument_inner(:)'./tau)).*(abs(argument_inner(:)')<tau);
end


