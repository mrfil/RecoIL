function result_f=k2image_we(uu,vv,k,ww,N,we,L,tt,W)
% function result_f=k2image_we(uu,vv,k,ww,N,we,L,tt,W)
% this function gets gridded k data from non cartesian k data
% then returns the NxN ifft image of it. 
% uu,vv is coordinates of the noncartesian k data
% ww is the data density precompensation weighting
%	it is obtained by weight_vor, with some modifications on 
%	its behavior at the edges of k-space.
% we is the fieldmap, in rad/sec.
% L is the number of time points to segment
% tt is the time vector for the k-space trajectory
% tau is the delay between echo times
% Beware the half-pixel shift, need to change flag in file if you
% do not want this.
% Maintainers:  Brad Sutton, Sangwoo Lee
% may 1, 2003  
  
%flag for shifting by half a pixel
flag_shft = 0;

if (nargin<9)
  W=6; % kernel width
end

%disp('new');
%xi=(-(N/2-0.25)):0.5:(N/2-0.25);         % doubles the matrix size and will be
%yi=xi;
%yi=(-(N/2-0.25)):0.5:(N/2-0.25);                       % subsampled after gridding
xi=(-N/2):0.5:(N/2-0.5);         % doubles the matrix size and will be
%yi=xi;
yi=(-(N/2)):0.5:(N/2-0.5);                       % subsampled after griddin


%Initial range check
if ((max(uu)>N/2) | (max(vv)>N/2) | (min(uu)<(-N/2)) | (min(vv)<(-N/2)))
     sprintf('k-space trajectory is out of range, functionality for k-space deviations outside of -N/2 to N/2 has not been taken into account')
end
    
% fmygrid2 : uses loc and output driven convolution
%loc=getloc(k,N,5,5);
% [temp lx ly bx by]=fmygrid2(uu,vv,k,xi,yi,2.5,2.5,'KB',loc,N);

% mygrid2 : uses input driven convolution
% [temp lx ly bx by]=mygrid2(uu,vv,k,xi,yi,2.5,2.5,'KB');

% ffmygrid : uses output driven convolution
% [temp lx ly bx by]=ffmygrid(uu,vv,k,xi,yi,2.5,2.5,'KB');

% fmygrid3 : uses loc and output driven convolution with weight correction
%loc=getloc(k,N,5,5);
%[temp lx ly bx by]=fmygrid3(uu,vv,k,xi,yi,2.5,2.5,'KB',loc,N);

% ffmygrid4 : uses output driven convolution with weight correction
%[temp lx ly bx by]=ffmygrid(uu,vv,ww.*k,xi,yi,W,W,'KB');
%Post gridding weight here.
post_wt = real(ffmygrid(uu,vv,ww.*ones(length(k),1),xi,yi,W,W,'KB'));


we = reshape(we,N,N);
min_tt = min(tt);
if ~(L==1)
  tau = (max(tt)-min_tt)/(L-1);
else
  tau = 0;
end 
          % to be deleted afterwhile
x=linspace(-0.5,0.5,N*2+1);
x=x(1:2*N);
%y=linspace(-0.5,0.5,N*2+1);
y=linspace(-0.5,0.5,N*2+1);
y=y(1:2*N);
if (flag_shft == 1)
  y = y+1/(4*N);
  x = x+1/(4*N);
end

for ll = 1:L
  Wo = exp(i*we*((ll-1)*tau+min_tt));
  aa = zeros(size(tt,1),1);
  jj = find(abs(tt-((ll-1)*tau+min_tt))<=tau);
  aa(jj) = 1/2 + 1/2*cos(pi*(tt(jj)-((ll-1)*tau+min_tt))/tau);
  [temp lx ly bx by] = ffmygrid(uu,vv,aa.*ww.*k,xi,yi,W,W,'KB');
  if ll == 1
	lx = 2*lx;   %THESE LINES HERE??
	ly = 2*ly;
	cx=sin(sqrt(pi^2*lx^2*x.^2-bx^2))./sqrt(pi^2*lx^2*x.^2-bx^2);
	cy=sin(sqrt(pi^2*ly^2*y.^2-by^2))./sqrt(pi^2*ly^2*y.^2-by^2);
	[ccx ccy]=meshgrid(cx,cy);
	c=ccy.*ccx;
	c=c/max(max(c));
  end
  l = find(~(post_wt==0));
  msk = logical(masking(ones(2*N),N));
  %temp(l) = temp(l)./mean(post_wt(msk));
  temp = temp./mean((post_wt(msk)));
  %l = find(isnan(temp));
  %temp(l)=0;

  if flag_shft
    shft = 1/2;
    phs = exp(i*2*pi*((yi./N)'*shft))*exp(i*2*pi*((xi./N)*shft));
    temp = temp.*phs; 
  end

  result=fftshift(ifft2(fftshift(temp)))./c;

	% to be deleted afterwhile
  result=result([N/2+1:3/2*N],[N/2+1:3/2*N]);

  if ll == 1
     result_f = Wo.*result;
  else
     result_f = result_f + Wo.*result;
  end            
end











