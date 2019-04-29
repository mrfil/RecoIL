function [result,tmp1,c] = k2image(uu,vv,k,ww,N,W)
% function [result,tmp1,c] = k2image(uu,vv,k,ww,N,W)
% this function gets gridded k data from non cartesian k data
% then returns the NxN ifft image of it. 
% uu,vv is coordinates of the noncartesian k data
% ww is the data density precompensation weighting
%	it is obtained by weight_vor, with some modifications on 
%	its behavior at the edges of k-space.
% Beware the half-pixel shift, need to change flag in file if you
% do not want this.
% Maintainers:  Brad Sutton, Sangwoo Lee
% May 1, 2003
  
%flag for shifting by half a pixel
flag_shft = 0;

if (nargin<6)
   W=6; % kernel width
end


%disp('new');
%xi=(-(N/2-0.25)):0.5:(N/2-0.25);         % doubles the matrix size and will be
%yi=xi;
%yi=(-(N/2-0.25)):0.5:(N/2-0.25);                       % subsampled after gridding
xi=(-N/2):0.5:(N/2-0.5);         % doubles the matrix size and will be
%yi=xi;
yi=(-(N/2)):0.5:(N/2-0.5);                       % subsampled after gridding

%Initial range check
if ((max(uu)>N/2) || (max(vv)>N/2) || (min(uu)<(-N/2)) || (min(vv)<(-N/2)))
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
[temp, lx, ly, bx, by]=ffmygrid(uu,vv,ww.*k,xi,yi,W,W,'KB');
%Post gridding weight here.
if 1
post_wt = real(ffmygrid(uu,vv,ww.*ones(length(k),1),xi,yi,W,W,'KB'));
%l = find(~(post_wt==0));
msk = logical(masking(ones(2*N),N));
temp = temp./mean((post_wt(msk)));
%temp(l) = temp(l)./(post_wt(l));
end

l = isnan(temp);
temp(l)=0;


if (nargout > 1)
    tmp1 = temp;
end


           % to be deleted afterwhile
x=linspace(-0.5,0.5,N*2+1);
x=x(1:2*N); %+1/(2*N);
%y=linspace(-0.5,0.5,N*2+1);
y=linspace(-0.5,0.5,N*2+1);
y=y(1:2*N); %+1/(2*N);
%if (flag_shft == 1)
%  y = y+1/(4*N);
%  x = x+1/(4*N);
%end
lx = 2*lx;   %THESE LINES HERE??
ly = 2*ly;
cx=sin(sqrt(pi^2*lx^2*x.^2-bx^2))./sqrt(pi^2*lx^2*x.^2-bx^2);
cy=sin(sqrt(pi^2*ly^2*y.^2-by^2))./sqrt(pi^2*ly^2*y.^2-by^2);
[ccx, ccy]=meshgrid(cx,cy);
c=ccy.*ccx;
c=c/max(max(c));

if flag_shft
  shft = 1/2;
  phs = exp(1i*2*pi*((yi./N)'*shft))*exp(1i*2*pi*((xi./N)*shft));
  temp = temp.*phs; 
end

result=fftshift(ifft2(fftshift(temp)))./c;
%keyboard

% to be deleted afterwhile
% Anh: floor is added to take care of the case when N is odd
result=result([floor(N/2)+1:floor(3/2*N)],[floor(N/2)+1:floor(3/2*N)]);







