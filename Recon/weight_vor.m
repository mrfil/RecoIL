function area_out=weight_vor(kx,ky,nl,skip_fix)
% function area=weight_vor(kx,ky,nl,skip_fix)
% nl is the number of interleaves
% For calculating sampling density function
%   for spiral trajectory
%  Sangwoo Lee and Brad Sutton
%   University of Michigan

area_out = zeros(1,length(kx));
[B,I,J] = unique([kx,ky],'rows','first');

[V,C]=voronoin(B);

if ~exist('skip_fix','var')
    skip_fix = 0;
end

len_used = length(B(:,1));
%len = length(kx);

area=zeros(1,len_used);

%keyboard
for ii=1:len_used
    xx=[];yy=[];s=[];
    %keyboard
    s=C{ii}(1:end);
    xx=V(s,1);
    yy=V(s,2);
    area(ii)=polyarea(xx,yy);
end;

%keyboard
area_out = area(J);
area = area_out;

%Now fix edges
if skip_fix
else
  for jj = 1:nl

  end_j = jj*(length(kx))/nl;

   ii = end_j;
   while (isnan(mean(area(ii-20:ii))) | (sum(abs(area(ii-20:ii)-mean(area(ii-20:ii))))>(1/nl)))
     ii = ii-10;
   end

   area(ii+1:end_j) = mean(area(ii-20:ii));

  end
end


area = area.';

% fix NAN's
area(find(isnan(area)))=1;

%FIX FOR SPIRAL
area(1) = 0;

area_out = area;
