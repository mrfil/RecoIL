function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

if a.is.empty
    error empty
end

AA = a.A;
sen = a.VS;
ncoils = size(sen,2);
nn = size(AA);
nx = nn(1);
ny = nn(2);


if ~a.is.transpose
    vo = zeros(ncoils*nx,1);
    for ii = 1:ncoils
        vo((ii-1)*nx+1:ii*nx) = AA*(sen(:,ii).*vi(:));
    end
else
    vo = zeros(ny,1);
    for ii = 1:ncoils
        vo = vo + (sen(:,ii)'.').*(AA'*vi((ii-1)*nx+1:ii*nx));
    end
end