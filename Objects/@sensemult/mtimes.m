 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

if a.is.empty
	error empty
end

image_mult = a.image_mult;



if ~a.is.transpose
       vo = image_mult(:).*vi(:);  
else
       vo = conj(image_mult(:)).*vi(:);
end
