 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

%Currently only good for body coil   sensitivity = 1!!
% Need to make loops over mm.
 
if a.is.empty
	error empty
end

A1 = a.a1;
A2 = a.a2;
A3 = a.a3;
A4 = a.a4;



if ~a.is.transpose
        vo = [A1*vi(:);A2*vi(:);A3*vi(:);A4*vi(:)]; 
	    
else 
            len = length(vi(:))/4;
	    vo = A1'*vi((1):(len))+A2'*vi((len+1):(2*len))+A3'*vi((2*len+1):(3*len))+A4'*vi((3*len+1):(4*len));
end





