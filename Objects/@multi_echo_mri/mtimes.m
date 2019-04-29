 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

%Currently only good for body coil   sensitivity = 1!!
% Need to make loops over mm.
 
if a.is.empty
	error empty
end


ntp = a.ntp;
A = a.A;
if ~a.is.transpose
      for kk = 1:ntp
	Atp = A{kk};
        
	if kk ==1
        vo = Atp*vi(:);
    else
	    vo = cat(1,vo,Atp*vi(:));
	end
      end  
else 
       start_index = 1;
       end_index = 0;
       
      for kk = 1:ntp
        Atp = A{kk};
        sizeA = size(Atp);
        len = sizeA(1);
        end_index = end_index +len;
     
	if kk ==1
        vo = Atp'*vi(start_index:end_index);
    else
        vo = vo + Atp'*vi(start_index:end_index);
    end
       tp =Atp'*vi(start_index:end_index);
       
       start_index = start_index+len;
%        name = [num2str(kk)];
%        eval(sprintf('save vo_%s vo',name))
%        eval(sprintf('save tp_%s tp',name))
      end
end







