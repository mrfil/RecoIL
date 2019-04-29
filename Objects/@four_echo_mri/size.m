function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class

dima = size(ob.a1);
dim = [4*dima(1) dima(2)];

