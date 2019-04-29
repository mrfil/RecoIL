function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class
A=ob.A;
length_total = 0;
for kk = 1:ob.ntp
    Atp = A{kk};
    sizeA = size(Atp);
    length_total = length_total + sizeA(1);
end
dim = [length_total,sizeA(2)] ;
