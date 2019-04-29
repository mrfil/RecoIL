function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class

image_mult = ob.image_mult;
tmp = length(image_mult(:));


dim = [tmp tmp];
