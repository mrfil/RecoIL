function ob = sensemult(image_mult)
%function ob = sensemult(image_mult)
%	Constructs object, which can do Ax and A'y operations
%           Where A is simply the image_mult image on the diagonal
%           Size is square, length of image_mult on each axis.


%	default object
ob.image_mult = 0;	
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'sensemult');
	return
end

if isa(image_mult, 'sensemult')
	ob = image_mult;
	return
end

if nargin ~= 1  %7
	help sensemult
	error nargin
end

	%	fill object
	ob.image_mult = image_mult;

	ob.is.empty	= logical(0);

%	ob.m = size(we);	% image size
%	ob.n = size(we);

	ob = class(ob, 'sensemult');



