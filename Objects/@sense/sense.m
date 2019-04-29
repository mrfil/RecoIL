function ob = sense(A,sen)
%function ob = sense(A,sen)
%	Construct MRI object, which can do Ax and A'y operations
%  sen is nptsxncoils

% This is an attempt to make MRI object valid for sensitivity
%   encoded runs also


%	default object
ob.A = 0;	% should these be []'s
ob.sen = 0;
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'sense');
	return
end

if isa(A, 'sense')
	ob = A;
	return
end

if nargin ~= 2  %7
	help sense
	error nargin
end

	%	fill object
	ob.A = A;
	ob.sen = sen;

	ob.is.empty	= logical(0);

%	ob.m = size(we);	% image size
%	ob.n = size(we);

	ob = class(ob, 'sense');



