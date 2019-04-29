function ob = multi_echo_mri(A,num_echoes)
     %function ob = four_echo_mri(A1,A2,A3,A4);
%function ob = four_echo_mri(st, we, Int, tt, ech_times, swt)
% This will be tim_seg_mri object for multiecho spiral,
%   specifically four shot
% st is structure from nufft_init_kb
% we is the field map, OR (field map -i* r2map) 
% Int is the temporal interpolator, L x number time pts.
% tt is the time vector
% ech_times is the echo times for each spiral arranged in
%    a column vector
% swt is either omitted or dr sinc(kx*dr)dr*sinc(ky*dr) 

%	default object
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
%ob.version = 1.0;
ob.A = 0;



if nargin == 0
	ob = class(ob, 'multi_echo_mri');
	return
end

if isa(A, 'multi_echo_mri')
	ob = A;
	return
end

if ((nargin < 2) | (nargin > 4))
	help multi_echo_mri
	error nargin
end



ob.A = A;
ob.ntp = num_echoes;

ob.is.empty	= logical(0);

 

ob = class(ob, 'multi_echo_mri');






