function ob = four_echo_mri(A1,A2,A3,A4)
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
ob.a1 = 0;
ob.a2 = 0;
ob.a3 = 0;
ob.a4 = 0;


if nargin == 0
	ob = class(ob, 'four_echo_mri');
	return
end

if isa(A1, 'four_echo_mri')
	ob = A1;
	return
end

if ((nargin < 4) | (nargin > 4))
	help four_echo_mri
	error nargin
end

%flag_swt = 0;
%if (nargin == 6)
%   flag_swt = 1;
%end

%if ~(size(ech_times,1)==4)
%     sprintf('Size of echo times vector must be 4x1')
%     keyboard
%end

%	%	fill object
%       ob.len = length(tt);
%       if flag_swt
%          ob.a1 = tim_seg_mri(st, we, Int, tt + ech_times(1),swt);
%          ob.a2 = tim_seg_mri(st, we, Int, tt + ech_times(2),swt);
%          ob.a3 = tim_seg_mri(st, we, Int, tt + ech_times(3),swt);
%          ob.a4 = tim_seg_mri(st, we, Int, tt + ech_times(4),swt);
%       else
%          ob.a1 = tim_seg_mri(st, we, Int, tt + ech_times(1));
%          ob.a2 = tim_seg_mri(st, we, Int, tt + ech_times(2));
%          ob.a3 = tim_seg_mri(st, we, Int, tt + ech_times(3));
%          ob.a4 = tim_seg_mri(st, we, Int, tt + ech_times(4));
%        end

ob.a1 = A1;
ob.a2 = A2;
ob.a3 = A3;
ob.a4 = A4;

	ob.is.empty	= logical(0);

 

	ob = class(ob, 'four_echo_mri');






