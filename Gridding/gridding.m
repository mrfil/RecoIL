function f = gridding(n, kxnl, sig_t)


% gridding.m
% Brad Sutton
% Aug 10 ,1999
% This will be a program to perform gridding for a 1-D 
%    k-space trajectory.  It will take a non-linear 
%    k-space and the values of the signal at each place
%    and return the values at each point of a 
%    gridded k-space, for Fourier Inversion.
%    Note that sig_t(i) is the value corresponding to
%    the sampled signal at kxnl(i)

% The variables are 
%   Input 
%	 n	is the number of k-space samples desired
%		 that is the size of kxdes 
%        kxnl	is the nonlinear k-space (1-D) trajectory
%	 sig_t  is the FID, or the k-space values
%	 kxdes  is the desired k-space trajectory, linear
%	 sigdes is the gridded k-space signal
%	 u	is an indexing variable

% First we initialize the desired k-space trajectory vector
kxdes=linspace(min(kxnl), max(kxnl), n);

% Now we need to assign values to these k-space points 
%    by interpolation

% Need to initialize 
sigdes = zeros(1,n);
u = zeros(1, n);

for i = 1:n
	if kxdes(i) == kxnl(i)
		sigdes(i)=sig_t(i);
	else 
		for k=1:n
			if kxnl(k)<kxdes(i)
				u(i)=k;
				% kxnl(u(i)+1)>kxdes(i)
			else
				;
			end
		end

		sigdes(i)=(sig_t(u(i)+1)*(kxdes(i)-kxnl(u(i)))...
			+sig_t(u(i))*(kxnl(u(i)+1)-kxdes(i)))...
			/(kxnl(u(i)+1)-kxnl(u(i)));
	end	
end

% Now to output the result
f=sigdes;