function vo = prepData(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

if a.is.empty
    error empty
end

AA = a.A;
%sen = a.VS;
V = a.V;
num_coils = a.num_coils; %size(sen,2);  % NOTE THIS IS RANK oF SVD Combination
rank_svd = a.rank_svd;
nn = size(AA);
nx = nn(1);
ny = nn(2);


if (size(vi,1) ~= (num_coils*nx))
    error('sense_svd/prepData: Error in prepData for coil combination SVD')
end


vi = reshape(vi,[nx num_coils]);


vo = vi*V(:,1:rank_svd);

vo = vo(:);

