 function vo = prepData(a, vi)
%	Trivial prepData function for non-coil-combination SENSE object

if a.is.empty
	error empty
end


vo = vi(:);
