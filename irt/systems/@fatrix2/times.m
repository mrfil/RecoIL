 function c = times(obj, b)
%function c = times(obj, b)
% "times" method for this class

if isempty(obj.handle_times)
	error 'no abs() method for this object'
end

c = feval(obj.handle_times, obj, b);
