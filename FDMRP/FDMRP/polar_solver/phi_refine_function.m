function [out] = phi_refine_function(in)
%refining function must start at 0 and end on 2*pi

A=in(1);
B=in(end);
out=in;%2*pi*(4/(B-A)^3)*((in-(B+A)/2).^3+((B-A)/2)^3);


end
