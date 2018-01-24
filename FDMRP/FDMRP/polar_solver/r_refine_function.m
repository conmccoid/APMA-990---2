function [out] = r_refine_function(in)
%refining function must start at 0 and end on R_max
%out=in;
B=in(end);
out=-2*(0.5/B*in.^2-in);
end