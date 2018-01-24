function [out] = r_refine_function(in)
   B=in(end);
   out=-2*(0.5/B*in.^2-in);
end