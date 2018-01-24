function [out] = theta_refine_function(in)
   A=in(1);
   B=in(end);
   out=pi*(4/(B-A)^3)*((in-(B+A)/2).^3+((B-A)/2)^3);
end