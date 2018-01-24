function [out]=neumann_boundary(in)
% This function defines what part of the boundary will be Neumann.
% The solver will check to see if the dirichlet boundary is set to NaN. If 
% there is no Dirichlet BC set at the boundary, then the solver applies a 
% Neumann boundary condition, else it applies a Dirichlet boundary condition.
% Thus if one wants a Neumann BC, then one must ensure that the Dirichlet
% BC at that point is set to NaN.

    for k=1:1:length(in)
            out(k)=3*sin(in(k))^2;  %2*sin(in(k))^2;
    end
end