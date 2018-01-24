function [out]=dirichlet_boundary(in)
% This function defines what part of the boundary will be Dirichlet.
% Anywhere that one does not want the boundary to be dirichlet, one must
% set equal to NaN (not a number). The solver will check to see if the
% dirichlet boundary is set to NaN. If it is, then the solver applies a
% Neumann boundary condition, else it applies a Dirichlet boundary
% condition.

out=zeros(1,length(in));
    for i=1:1:length(in)
            out(i)=NaN;
    end
end