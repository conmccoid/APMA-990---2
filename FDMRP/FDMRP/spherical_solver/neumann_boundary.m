function [neuBC]=neumann_boundary(phi,theta)
    N=length(phi);   
    M=length(theta); 
    neuBC=zeros(N,M);

    for j=1:1:M
        for i=1:1:N
            neuBC(i,j)=0;
            
        end
    end
end