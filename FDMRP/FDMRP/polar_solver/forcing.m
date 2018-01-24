function [f_mat]=forcing(phis,rs)
% This function calculates the forcing term in Poisson's Equation    

    l_r=length(rs);                      %calculate number of rs
    l_phi=length(phis);                  %calculate number of phis
    f_mat=zeros(l_phi,l_r);             %initialize arrays to store forcing terms in
    
    for k=1:1:l_r
        r=rs(k);
        for i=1:1:l_phi
            phi=phis(i);
            f_mat(i,k)=-r*(5*cos(phi)^2-7); %2;
        end
    end
    
end