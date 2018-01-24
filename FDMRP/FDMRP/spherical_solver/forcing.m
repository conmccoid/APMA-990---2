function [f_mat]=forcing(phi,theta,r)
    l_r=length(r);   
    l_theta=length(theta); 
    l_phi=length(phi);     
    f_mat=zeros(l_phi,l_theta,l_r);
    
    for k=1:1:l_r
        for j=1:1:l_theta
            for i=1:1:l_phi
                f_mat(i,j,k)=0;
            end
        end
    end
    
end