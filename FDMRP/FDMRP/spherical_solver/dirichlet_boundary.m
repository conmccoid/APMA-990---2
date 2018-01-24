function [diriBC]=dirichlet_boundary(phi,theta)
    N=length(phi);
    M=length(theta);
    diriBC=zeros(N,M); 

    root1=pi/3;  root2=pi;
    epsilon=0.1;
    for j=1:1:M
        for i=1:1:N

            if ((theta(j)<pi/2+epsilon)&&(theta(j)>pi/2-epsilon)&&((phi(i)<root1+epsilon)&&(phi(i)>root1-epsilon)))
                diriBC(i,j)=1;
            elseif ((theta(j)<pi/2+epsilon)&&(theta(j)>pi/2-epsilon)&&((phi(i)<root2+epsilon)&&(phi(i)>root2-epsilon)))
                diriBC(i,j)=-10;
            else
                diriBC(i,j)=NaN;
            end
        end
    end
end