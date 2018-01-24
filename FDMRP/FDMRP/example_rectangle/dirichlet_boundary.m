function [FRONT, BACK, LEFT, RIGHT]=dirichlet_boundary(xs,ys)
        N=length(xs);
        M=length(ys);
        LEFT=zeros(1,M);
        RIGHT=zeros(1,M);       
        FRONT=zeros(1,N);
        BACK=zeros(1,N);

        %FRONT DIRICHLET BOUNDARY (y=0)
        for i=1:1:N
            x=xs(i);
            FRONT(i)=NaN;
        end
        
        %BACK DIRICHLET BOUNDARY (y=Ly)
        for i=1:1:N
            x=xs(i);
            BACK(i)=cos(2*(x-.5)*pi)*(3*exp(-100*(x-.95)^2-10)+1);
        end
        
        %LEFT DIRICHLET BOUNDARY (x=0)
        for j=1:1:M
            y=ys(j);
            LEFT(j)=NaN;
        end
 
        %RIGHT DIRICHLET BOUNDARY (x=Lx)
        for j=1:1:M
            y=ys(j);
            RIGHT(j)=cos(1.0*(y-1)*pi)*(3*exp(-.2500-10*(y-1)^2)+1);
        end
end