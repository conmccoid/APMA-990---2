function [FRONT, BACK, LEFT, RIGHT]=neumann_boundary(xs,ys)
        N=length(xs);           
        M=length(ys);           
        LEFT=zeros(1,M);        
        RIGHT=zeros(1,M);
        FRONT=zeros(1,N);
        BACK=zeros(1,N);

        %FRONT NEUMANN BOUNDARY (y=0)
        for i=1:1:N
            x=xs(i);
            FRONT(i)=-2*sin(-2*(x-.5)*pi)*(x-.5)*pi*(3*exp(-100*(x-.95)^2-10)+1)+60*cos(-2*(x-.5)*pi)*exp(-100*(x-.95)^2-10);
        end
        
        %BACK NEUMANN BOUNDARY (y=Ly)
        for i=1:1:N
            x=xs(i);
            BACK(i)=0;
        end
        
        %LEFT NEUMANN BOUNDARY (x=0)
        for j=1:1:M
            y=ys(j);
            LEFT(j)=-2*sin(-1.0*(y-1)*pi)*(y-1)*pi*(3*exp(-90.2500-10*(y-1)^2)+1)+570.00*cos(-1.0*(y-1)*pi)*exp(-90.2500-10*(y-1)^2);
        end
 
        %RIGHT NEUMANN BOUNDARY (x=Lx)
        for j=1:1:M
            y=ys(j);
            RIGHT(j)=0;
        end
end