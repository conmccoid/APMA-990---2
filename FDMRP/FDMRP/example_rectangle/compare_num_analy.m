function [xs ys u_diff u_a] = compare_num_analy(xs,ys,u)
[N M]=size(u);
u_a=zeros(N,M);
for i=1:1:N
    for j=1:1:M
        x=xs(i,j);
        y=ys(i,j);
        u_a(i,j)=cos(2*(y-1)*(x-.5)*pi)*(3*exp(-100*(x-.95)^2-10*(y-1)^2)+1);
    end
end

u_diff=abs(u_a-u);



end
