%% DiskFun runthrough

% Laplace on disk

% RHS
f = @(t,r) 0;

% Exact solution
UUk = @(t,r,k) (r.^k).*cos(k*t)/(k*(k-1));
UU = @(t,r) 0;
for k = 2:100
    UU = @(t,r) UU(t,r) + UUk(t,r,k);
end
uu = diskfun(UU,'polar');

% BC
g = @(t) UU(t,1);

% Solution
    % f is the RHS function, g is the Dirichlet BC, last entries is number of
    % points; f and g can be function handles in polar coord.s or
    % matrix/vector evaluated on the grid [-pi,pi]x[0,1]
u = diskfun.poisson( f, g, 100 );

figure(1)
surf(u)
figure(2)
surf(uu)
figure(3)
surf(uu - u)