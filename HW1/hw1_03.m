% p34.m - Allen-Cahn eq. u_t = eps*u_xx+u-u^3, u(-1)=-1, u(1)=1
%         (compare p6.m and p32.m)

%% HW1 - 3
% u_t - mu*u_xx = f, u=0 on bdry (1D)
% Euler: u_n+1 = u_n + dt*(f + mu*u_xx)

g = @(x,x0,a,b) a*exp(-b*(x-x0).^2);

%---Problem parameters---%
N = 20;
mu = 0.01;
dt = 0.01;
Tf = 10;

%---Differentiation---%
[D,x] = cheb(N); D2 = D^2;
I = eye(N+1);
A = I + dt*mu*D2;
A([1 N+1],:) = 0;

%---Forcing---%
f = @(x,t) sin(t).*sin(pi*x);

%---Iterations---%
% Initialization
u = g(x,0,1,10) + g(x,-0.5,0.5,20) + g(x,0.25,1.25,40) + g(x,0.75,0.1,100);
t = 0;

% % Solve PDE by Euler formula and plot results:
  tmax = 100;
  tplot = 0.5;
  nplots = round(tmax/tplot);
  plotgap = round(tplot/dt);
  dt = tplot/plotgap;
  
  % Initial plot data
  xx = -1:.025:1;
  uu = polyval(polyfit(x,u,N),xx);
  plotdata = [uu; zeros(nplots,length(xx))];
  tdata = t;
  
  % Iterations
  for i = 1:nplots
    for n = 1:plotgap
        F = f(x,t); F(1) = 0; F(N+1) = 0;
        u = dt*F + A*u;    
        t = t+dt;    % Euler
    end
    uu = polyval(polyfit(x,u,N),xx);
    plotdata(i+1,:) = uu; tdata = [tdata; t];
  end
  clf, subplot('position',[.1 .4 .8 .5])
  mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -3 3]),
  view(-60,55), colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u
