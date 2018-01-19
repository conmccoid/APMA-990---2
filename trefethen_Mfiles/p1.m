% p1.m - convergence of fourth-order finite differences

% For various N, set up grid in [-pi,pi] and function u(x):
  Nvec = 2.^(3:12);
  figure(1)
  clf, subplot('position',[.1 .4 .8 .5])
  for N = Nvec
    h = 2*pi/N; x = -pi + (1:N)'*h;
%     u = exp(sin(x)); uprime = cos(x).*u; 
%     u = x.^2; uprime = 2*x;
%     u = abs(x); uprime = ones(size(x)); uprime(x<0) = -1; uprime(x==0) = 0;
%     u = x; uprime = ones(size(x));
%     u = abs(sin(x)); uprime = cos(x).*sign(x);
%     u = exp(abs(sin(x))); uprime = sign(x).*cos(x).*u;
%     u = abs(sin(x/2)); uprime = sign(x)/2.*cos(x/2); % note: this one is
%     diff on bdry but not at zero
    w = exp(-abs(x - 0).^2/0.1); wprime = -2*(x-0)/0.1.*w;
    v = abs(x); vprime = sign(x);
    u = w.*v; uprime = w.*vprime + wprime.*v;

    % Construct sparse fourth-order differentiation matrix:
    e = ones(N,1);
    D =   sparse(1:N,[2:N 1],2*e/3,N,N)...
        - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    % Plot max(abs(D*u-uprime)):
    figure(1)
    error = norm(D*u-uprime,inf);
    loglog(N,error,'.','markersize',15), hold on
    
    figure(2)
    plot(x,D*u,x,uprime)
    pause(0.5)
  end
  figure(1)
  grid on, xlabel N, ylabel error
  title('Convergence of fourth-order finite differences')
  semilogy(Nvec,Nvec.^(-4),'--') 
  text(105,5e-8,'N^{-4}','fontsize',18)
  
  figure(2)
  plot(x,D*u,x,uprime)