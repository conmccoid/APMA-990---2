% p2.m - convergence of periodic spectral method (compare p1.m)

% For various N (even), set up grid as before:
  clf, subplot('position',[.1 .4 .8 .5])
  logN = round(logspace(log10(2),2));
%   for N = 2:2:100;
for N = logN;
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
%     u = exp(sin(x)); uprime = cos(x).*u;
%     u = x.^2; uprime = 2*x;
    w = exp(-abs(x - 0).^2/0.1); wprime = -2*(x-0)/0.1.*w;
    v = abs(x); vprime = sign(x);
    u = w.*v; uprime = w.*vprime + wprime.*v;

    % Construct spectral differentiation matrix:
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    D = toeplitz(column,column([1 N:-1:2]));

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf);
    loglog(N,error,'.','markersize',15), hold on
  end
  grid on, xlabel N, ylabel error
  title('Convergence of spectral differentiation')

  figure(2)
  plot(x,D*u,x,uprime)