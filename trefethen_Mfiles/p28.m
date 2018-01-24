% p28.m - eigenmodes of Laplacian on the disk (compare p22.m)

% r coordinate, ranging from -1 to 1 (N must be odd):
  N = 25; N2 = (N-1)/2;
  [D,r] = cheb(N);
  D2 = D^2;
  D1 = D2(1:N2+1,1:N2+1); D2 = D2(1:N2+1,N+1:-1:N2+2);
  E1 =  D(1:N2+1,1:N2+1); E2 =  D(1:N2+1,N+1:-1:N2+2);
  
% t = theta coordinate, ranging from 0 to 2*pi (M must be even):
  M = 20; dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2;
  D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
                 .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
             
  % BCs
  D1(1,:) = 0; D1(1,1) = 1;
  D2(1,:) = 0;
  k = 2:100;
  coeffs = 1./(k'.*(k'-1));
  G = cos(t*k)*coeffs;

% Laplacian in polar coordinates:
  R = diag(1./r(1:N2+1)); R(1,1) = 0;
  Z = zeros(M2); I = eye(M2);
  L = kron(D1+R*E1,eye(M)) + kron(D2+R*E2,[Z I;I Z]) ...
                           + kron(R^2,D2t);       

% Plot eigenmodes with nodal lines underneath:
  [rr,tt] = meshgrid(r(1:N2+1),[0;t]);
  [xx,yy] = pol2cart(tt,rr);
  
  [rr,tt] = meshgrid(r(1:N2+1),t);
  rr = rr(:); tt = tt(:);
  f = zeros(size(rr));
  f(1:M) = G;
  u = L\f;
  
  % Exact solution
  uu = bsxfun(@power,rr,k);
  uu = uu.*cos(tt*k);
  uu = uu * coeffs;
  
  % Reshape
  u = reshape(u,M,N2+1);
  u = u([M 1:M],:);
  mesh(xx,yy,u)
  
  uu = reshape(uu,M,N2+1);
  uu = uu([M 1:M],:);
%   mesh(xx,yy,abs(uu - u))
