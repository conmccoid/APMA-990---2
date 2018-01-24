%% DiskFun runthrough

% Laplace on disk

% Problem parameters
N = 24; % dimension in r
M = 24; % dimension in t
r = linspace(0,1,N);
t = linspace(0,2*pi,M);

gk = @(t,k) cos(k.*t)./(k.*(k-1)); % pieces of BC

% [tt,kk] = meshgrid(2:100,t);
% G = gk(tt,kk);
% g = G*ones(99,1);