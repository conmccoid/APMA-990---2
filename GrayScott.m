%% Gray-Scott equations in 2D
% Nick Trefethen, April 2016

%%
% (Chebfun Example pde/GrayScott.m)
% [Tags: #Gray-Scott, #spin2]

% ut = eps1 Del(u) + b(1-u) - uv^2
% vt = eps2 Del(v) - dv + uv^2

%---Parameters---%
ep1 = 0.00002;
ep2 = 0.00001;
b = 0.04;
d = 0.1;

%---Domain---%
dom = [-1 1 -1 1];
x = chebfun('x',dom(1:2));
tspan = [0 3500];

%---Operator---%
S = spinop2(dom,tspan);
S.lin = @(u,v) [ep1*lap(u); ep2*lap(v)];
S.nonlin = @(u,v) [b*(1-u)-10*u.*v.^2;-d*v+10*u.*v.^2];
S.init = chebfun2v(@(x,y) 1-exp(-80*((x+.05).^2+(y+.02).^2)), ...
                   @(x,y) exp(-80*((x-.05).^2+(y-.02).^2)),dom);
           
%---Solve---%
tic
u = spin2(S,200,2,'plot','off');
plot(u{2}), view(0,90), axis equal, axis off
time_in_seconds = toc

%% References
%
% [1] P. Gray and S. K. Scott, _Chemical Oscillations and
% Instabilities: Non-linear Chemical Kinetics_, v. 21 of
% International Series of Monographs on Chemistry, 1994.
%
% [2] L. N. Trefethen and K. Embree, editors, article 23 on
% "The Gray-Scott equations",
% _The (Unfinished) PDE Coffee Table Book_,
% `https://people.maths.ox.ac.uk/trefethen/pdectb.html`.
%
% [3] H. Montanelli and N. Bootland, _Solving periodic semilinear stiff PDEs
% in 1D, 2D and 3D with exponential integrators_, submitted, 2016.