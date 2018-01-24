function [xs ys u relres iter resvec] = square_2d_poisson(x_max,y_max,num_xs,num_ys,ref_xs,ref_ys,DiriBC,NeuBC,forcing,useiter)
%   square_2d_poisson - solve Poisson's equation in 2d Cartesian coordinates
%
%   [xs ys u relres iter resvec] = square_2d_poisson(x_max,y_max,num_xs,num_ys,ref_xs,ref_ys,DiriBC,NeuBC,forcing,useiter)
%   is used to solve the boundary value problem in Cartesian coordinates in a square domain. 
%   An large system of equations is solved using the matrix equation where the system of 
%   equations is composed of one equation for every point in the domain. Note: Dirichlet 
%   Boundary Conditions overwrite Neumann Boundary Conditions! ie) Dirichlet BC at x=1 
%   that is not set to NaN in the DiriBC function will overwrite a Neumann BC at x=1 that 
%   is set by the NeuBC function.
%
% INPUTS:
%   -x_max is the length of the x side
%
%   -y_max is the length of the y side
%
%   -num_xs is the number of points to take in the x direction
%
%   -num_ys is the number of points to take in the y direction
%
%   -ref_xs is the function that defines the location of points in x. It
%   must be passed as a string ex) 'x_refine_function' (with quotes)
%
%   -ref_ys is the function that defines the location of points in y. It
%   must be passed as a string ex) 'y_refine_function' (with quotes)
%   
%   -DiriBC is a function that generates the Dirichlet Boundary Conditions. It
%   is passed as a string ex) 'dirichlet_boundary' (with quotes)
%   
%   -NeuBC is a function that generates the Neumann Boundary Conditions. It
%   is passed as a string ex) 'neumann_boundary' (with quotes)
%
%   -forcing is a function that generates the forcing term. It is passed as a 
%   string ex) 'forcing' (with quotes)
%
%   -useiter is an optional parameter. If not set, it defaults to 0 meaning
%   mldivide is used to solve the matrix equation. Any value other than 0
%   tells the solver to use an iterative method.
%
% OUTPUTS
%   xs and ys are the cartesian coordinates for the u values.
% 
%   if useiter~=0, the iterative solver is used and relres final residual, iter 
%   is the number of iterations needed, and resvec is a vector of the history of 
%   the residual
%
% Run Example
%   [xs ys u relres iter resvec] = square_2d_poisson(1,2,100,100,'x_refine_function','y_refine_function','dirichlet_boundary','neumann_boundary','forcing');
%   figure(1), surf(xs,ys,u,'EdgeColor', 'none'), view(-35,40);%, daspect([1 1 5]);

%   [xs ys u_diff u_a] = compare_num_analy(xs,ys,u);
%   figure(2), surf(xs,ys,u_diff,'EdgeColor', 'none'), view(-35,40);%, daspect([1 1 5]);

%   Ashton S. Reimer & Alexei F. Cheviakov
%   Copyright 2011
%   Permission is granted to copy, distribute and/or modify this document
%   under the terms of the GNU General Public License, Version 3.0
%   or any later version published by the Free Software Foundation;
%   with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

%   $Version: 1.0 $  $Date: 2011/02/15 21:10 $

if ( nargin < 10 || isempty(useiter) )
    useiter=0;
end

    tic;                                %start timer
    J=x_max;                            %simplify notation for later (max x value)
    K=y_max;                            %simplify notation for later (max y value)
    N=num_xs;                           %number of points in x     
    M=num_ys;                           %number of points in y
    xx=((1:1:N)-1)/(N-1)*J;             %points used to calculate space step in x
    yy=((1:1:M)-1)/(M-1)*K;             %points used to calculate space step in y

% *************************************************************************
% DEFINE POINT LOCATIONS

    ref_x_handle=str2func(ref_xs);              %define handler for ref_xs function
    ref_y_handle=str2func(ref_ys);              %define handler for ref_ys function
    ys=ref_y_handle(yy);                        %call ref_xs function
    xs=ref_x_handle(xx);                        %call ref_ys function

% *************************************************************************
% CHANGE THESE TO ALTER THE SPACE STEP.             
    
%Initialize arrays to hold space steps
    hx=zeros(N-1,1);
    hy=zeros(M-1,1);
%Calculate the space steps in x
    for i=1:1:N-1
        hx(i)=xs(i+1)-xs(i);
    end
%and in y
    for j=1:1:M-1
        hy(j)=ys(j+1)-ys(j);
    end
    
% *************************************************************************
% SET UP THE A MATRIX FOR MATRIX EQUATION A*x=R   

    w=zeros(N,1);           %coefficient variables
    A_m=zeros(N,1);
    A_p=zeros(N,1);
    B_m=zeros(N,1);
    B_p=zeros(N,1);
    i=1:1:N;
    hx_i(1:1:N-1)=hx(1:1:N-1);
    hx_i(N)=hx(N-1);
    hx_im1(2:1:N)=hx(1:1:N-1);
    hx_im1(1)=hx(1);
        for j=1:1:M
            ind=i+(j-1)*N;
            if (j==M)
                hy_j=hy(j-1);
            else
                hy_j=hy(j);
            end
            if (j==1)
                hy_jm1=hy(j);
            else
                hy_jm1=hy(j-1);
            end
            %Constants for finite difference formula
            w(ind)=(2./(hx_i+hx_im1)).*(-1./hx_i+-1./hx_im1)+(2/(hy_j+hy_jm1))*(-1/hy_j+-1/hy_jm1);
            A_m(ind)=(2./(hx_i+hx_im1)).*(1./hx_im1);
            A_p(ind)=(2./(hx_i+hx_im1)).*(1./hx_i);
            B_m(ind)=(2/(hy_j+hy_jm1))*(1/hy_jm1);
            B_p(ind)=(2/(hy_j+hy_jm1))*(1/hy_j);
            
        end
    
    %Initialize variables for setting up A*x=R
    u=zeros(N,M);               %Grid of solution
    FRONT_N=u;                  %Neumann BC grid
    BACK_N=FRONT_N;             %Neumann BC grid
    LEFT_N=u;                   %Neumann BC grid
    RIGHT_N=LEFT_N;             %Neumann BC grid
    FRONT_D=u;                  %Dirichlet BC grid
    FRONT_D(1:1:N,1:1:M)=NaN;   %Make all Dirichlet
    BACK_D=FRONT_D;             %Dirichlet BC grid
    LEFT_D=BACK_D;              %Dirichlet BC grid
    RIGHT_D=LEFT_D;             %Dirichlet BC grid
    P=N*M;                      %Total number of points
    R=zeros(P,1);               %R vector for equation A*x=R

% *************************************************************************
% CHANGE THESE TO ALTER THE SPACE STEP.       

    diriBC_handle=str2func(DiriBC);                     %create handler for DiriBC function
    [FRONT, BACK, LEFT, RIGHT]=diriBC_handle(xs,ys);    %call DiriBC function
    FRONT_D(:,1)=FRONT;                                 %Assign output to appropriate locations
    BACK_D(:,M)=BACK;
    LEFT_D(1,:)=LEFT;
    RIGHT_D(N,:)=RIGHT;

% *************************************************************************
% CHANGE THESE TO ALTER THE SPACE STEP.       

    NeuBC_handle=str2func(NeuBC);                       %create handler for NeuBC function
    [FRONT, BACK, LEFT, RIGHT]=NeuBC_handle(xs,ys);     %call NeuBC function
    FRONT_N(:,1)=FRONT;                                 %Assign output to appropriate locations
    BACK_N(:,M)=BACK;
    LEFT_N(1,:)=LEFT;
    RIGHT_N(N,:)=RIGHT;
    
% *************************************************************************
    
    disp(['BCs set up in ',num2str(toc),'s']);
    
    %Place coefficients into appropriate locations in A matrix
    A=sparse(P,P);
    diag0=(1:1:P);
    A=A+sparse(diag0,diag0,w,P,P);          %MAIN DIAGONAL OF A
    diag1=1:1:P-1;
    diag2=1:N:P-1;
    diag3=N:N:P-1;
    diag4=N-1:N:P-1;

    A=A+sparse(diag1,diag1+1,A_p(diag1),P,P);               %X+1 Direction
    A=A+sparse(diag2,diag2+1,A_m(diag2),P,P);               %X+1 Direction NEUMANN MOD
    A=A+sparse(diag3,diag3+1,-A_p(diag3),P,P);              %X+1 Direction Entry Corrector
    A=A+sparse(diag1+1,diag1,A_m(diag1+1),P,P);             %X-1 Direction
    A=A+sparse(diag4+1,diag4,A_p(diag4+1),P,P);             %X-1 Direction NEUMANN MOD
    A=A+sparse(diag3+1,diag3,-A_m(diag3+1),P,P);            %X-1 Direction Entry Corrector
%**************************************************************************
%2nd Outer Diagonals
% These two diagonals take care of the Y+1 direction (right diagonal in A 
% matrix) and the Y-1 direction (left diagonal in A matrix)
%**************************************************************************
    diag5=N+1:1:P;
    diag6=N+1:1:2*N;

    A=A+sparse(diag5-N,diag5,B_p(diag5-N),P,P);                        %Y+1 Direction
    A=A+sparse(diag6-N,diag6,B_m(diag6-N),P,P);                        %Y+1 Direction NEUMANN MOD
    A=A+sparse(diag5,diag5-N,B_m(diag5),P,P);                          %Y-1 Direction
    A=A+sparse(P-(diag6-N-1),P-diag6+1,B_p(P-(diag6-N-1)),P,P);        %Y-1 Direction NEUMANN MOD

%SET UP R vector for the system A*x=R.
    R=R+2*A_m(1)*hx(1)*reshape(LEFT_N,P,1);
    R=R-2*A_p(end)*hx(end)*reshape(RIGHT_N,P,1);
    R=R-2*B_p(end)*hy(end)*reshape(BACK_N,P,1);
    R=R+2*B_m(1)*hy(1)*reshape(FRONT_N,P,1);
    disp(['A matrix done in ',num2str(toc)]);
% *************************************************************************
% FORCING

    forcing_handle=str2func(forcing);           %create handler for forcing function
    f=forcing_handle(xs,ys);                    %call forcing function

% *************************************************************************

    F=reshape(f,P,1);
    R=R+F;
    disp(['Forcing calculated in ',num2str(toc)]);

% *************************************************************************
% Apply Boundary Conditions
    any_dirichlet=0;
    ind=find(~isnan(FRONT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=FRONT_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    ind=find(~isnan(BACK_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=BACK_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    ind=find(~isnan(LEFT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=LEFT_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    ind=find(~isnan(RIGHT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=RIGHT_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    
% *************************************************************************
% SOLVE MATRIX EQUATION

    disp('Solving...');
    if (useiter==0)
        x=A\R;
        relres=max(abs(R-A*x));
        iter='mat';
        resvec='mat';
    else
        setup.type = 'nofill';
        [L,U] = ilu(A,setup);           % GENERATES PRECONDITIONERS TO SPEED UP bicg ITERATIONS
        [x flag relres iter resvec]=bicg(A,R,0.000001,10000,L,U);
    end
    u=reshape(x,N,M);
    if (any_dirichlet == 0)
        u=u-min(min(u));
    end
    [xs ys]=meshgrid(xs,ys);
    xs=xs';
    ys=ys';
    disp(['Complete in ',num2str(toc), ' secs']);
    
end
