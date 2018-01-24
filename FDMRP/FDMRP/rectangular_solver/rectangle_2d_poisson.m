function [xs ys u relres iter resvec] = rectangle_2d_poisson(x_max,y_max,num_xs,num_ys,ref_xs,ref_ys,DiriBC,NeuBC,forcing,useiter)
%   rectangle_2d_poisson - solve Poisson's equation in 2d Cartesian coordinates
%
%   [xs ys u relres iter resvec] = rectangle_2d_poisson(x_max,y_max,num_xs,num_ys,ref_xs,ref_ys,DiriBC,NeuBC,forcing,useiter)
%   is used to solve the boundary value problem in Cartesian coordinates in a rectangular
%   domain. An large system of equations is solved using the matrix equation where the 
%   system of equations is composed of one equation for every point in the domain. 
%   Note: Dirichlet Boundary Conditions overwrite Neumann Boundary Conditions! ie) Dirichlet 
%   BC at x=1 that is not set to NaN in the DiriBC function will overwrite a Neumann BC at 
%   x=1 that is set by the NeuBC function.
%
% INPUTS:
%   -x_max is the length of the x side referenced from the origin
%
%   -y_max is the length of the y side referenced from the origin
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
%   if useiter~=0, the iterative solver is used and relres is the final residual, iter 
%   is the number of iterations needed, and resvec is a vector of the history of 
%   the residual
%
% Run Example
%   [xs ys u relres iter resvec] = rectangle_2d_poisson(1,2,100,100,'x_refine_function','y_refine_function','dirichlet_boundary','neumann_boundary','forcing');
%   surf(xs,ys,u,'EdgeColor', 'none')
%
%   Ashton S. Reimer & Alexei F. Cheviakov
%   Copyright 2011
%   Permission is granted to copy, distribute and/or modify this document
%   under the terms of the GNU General Public License, Version 3.0
%   or any later version published by the Free Software Foundation;
%   with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

%   $Version: 1.0 $  $Date: 2011/08/16 00:04 $

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

    %***********
    %First we are going to generate the coefficients for the matrix
    %elements as dictated by the finite-difference equation
    %***********
    %Initialize coefficient variables
    w=zeros(N,1);          
    A_m=zeros(N,1);
    A_p=zeros(N,1);
    B_m=zeros(N,1);
    B_p=zeros(N,1);
    i=1:1:N;                                        %Use this as an index 
    hx_i(1:1:N-1)=hx(1:1:N-1);                      %space step for points i
    hx_i(N)=hx(N-1);                                %Make the last two space steps equivalent (then we have same # of space steps as points)
    hx_im1(2:1:N)=hx(1:1:N-1);                      %space step for points i-1
    hx_im1(1)=hx(1);                                %Make the last two space steps equivalent (then we have same # of space steps as points)
        for j=1:1:M
            ind=i+(j-1)*N;
            if (j==M)                               %Make the last two space steps equivalent (then we have same # of space steps as points)
                hy_j=hy(j-1);                       %space step for points j
            else
                hy_j=hy(j);
            end
            if (j==1)                               %Make the last two space steps equivalent (then we have same # of space steps as points)
                hy_jm1=hy(j);                       %space step for points j-1
            else
                hy_jm1=hy(j-1);
            end
            %Constants for finite difference formula
            w(ind)=(2./(hx_i+hx_im1)).*(-1./hx_i+-1./hx_im1)+(2/(hy_j+hy_jm1))*(-1/hy_j+-1/hy_jm1);     %Main diagonal of A matrix
            A_m(ind)=(2./(hx_i+hx_im1)).*(1./hx_im1);                                                   %Diagonal just below main
            A_p(ind)=(2./(hx_i+hx_im1)).*(1./hx_i);                                                     %Diagonal just above main
            B_m(ind)=(2/(hy_j+hy_jm1))*(1/hy_jm1);                                                      %Lower Interaction layer Diagonal
            B_p(ind)=(2/(hy_j+hy_jm1))*(1/hy_j);                                                        %Upper Interaction layer Diagonal
            
        end
   
    %********
    %Next we initialize variables to store the Boundary Conditions
    %********
    %Initialize variables for setting up R in A*x=R
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
% Get Dirichlet boundary conditions from passed function DiriBC       

    diriBC_handle=str2func(DiriBC);                     %create handler for DiriBC function
    [FRONT, BACK, LEFT, RIGHT]=diriBC_handle(xs,ys);    %call DiriBC function
    FRONT_D(:,1)=FRONT;                                 %Assign output to appropriate locations
    BACK_D(:,M)=BACK;
    LEFT_D(1,:)=LEFT;
    RIGHT_D(N,:)=RIGHT;

% *************************************************************************
% Get Neumann boundary conditions from passed function NeuBC

    NeuBC_handle=str2func(NeuBC);                       %create handler for NeuBC function
    [FRONT, BACK, LEFT, RIGHT]=NeuBC_handle(xs,ys);     %call NeuBC function
    FRONT_N(:,1)=FRONT;                                 %Assign output to appropriate locations
    BACK_N(:,M)=BACK;
    LEFT_N(1,:)=LEFT;
    RIGHT_N(N,:)=RIGHT;
    
    
    disp(['BCs set up in ',num2str(toc),'s']);          %Let the user know where we are/how long it has taken
% *************************************************************************
    

    %*******
    %Now that we have the BCs and the A matrix elements, let's set up A and
    %R in A*x=R
    %******* 
    %Place coefficients into appropriate locations in A matrix
    A=sparse(P,P);                          %Initialize A as a sparse matrix (Thankfully we can do this!)
    
    %MAIN DIAGONAL OF A
    diag0=(1:1:P);
    A=A+sparse(diag0,diag0,w,P,P);
    
    %THE INNER MOST TWO DIAGONALS (the diags closest to the main diagonal)
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

    %2nd Outer Diagonals
    % These two diagonals take care of the Y+1 direction (right diagonal in A 
    % matrix) and the Y-1 direction (left diagonal in A matrix)
    diag5=N+1:1:P;
    diag6=N+1:1:2*N;
    A=A+sparse(diag5-N,diag5,B_p(diag5-N),P,P);                        %Y+1 Direction
    A=A+sparse(diag6-N,diag6,B_m(diag6-N),P,P);                        %Y+1 Direction NEUMANN MOD
    A=A+sparse(diag5,diag5-N,B_m(diag5),P,P);                          %Y-1 Direction
    A=A+sparse(P-(diag6-N-1),P-diag6+1,B_p(P-(diag6-N-1)),P,P);        %Y-1 Direction NEUMANN MOD

    %SET UP R vector for the system A*x=R using only the Neumann BCs for
    %now
    R=R+2*A_m(1)*hx(1)*reshape(LEFT_N,P,1);
    R=R-2*A_p(end)*hx(end)*reshape(RIGHT_N,P,1);
    R=R-2*B_p(end)*hy(end)*reshape(BACK_N,P,1);
    R=R+2*B_m(1)*hy(1)*reshape(FRONT_N,P,1);
    disp(['A matrix done in ',num2str(toc)]);

    %********
    % Now we add the FORCING to the R vector, so first we must get it from
    % the passed forcing function
    %********
    forcing_handle=str2func(forcing);           %create handler for forcing function
    f=forcing_handle(xs,ys);                    %call forcing function

    F=reshape(f,P,1);                           %Shape the forcing matrix into a vector
    R=R+F;                                      %add forcing to the R vector
    disp(['Forcing calculated in ',num2str(toc)]);

    %********
    % Apply Dirichlet boundary conditions. To do this, we look for entries
    % in FRONT_D, BACK_D, LEFT_D, and/or RIGHT_D that are not NaN. Then we
    % set the corresponding row in the A matrix to be all zeros except for
    % a value of 1 on the main diagonal entry. Then we make the R vector at
    % the same location equal to the Dirichlet boundary condition.
    %********
    any_dirichlet=0;                            %Counts to see if any BCs are Dirichlet (if not we need to subtract 
                                                %arbitrary constant from solution to Neumann problem.
    %FRONT_D
    ind=find(~isnan(FRONT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=FRONT_D(ind);
    A=A';                                       %It is faster to perform COLUMN operations on matricies in Matlab, so
    A(:,ind)=0;                                 %first we transpose A, then we blank the appropriate column (same as row
    A=A';                                       %of A) and then we transpose back
    A=A+sparse(ind,ind,1,P,P);
    %BACK_D
    ind=find(~isnan(BACK_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=BACK_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    %LEFT_D
    ind=find(~isnan(LEFT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=LEFT_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    %RIGHT_D
    ind=find(~isnan(RIGHT_D));
    any_dirichlet=any_dirichlet+~isempty(ind);
    R(ind)=RIGHT_D(ind);
    A=A';
    A(:,ind)=0;
    A=A';
    A=A+sparse(ind,ind,1,P,P);
    
    %*******
    %Now we can solve the matrix equation
    %*******
    disp('Solving...');
    if (useiter==0)                     %Use mldivide (pretty quick!)
        x=A\R;
        relres=max(abs(R-A*x));         %Calculate the residual
        iter='mat';                     %Since we didn't use iterative method, let
        resvec='mat';                   %these variables reflect that.
    else                                %Use iterative method
        setup.type = 'nofill';          %variable for preconditioners
        [L,U] = ilu(A,setup);           % GENERATES PRECONDITIONERS TO SPEED UP bicg ITERATIONS
        [x flag relres iter resvec]=bicg(A,R,0.000001,10000,L,U);
    end
    
    %********
    %Now format the solution for plotting
    %********
    u=reshape(x,N,M);                   %reshape into matrix
    if (any_dirichlet == 0)             %if no Dirichlet BCs, then remove arbitrary constant 
        u=u-min(min(u));                %from solution, force minimum value in solution to be 0
    end
    [xs ys]=meshgrid(xs,ys);            %Generate x and y values for plotting
    xs=xs';                             %format them for plotting
    ys=ys';
    disp(['Complete in ',num2str(toc), ' secs']);   %Show time needed to finish
    
end
