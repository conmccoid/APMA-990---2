function [xs ys u rs phis relres iter resvec] = polar_2d_poisson(R_max,num_rs,num_phis,ref_rs,ref_phis,DiriBC,NeuBC,forcing,useiter)
%   polar_2d_poisson - solve Poisson's equation in 2d polar coordinates
%
%   [xs ys u rs phis relres iter resvec] = polar_2d_poisson(R_max,num_rs,num_phis,ref_rs,ref_phis,DiriBC,NeuBC,forcing,useiter)
%   is used to solve the boundary value problem in polar coordinates in a circular domain. 
%   An large system of equations is solved using the matrix equation where the system of 
%   equations is composed of one equation for every point in the domain. Note: Dirichlet 
%   Boundary Conditions overwrite Neumann Boundary Conditions! ie) Dirichlet BC at theta=pi 
%   that is not set to NaN in the DiriBC function will overwrite a Neumann BC at theta=pi
%   that is set by the NeuBC function.
%
% INPUTS:
%   -R_max is the radius of the circle
%
%   -num_rs is the number of points to take in radius r
%
%   -num_phis is the number of points to take in polar angle phi
%
%   -ref_rs is the function that defines the location of points in r. It
%   must be passed as a string ex) 'r_refine_function' (with quotes)
%
%   -ref_phis is the function that defines the location of points in phi. It
%   must be passed as a string ex) 'phi_refine_function' (with quotes)
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
%   xs and ys are the cartesian coordinates for the u values, rs and phis are 
%   the polar coordinates for the u values.
% 
%   if useiter~=0, the iterative solver is used and relres is the final residual, iter 
%   is the number of iterations needed, and resvec is a vector of the history of 
%   the residual
%
% Run Example
%   [xs ys v relres iter resvec] = polar_2d_poisson(1,100,200,'r_refine_function','phi_refine_function','dirichlet_boundary','neumann_boundary','forcing');
%   surf(xs,ys,v,'EdgeColor','none')

%   Ashton S. Reimer & Alexei F. Cheviakov
%   Copyright 2011
%   Permission is granted to copy, distribute and/or modify this document
%   under the terms of the GNU General Public License, Version 3.0
%   or any later version published by the Free Software Foundation;
%   with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

%   $Version: 1.0 $  $Date: 2011/08/21 18:48 $

if ( nargin < 9 || isempty(useiter) )
    useiter=0;
end

    tic;                                %start timer
    N=num_phis;                         %number of points to use in phi
    O=num_rs;                           %number of points to use in r
    P=N*O;                              %number of points in total
    phi=((1:1:N)-1)/(N)*2*pi;           %points used to calculate space step in phi
    r=((1:1:O)-1)/(O-1)*R_max;          %points used to calculate space step in r

% *************************************************************************
% DEFINE POINT LOCATIONS

    ref_r_handle=str2func(ref_rs);      %define handler for ref_rs function
    ref_phi_handle=str2func(ref_phis);  %call ref_rs function
    phis=ref_phi_handle(phi);           %define handler for ref_phis function
    rs=ref_r_handle(r);                 %call ref_phis function

% *************************************************************************
%  CALCULATE SPACE STEP

%Initialize arrays to hold space steps
    hphi=zeros(N-1,1);
    hr=zeros(O-1,1);
%Calculate the space steps in phi
    for j=1:1:N-1
        hphi(j)=phis(j+1)-phis(j);
    end
%and in r
    for k=1:1:O-1
        hr(k)=rs(k+1)-rs(k);
    end
    
% *************************************************************************

    A=sparse(P,P);                      %A matrix initialization
    R=zeros(P,1);                       %R vector for A*x=R
    
% *************************************************************************
% FORCING

    forcing_handle=str2func(forcing);   %define forcing function handle
    f_mat=forcing_handle(phis,rs);      %call forcing function
    f=reshape(f_mat,P,1);               %convert f_mat matrix into column vector f
    
    disp(['Initialized Points, Space Steps, and Forcing in ',num2str(toc)]);
% *************************************************************************    


    disp(['Setting up A matrix...']);

% *************************************************************************
% SET UP THE A MATRIX FOR MATRIX EQUATION A*x=R   


    %***********
    %First we are going to generate the coefficients for the matrix
    %elements as dictated by the finite-difference equation
    %***********
    % 1) Inside: 0<r<R_max 
    %**********************

    ind=0;                      %index variables (locations to place
    I_00=zeros((O-2)*N,1);      %coefficient variables)
    I_p0=zeros((O-2)*N,1);
    I_m0=zeros((O-2)*N,1);
    I_0p=zeros((O-2)*N,1);
    I_0m=zeros((O-2)*N,1);
    w=zeros((O-2)*N,1);         %coefficient variables
    A_p=zeros((O-2)*N,1);
    A_m=zeros((O-2)*N,1);
    C_p=zeros((O-2)*N,1);
    C_m=zeros((O-2)*N,1);
    r=zeros((O-2)*N,1);
    for k=2:1:O-1
        for i=1:1:N
            ind=ind+1;
            r(ind)=rs(k);
            I_00(ind)=i+(k-1)*N;
            hr_k=hr(k);                                     %space step at r
            hr_km1=hr(k-1);                                 %space step at r-1
            if(i==N) 
                I_p0(ind)=1+(k-1)*N;                        %Makes an index for upper diagonal nearest main
                hphi_i=hphi(i-1);                           %Make the last two space steps equivalent (then we have same # of space steps as points)
            else
                I_p0(ind)=i+1 +(k-1)*N;                     %Makes an index for upper diagonal nearest main
                hphi_i=hphi(i);                             %Space step in polar angle
            end

            if(i==1)
                I_m0(ind)=N+(k-1)*N;                        %Makes an index for lower diagonal nearest main
                hphi_im1=hphi(i);                           %Make the first two space steps equal
            else
                I_m0(ind)=i-1+(k-1)*N;                      %Makes an index for lower diagonal nearest main
                hphi_im1=hphi(i-1);                         %Space step in polar angle
            end

            I_0p(ind)=i+(k)*N;
            I_0m(ind)=i+(k-2)*N;

            %Constants for finite difference formula
            w(ind)=-((2/(hr_k+hr_km1))*(1/hr_k+1/hr_km1)*r(ind)^2 + (2/(hphi_i+hphi_im1))*(1/hphi_i+1/hphi_im1));       %Main diagonal of A matrix
            A_p(ind)=(2/(hr_k+hr_km1))*(1/hr_k)*r(ind)^2 + 1/(hr_k+hr_km1)*r(ind);                                      %Diagonal just above main
            A_m(ind)=(2/(hr_k+hr_km1))*(1/hr_km1)*r(ind)^2- 1/(hr_k+hr_km1)*r(ind);                                     %Diagonal just below main
            C_p(ind)=(2/(hphi_i+hphi_im1))*(1/hphi_i);                                                                  %Upper Interaction layer Diagonal
            C_m(ind)=(2/(hphi_i+hphi_im1))*(1/hphi_im1);                                                                %Lower Interaction layer Diagonal

        end        
    end    
    
    %Place coefficients into proper locations in A matrix
    A=A+sparse(I_00,I_00,w,P,P);
    A=A+sparse(I_00,I_p0,C_p,P,P);
    A=A+sparse(I_00,I_m0,C_m,P,P);
    A=A+sparse(I_00,I_0m,A_m,P,P);
    A=A+sparse(I_00,I_0p,A_p,P,P);
    R(I_00)=r.^2.*f(I_00);
    disp(['  Inside done in ',num2str(toc)]);

    
    % 2) Middle: r=0
    % This requires special boundary conditions due to the Jacobian being
    % singular at r=0 (see paper for more details)
    %**********************
    r=0;                        %This is where we are now.
    hphiss=hphi;                %polar angle space step
    hphiss(N)=hphi(end);        %set the last space step to be the same as the second last
    I_00=1+(1-1)*N;             %only for phi=0. Set others equal later.
	
    dA=pi*(0.5*hr(1))^2;        %set up one of the r=0 coefficients
    dL=(hr(1)/2)*hphiss;        %set up one of the r=0 coefficients
    coeff_u00=-sum(dL)/(1);     %set up one of the r=0 coefficients
    coeff_etc=dL/(1);           %set up one of the r=0 coefficients
    
	A=A+sparse(I_00,I_00,coeff_u00,P,P);    %now put one of the coefficients into the A matrix (A*x=R)

    i=1:1:N;                    %array to generate etc coefficients
    I_etc=i+(2-1)*N;            %indicies for etc coefficients
    
    A=A+sparse(I_00,I_etc,coeff_etc,P,P);   %now put the other coefficients into the A matrix (A*x=R)
	
    R(I_00)=dA*hr(1)*f(I_00);   %assemble the R part of A*x=R for r=0 points

	I_axis=1+(1-1)*N;           %at r=0, phis do not matter, so set them all equal
	i=2:1:N;                    
    I_00=i+(1-1)*N;             %indicies for phis at r=0
    A=A+sparse(I_00,I_00,1,P,P);        %main diagonal is 1 and the location of r=0
    A=A+sparse(I_00,I_axis,-1,P,P);     %is -1 to force equality of two things
    R(I_00)=0.0;                        %but only if forcing is also 0.
        
    
    %********
    %Next we initialize variables for the coefficients at r=R_max and 
    %variables to store the Boundary Conditions
    %********
    k=O;
    ind=0;
    I_curr=zeros(N,1);          %index variables (store where to place coefficient
    I_p0=zeros(N,1);            %variables)
    I_m0=zeros(N,1);
    I_0m=zeros(N,1);
    I_0p=zeros(N,1);
    DirBC_val=zeros(N,1);
    w=zeros(N,1);               %coefficient variables
    A_p=zeros(N,1);
    A_m=zeros(N,1);
    C_p=zeros(N,1);
    C_m=zeros(N,1);
    hr_k=hr(end);
    hr_km1=hr(end);
    r=rs(end);
    DiriBC_handle=str2func(DiriBC);     %create handler for DiriBC function
    DiriBC=DiriBC_handle(phis);         %call DiriBC function
    NeumBC_handle=str2func(NeuBC);      %create handler for NeuBC function
    NeumBC=NeumBC_handle(phis);         %call NeuBC function
    any_dirichlet=0;                    %constant that ends up zero if all BC are neumann
    for i=1:1:N
        ind=ind+1;
        I_curr(ind)=i+(O-1)*N;          %Index for main diagonal
        if(i==N) 
            I_p0(ind)=1+(k-1)*N;        %Index for upper diagonal nearest main
            hphi_i=hphi(i-1);           %space step in polar angle
        else
            I_p0(ind)=i+1 +(k-1)*N;     %Index for upper diagonal nearest main
            hphi_i=hphi(i);             %space step in polar angle
        end
        if(i==1) 
            I_m0(ind)=N+(k-1)*N;        %Index for lower diagonal nearest main
            hphi_im1=hphi(i);           %space step in polar angle
        else
            I_m0(ind)=i-1+(k-1)*N;      %Index for lower diagonal nearest main
            hphi_im1=hphi(i-1);         %space step in polar angle
        end

        I_0m(ind)=i+(k-2)*N;
        I_0p(ind)=I_0m(ind);            

       if (~isnan(DiriBC(ind))) 
            DirBC_val(ind)=DiriBC(ind);
            any_dirichlet=any_dirichlet + 1;
            %Constants if Dirichlet BC is found. (Dirichlet overrides Neumann)
            w(ind)=1;
            A_p(ind)=0;
            A_m(ind)=0;
            C_p(ind)=0;
            C_m(ind)=0;
       else
            %Constants if Neumann BC is found.
            w(ind)=-((2/(hr_k+hr_km1))*(1/hr_km1)*r^2 -1/(hr_k+hr_km1)*r + (2/(hphi_i+hphi_im1))*(1/hphi_i+1/hphi_im1));
            A_p(ind)=0;
            A_m(ind)=(2/(hr_k+hr_km1))*(1/hr_km1)*r^2 - 1/(hr_k+hr_km1)*r;
            C_p(ind)=(2/(hphi_i+hphi_im1))*(1/hphi_i);
            C_m(ind)=(2/(hphi_i+hphi_im1))*(1/hphi_im1);
            DirBC_val(ind)=r^2*f(I_curr(ind))-((2/(hr_k+hr_km1))*r^2+hr_k/(hr_k+hr_km1)*r)*NeumBC(ind);

        end
    end
    %Add new coefficients to A matrix
    A=A+sparse(I_curr,I_curr,w,P,P);
    A=A+sparse(I_curr,I_p0,C_p,P,P);
    A=A+sparse(I_curr,I_m0,C_m,P,P);
    A=A+sparse(I_curr,I_0m,A_m,P,P);
    A=A+sparse(I_curr,I_0p,A_p,P,P);
    R(I_curr)=DirBC_val;
    disp(['A Matrix done in ',num2str(toc)]); 

    %*******
    %Now we can solve the matrix equation
    %*******
    disp(['Solving...']);  
    if (useiter==0)                                     %Use mldivide (pretty quick!)
        x=A\R;
        relres=max(abs(R-A*x));                         %Calculate the residual
        iter='mat';                                     %Since we didn't use iterative method, let
        resvec='mat';                                   %these variables reflect that.
    else                                                %Use iterative method
        disp('     Using qmr with tolerance 0.000001');
        setup.type = 'nofill';                          %variable for preconditioners
        [L,U] = ilu(A,setup);                           % GENERATES PRECONDITIONERS TO SPEED UP bicg ITERATIONS
        [x flag relres iter resvec]=bicg(A,R,0.000001,10000,L,U);
    end
    
    %********
    %Now format the solution for plotting
    %********
    u=(reshape(x,N,O));             %reshape into matrix
    
    %Add extra point at phi=2*pi so that when plotting, one sees completed circle
    phis1=phis;
    phis1(N+1)=2*pi;
    u1=zeros(N+1,O);
    for ii=1:1:N
        u1(ii,:)=u(ii,:);
    end
	u1(N+1,:)=u(1,:);
    u=u1;
    if (any_dirichlet == 0)         %if no Dirichlet BCs, then remove arbitrary constant 
        u=u-(min(min(u)));          %from solution, force minimum value in solution to be 0
    end
    [rr tt]=meshgrid(rs,phis1);     %Create a mesh of points then
    [xs ys]=pol2cart(tt,rr);        %generate xs and ys for plotting

	disp(['Completed in ',num2str(toc),' seconds'])     %Show time needed to finish
    
end


