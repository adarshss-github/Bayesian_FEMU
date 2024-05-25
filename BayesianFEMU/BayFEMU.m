function [lamstar,Phistar,PhiLo,thetastar,thetait,COVtheta,eigerr] = BayFEMU(K,Kmat3d,M,Nm,Nit,lamhat,dof,Psihat,thetan,Covmattheta,sigeq2,Covmateps)


%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[lamstar,Phistar,thetastar,thetait,COVtheta] = BayFEMU(K,Kmat3d,M,Nm,Nit,lamhat,dof,Psihat,thetan,Covmattheta,sigeq2,Covmateps)
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description:
%------------
%This function calculates the optimal stiffness parameters and their uncertainities, from a user defined finite element model class, by Bayesian Finite Element Model
%Updating, using an iterative algorithm.The associated non-linear optimization problem is converted to three coupled linear optimization problems first.
%This function also includes the expansion of incomplete modeshapes using Bayesian inference. The data used for Bayesian
%inference includes the squared natural frequencies and incomplete modeshapes.
%
%Input arguments:
%----------------
%1)K: is the function handle of the parametrized stiffness matrix defined
%by the user
%2)Kmat3d: is the 3D matrix containing the stiffness matrix components of
%substructures corresponding to the various model parameters
%3) M: is the user defined mass matrix (unparametrized)
%4)Nm: is the number of modes which were detected
%5)Nit: is the number of iterations
%6)lamhat: is the column vector of measured squared natural frequencies
%7) dof: is the vector containing the degrees of freedom corresponding to
%the FE model, which were used for measurements
%8) Psihat: is the column vector containing the measured partial mode shapes
%9)thetan is the vector containing the nominal values of the stiffness matrix parameters
%10) Covmattheta: is the covariance matrix of the stiffness matrix parameters
%11)sigeq2: is the variance for the error of the eigen value equation
%12) Covmateps: is the covariance matrix of the measured modal data
%13) eigerr: is the cumilative eigen equation error
%
%Output arguments:
%----------------
%1)lamstar: is the updated sqaured natural frequencies
%2)Phistar: is the expanded/updated mode shapes
%3) PhiLo: is the value of the updated modeshapes at the measured DOFs
%4)thetastar: is the updated parameters
%5)thetait: is the matrix comtaining he values of each parameters during
%iteration, arranged in rows
%6)COVtheta: is the coefficient of variation of the parameters in percentage
%
%Ex:[lamstar,Phistar,thetastar,thetait,COVtheta] = BayFEMU(K,Kmat3d,M,Nm,2000,lamhat,dof,Psihat,thetan,Covmattheta,sigeq2,Covmateps) ;
%
%
%Notes:
%------
%Comateps can be found out by using Bayesian modal identification
%algorithms from recorded data
%thetan can be taken from experiments (for small models) or FEM
%Covmattheta should be large enough to account for larger range of possible
%values
%Use high values of Covmattheta and sigeq2 for improper priors
%Warnings are turned off to keep the command window clear, you may switch
%the warnings on
%
%Warning: The Hessian and the coefficient matirces of the linear optimization problems may be ill-conditioned
%--------
%
%Reference: Yuen, K.V., Beck, J.L. and Katafygiotis, L.S., 2006. Efficient model updating and health monitoring methodology using incomplete modal
%data without mode matching. Structural Control and Health Monitoring: The Official Journal of the International Association for Structural Control
%and Monitoring and of the European Association for the Control of Structures, 13(1), pp.91-107.
%
%                                     -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; May, 2020 || *********
%                                     -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

tic
[Nd,~] = size(M) ;
No = length(dof) ;
Ntheta = length(thetan) ;
Lo = zeros(Nm*No,Nm*Nd) ;
dum = 0 ;
warning('off')

%                                  Formation of the gather matrix
%=======================================================================================

for i = 1:1:Nm
    k = 1 ;
    for j = No*(i-1)+1:1:No*(i-1)+No
        Lo(j,dof(k)+dum) = 1 ;
        k = k+1 ;
    end
    dum = dum + Nd ;
    
end

%=======================================================================================

clear dum

lamstar = lamhat ;
Phistar = zeros(Nm*Nd,1) ;
thetastar = thetan ;
Kstar = K(thetastar) ;
Gphi = zeros(Nd*Nm,Nd*Nm) ;
Glam1 = zeros(Nm,Nm) ;
Glam2 = zeros(Nm,1) ;
Gtheta = zeros(Nd*Nm,Ntheta) ;
b = zeros(Nd*Nm,1) ;
L1 = zeros(Nm,Nd*Nm) ;
L2 = zeros(Nm,Ntheta) ;
L3 = zeros(Nd*Nm,Ntheta) ;
invCovmateps = inv(Covmateps) ;
invCovmattheta = inv(Covmattheta) ;
invCovmateps11 = invCovmateps(1:Nm,1:Nm) ;
invCovmateps12 = invCovmateps(1:Nm,Nm+1:end) ;
invCovmateps21 = invCovmateps(Nm+1:end,1:Nm) ;
invCovmateps22 = invCovmateps(Nm+1:end,Nm+1:end) ;
clear invCovmateps Covmattheta
thetait = zeros(Ntheta,Nit) ;

%                   Iterations for the coupled optimization problems
%=====================================================================================

for i = 1:1:Nit
    
    thetait(:,i) = thetastar ;
    
    for j = 1:1:Nm
        Gphi((j-1)*Nd+1:j*Nd,(j-1)*Nd+1:j*Nd) = (lamstar(j)*M - Kstar)^2 ;
    end
    
    Phistar = ((1/sigeq2)*Gphi + Lo'*invCovmateps22*Lo)\(Lo'*(invCovmateps21*(lamhat-lamstar) + invCovmateps22*Psihat)) ;
    
    for j = 1:1:Nm
        Glam1(j,j) = (1/sigeq2)*Phistar((j-1)*Nd+1:j*Nd)'*M^2*Phistar((j-1)*Nd+1:j*Nd) ;
        Glam2(j,1) = (1/sigeq2)*Phistar((j-1)*Nd+1:j*Nd)'*M*Kstar*Phistar((j-1)*Nd+1:j*Nd) ;
    end
    
    lamstar = (Glam1 + invCovmateps11)\( Glam2 + invCovmateps11 *lamhat + invCovmateps12*(Psihat - Lo*Phistar) ) ;
    
    for j = 1:1:Nm
        for k = 1:1:Ntheta
            Gtheta((j-1)*Nd+1:j*Nd,k) = reshape(Kmat3d(:,:,k+1),Nd,Nd)*Phistar((j-1)*Nd+1:j*Nd) ;
        end
        b((j-1)*Nd+1:j*Nd) = (lamstar(j)*M-reshape(Kmat3d(:,:,1),Nd,Nd))*Phistar((j-1)*Nd+1:j*Nd) ;
    end
    
    thetastar = ((1/sigeq2)*Gtheta'*Gtheta + invCovmattheta)\((1/sigeq2)*Gtheta'*b + invCovmattheta*thetan) ;
    
    Kstar = K(thetastar) ;
    
end

%=====================================================================================

%                        Covariance matrix calculation
%=====================================================================================

for i = 1:1:Nm
    
    L1(i,(i-1)*Nd+1:i*Nd) = Phistar((i-1)*Nd+1:i*Nd)'*M*(Kstar-lamstar(i)*M) ;
    
end

for i = 1:1:Nm
    for j = 1:1:Ntheta
        L2(i,j) = Phistar((i-1)*Nd+1:i*Nd)'*M*reshape(Kmat3d(:,:,j+1),Nd,Nd)*Phistar((i-1)*Nd+1:i*Nd) ;
    end
end

for i = 1:1:Nm
    for j = 1:1:Ntheta
        L3((i-1)*Nd+1:i*Nd,j) = ( ( Kstar-lamstar(i)*M )*reshape(Kmat3d(:,:,j+1),Nd,Nd) + reshape(Kmat3d(:,:,j+1),Nd,Nd)*( K(thetastar)-lamstar(i)*M ) )*Phistar((i-1)*Nd+1:i*Nd) ;
    end
end

H = [ Glam1+invCovmateps11 -2*(1/sigeq2)*L1+invCovmateps12*Lo -(1/sigeq2)*L2;
    (-2*(1/sigeq2)*L1+invCovmateps12*Lo)' (1/sigeq2)*Gphi+Lo'*invCovmateps22*Lo (1/sigeq2)*L3;
    (-(1/sigeq2)*L2)' ((1/sigeq2)*L3)' (1/sigeq2)*Gtheta'*Gtheta+invCovmattheta;] ;
Covmat = inv(H) ;
Covmattheta = Covmat(end-Ntheta+1:end,end-Ntheta+1:end) ;
COVtheta = (sqrt(diag(Covmattheta))./thetastar)*100 ;

%=====================================================================================

PhiLo = Lo*Phistar ;

eigerr = 0 ;

%                        Eigen equation error computation
%=====================================================================================

for i = 1:1:Nm
    eigerr = eigerr + (norm((K(thetastar)-lamstar(i)*M)*Phistar((i-1)*Nd+1:i*Nd)))^2 ;
end

%=====================================================================================

%                    Plot of parameter values vs Iteration numbers
%=====================================================================================

figure(1)
plot(thetait') ;
title('Parameter values vs iteration number','fontsize',14) ;
xlabel('Iteration Number ','fontsize',14)
ylabel('Parameter values','fontsize',14)

%=====================================================================================

clear H Covmat Covmattheta Kstar Lo Gphi Glam1 Glam2 Gtheta ...
    K Kmat3d M Nm Nit lamhat dof Psihat thetan Covmattheta sigeq2 Covmateps L1 L2 L3 ...
    nvCovmateps11 nvCovmateps12 nvCovmateps21 nvCovmateps22
toc

end