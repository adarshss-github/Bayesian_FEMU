function Theta = MH_MCMC_GAUSSPROP(postpdf,proppdfSigma,theta0,Ntheta,burnVal)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%Theta = MH_MCMC_GAUSSPROP(postpdf,proppdfMu,proppdfSigma,theta0,Ntheta,burnVal) ;
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Function description:
%------------------------
%Generates samples from a multivariable posterior PDF corresponding by using the
%Metropolis-Hastings Markov Chain Monte Carlo algorithm and a Gaussian
%proposal PDF
%
%%Input Arguments:
%----------------
%1) postpdf: The function handle of the posterior PDF
%2) proppdfSigma : Covariance matrix of the proposal PDF
%3) thet0: Column vector containing the initial value of the parameter
%4) Ntheta: Number of iterations to be run
%5) burnVal: Number of burn-in values to be discarded
%
%Output Arguments:
%-----------------
%1) Theta: Samples generated from the posterior PDF
%
%Example:MH_MCMC_GAUSSPROP(pdf,Covmatprop,[0.9;0.9;0.9;0.9],50000,100) ; ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%

%===========================================================================================================================================================
%***********************************************************************************************************************************************************

[mtheta] = length(theta0) ;

Theta(:,1) = theta0 ;

for i = 2:1:Ntheta

    U = rand(1) ;
    thetadum = mvnrnd(Theta(:,i-1)',proppdfSigma) ;
    r = ( mvnpdf(Theta(:,i-1)',thetadum,proppdfSigma)*postpdf(thetadum) )...
        /(  mvnpdf(thetadum',Theta(:,i-1),proppdfSigma)*postpdf(Theta(:,i-1))) ;

    if U<r

        Theta(:,i) = thetadum ;

    else

        Theta(:,i) = Theta(:,i-1) ;

    end

end

Theta(:,1:burnVal) = [] ;

end