for i = 1:1:65
    Sigeps(i,i) = ((1/100)*data1(i))^2 ;
end

data1n = diag(Sigeps) ;
Covmateps = Sigeps ;

lamhat = data1n(1:5) ;
psihat = data1n(6:end) ;
L0 = [eye(60,60)] ;
Covmattheta = (1/100^2)*eye(12,12) ;
sig2 = 1 ;