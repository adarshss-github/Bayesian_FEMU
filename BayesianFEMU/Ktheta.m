function [K] = Ktheta(theta)

K = zeros(12,12) ;

for i = 1:1:11
    
    Kdum = zeros(12,12) ;
    Kdum(i,i) = 10^6 ;
    Kdum(i,i+1) = -10^6 ;
    Kdum(i+1,i+1) = 10^6 ;
    Kdum(i+1,i) = -10^6 ;
    K = theta(i)*Kdum + K ;

end

Kgr = zeros(12,12) ;
Kgr(12,12) = theta(12)*10^6 ;
K = Kgr + K ;


end