function [Kmat3d] = Kmat(nf)

Kmat3d = zeros(nf,nf,nf+1) ;


for i = 1:1:nf-1
    
    Kdum = zeros(nf,nf) ;
    Kdum(i,i) = 1 ;
    Kdum(i,i+1) = -1 ;
    Kdum(i+1,i+1) = 1 ;
    Kdum(i+1,i) = -1 ;
    Kmat3d(:,:,i+1) = Kdum ;
    

end

Kgr = zeros(nf,nf) ;
Kgr(nf,nf) = 1 ;
Kmat3d(:,:,nf+1) = Kgr ;


end