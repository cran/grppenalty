cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Fortran code for grouped concave penalties, including the L1 
!     and L2 GSCAD and GMCP. 
!     The L2 GMCP and GSCAD includes the group Lasso as a special case.
!     This is both for linear and logistic models.
!     Version 1.0
!     June 6, 2012
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL I FUNCTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     index of start & end  for each group
c     dfgrp: x1, x2,    x3, ...
c     sinx : 1 , x1+1,  x1+x2+1,
c     einx : x1, x1+x2, x1+x2+x3
c***********************************************************************
      subroutine fseinx(sinx,einx,dfgrp,ngrp)
      integer ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),i
      sinx(1)=1
      einx(1)=dfgrp(1)
      do 00004 i=2,ngrp
         sinx(i)=sinx(i-1)+dfgrp(i-1)
         einx(i)=einx(i-1)+dfgrp(i)
00004 continue
      end

C     perform column-wise standardization
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine standard(as,sz,z,n,p)
      integer n,p
      double precision as(p),sz(n,p),z(n,p)
      integer j
      do 10000 j=1,p
         as(j)=sqrt( dble(n)/dot_product(z(:,j),z(:,j)) )
         sz(:,j)=as(j)*z(:,j)
10000 continue
      end

!     perform group-wise standardization
c***********************************************************************
      subroutine grpstd(cholinv,szj,zj,n,dj)
      integer n,dj
      double precision cholinv(dj,dj),szj(n,dj),zj(n,dj)
!     local vars
      integer info
      call dsyrk("U","T",dj,n,1.d0,zj,n,0.d0,cholinv,dj)
      cholinv=cholinv/dble(n)
!     zjtzj=R'R
      call dpotrf("U",dj,cholinv,dj,info)
!     get inverse of R
      call dtrtri("U","N",dj,cholinv,dj,info)
      szj(:,:)=zj(:,:)
      call dtrmm("R","U","N","N",n,dj,1.d0,cholinv,dj,szj,n)
      end

c     convergence check, using L2 norm
c***********************************************************************
      subroutine converge1(tag,newx,oldx,lng,epsilon)
      integer tag,lng
      double precision newx(lng),oldx(lng), epsilon
c     local var
      integer i
      double precision nnorm(lng),onorm(lng),nsum,osum
      tag=1
      do 00003 i=1,lng
         nnorm(i)=newx(i)*newx(i)
         onorm(i)=oldx(i)*oldx(i)
00003 continue
      nsum=sqrt(sum(nnorm(:)))
      osum=sqrt(sum(onorm(:)))
      if ( abs(nsum-osum)/(osum+0.01) .gt. epsilon) tag=0
      end
c     convergence check, using L2 norm of both alpha and beta
c***********************************************************************
      subroutine converge2(tag,newx,oldx,lgx,newz,oldz,lgz,epsilon)
      integer tag,lgx,lgz
      double precision newx(lgx),oldx(lgx),newz(lgz),oldz(lgz),epsilon
c     local var
      integer i
      double precision nx(lgx),ox(lgx),nz(lgz),oz(lgz),sn,so
      tag=1
      do 00005 i=1,lgx
         nx(i)=newx(i)*newx(i)
         ox(i)=oldx(i)*oldx(i)
00005 continue
      do 00006 i=1,lgz
         nz(i)=newz(i)*newz(i)
         oz(i)=oldz(i)*oldz(i)
00006 continue
      sn=sqrt( sum(nx(:)) + sum(nz(:)) )
      so=sqrt( sum(ox(:)) + sum(oz(:)) )
      if ( abs(sn-so)/(so+0.01) .gt. epsilon) tag=0
      end

c     soft threshold
c***********************************************************************
      subroutine soft(out,m,lambda)
      double precision out,m,lambda
      if (lambda .le. 0.d0) lambda=0.d0
      if (abs(m) .le. lambda) then
         out=0.d0
      else if (m .gt. lambda) then
         out=m-lambda
      else if (m .lt. -lambda) then
         out=m+lambda
      endif
      end

!     Calculate the AUC of ROC given cvpy,cvpi
!***********************************************************************
      subroutine aucroc(auc,y,pi,n)
      integer n
      double precision auc,y(n),pi(n)
!     local vars
      integer i,j,n0,grv(n),ltv(n),eqv(n),tmp,l
      double precision rankv(n),u0
      n0=0
!     get the group size, n0, n1
      do 1000 i=1,n
         if (y(i) .eq. 0.d0) n0=n0+1
 1000 continue
!     get the ranks of data
      do 1001 i=1,n
         grv(i)=0
         ltv(i)=0
         eqv(i)=0
         do 1002 j=1,n
            if ( pi(i) .lt. pi(j) ) then
               grv(i)=grv(i)+1
            else if ( pi(i) .eq. pi(j) ) then
               eqv(i)=eqv(i)+1
            else
               ltv(i)=ltv(i)+1
            end if
 1002    continue
         grv(i)=n-grv(i)
         ltv(i)=1+ltv(i)
         if ( eqv(i) .eq. 1) then
            rankv(i)=dble( grv(i) )
         else
            tmp=0
            do 1003 l=ltv(i),grv(i)
               tmp=tmp+l
 1003       continue
            rankv(i)=dble(tmp)/dble(eqv(i))
         end if
 1001 continue
!     calculate u statistics
      u0= - dble( n0*(n0+1)/2 )
      do 1004 i=1,n
         if ( y(i) .eq. 0.d0) u0=u0+rankv(i)
 1004 continue
      auc=dble(u0)/dble( n0*(n-n0) )
      if (auc .lt. 0.5) auc=1.d0-auc
      end






cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL II FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     L1 group
c     Calculate lambda_max, byproduct:alpha,beta,xinvxtx,r
c***********************************************************************
      subroutine l1maxga(lmdamax,alpha,beta,xinvxtx,r,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision lmdamax,alpha(q),beta(p),xinvxtx(n,q),
     +     r(n),y(n),x(n,q),z(n,p)
c     local vars
      integer info,i,j
      double precision xtx(q,q),px(n,n),qx(n,n),tmp(p),tmpgrp(ngrp)
      call dsyrk("U","T",q,n,1.d0,x,n,0.d0,xtx,q)
      call dpotrf("U",q,xtx,q,info)
      call dpotri("U",q,xtx,q,info)
c     xinvxtx=x(x'x)^{-1}
      call dsymm("R","U",n,q,1.d0,xtx,q,x,n,0.d0,xinvxtx,n)
c     alpha=(x'x)^{-1}x'y, beta=0
      call dgemv("T",n,q,1.d0,xinvxtx,n,y,1,0.d0,alpha,1)
      beta(:)=0.d0
c     calculate residual r= ( I - x(x'x)^{-1}x')) Y 
c     px=x(x'x)^{-1}x'
      call dgemm("N","T",n,n,q,1.d0,xinvxtx,n,x,n,0.d0,px,n)
c     qx=I-px
      do 1200 j=1,n
         do 1201 i=1,j
            if (i .eq. j) then
               qx(i,j)=1.d0-px(i,j)
            else
               qx(i,j)=-0.5*(px(i,j)+px(j,i))
            endif
 1201    continue
 1200 continue
c     r=qx*y
      call dsymv("U",n,1.d0,qx,n,y,1,0.d0,r,1)
c     tmp=z'_{ij}r
      do 1202 i=1,p
         tmp(i)=abs( (dot_product(z(:,i),r))/dble(n) )
 1202 continue
      do 1203 j=1,ngrp
         tmpgrp(j)=maxval(tmp(sinx(j):einx(j)))/dble(dfgrp(j))
 1203 continue
      lmdamax=maxval(tmpgrp)
      end

c     Lasso solution as initial for l1MCP and l1scad
c***********************************************************************
      subroutine l1iniga(alpha,beta,r,xinvxtx,lmda,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),r(n),xinvxtx(n,q),lmda,
     +     y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),newlmda,
     +     m,numor
c     begin of iteration
      count=0
 1204 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1205 j=1,ngrp
         newlmda= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            m=beta(idx) + dot_product(z(:,idx),r)/dble(n)
            call soft(numor,m,newlmda)
            beta(idx)=numor
c     .... update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1205 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 1209 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
 1209 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1210
      if (count .ge. maxit) call rexit("Lasso solution diverges! \n")
      if (count .lt. maxit) go to 1204
 1210 continue
c     end of iteration loop
      end            

c     MCP solution given lambda, kappa
c***********************************************************************
      subroutine l1lmga(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),thresh,m,sumabs,newlmda,
     +     numor
c     begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         thresh= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            m=beta(idx) + dot_product(z(:,idx),r)/dble(n)
            sumabs=sum(abs(beta(sinx(j):einx(j))))
            newlmda=thresh - ka*( sumabs-abs(beta(idx)) )
            if (abs(m) .lt. newlmda/ka) then
               call soft(numor,m,newlmda)
               beta(idx)=numor/(1.d0-ka)
            else
               beta(idx)=m
            end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1001 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end

c     SCAD solution given lambda, kappa
c***********************************************************************
      subroutine l1scadga(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),thresh,m,sumabs,newlmda,
     +     numor,gamma,tuna
c     begin of iteration
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     call dblepr("beta",-1,beta,p)
c     .. update of beta
      do 1001 j=1,ngrp
         thresh= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            m=beta(idx) + dot_product(z(:,idx),r)/dble(n)
            sumabs=sum(abs(beta(sinx(j):einx(j))))
            newlmda= sumabs-abs(beta(idx)) 
            if (abs(m) .le. (2.d0*thresh-newlmda) ) then
               call soft(numor,m,thresh)
               beta(idx)=numor
            else if (abs(m) .ge. gamma*thresh- newlmda) then
               beta(idx)=m
            else
               tuna=(gamma*thresh-newlmda)/(gamma-1)
               call soft(numor,m,tuna)
               beta(idx)=numor*(gamma-1)/(gamma-2)
            end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1001 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end

c     L2 group 
c***********************************************************************
c     Calculate lambda_max, byproduct:alpha,beta,xinvxtx,r
c***********************************************************************
      subroutine l2maxga(lmdamax,alpha,beta,xinvxtx,r,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision lmdamax,alpha(q),beta(p),xinvxtx(n,q),
     +     r(n),y(n),x(n,q),z(n,p)
c     local vars
      integer info,i,j
      double precision xtx(q,q),px(n,n),qx(n,n),tmp(p),tmpgrp(ngrp)
      call dsyrk("U","T",q,n,1.d0,x,n,0.d0,xtx,q)
      call dpotrf("U",q,xtx,q,info)
      call dpotri("U",q,xtx,q,info)
c     xinvxtx=x(x'x)^{-1}
      call dsymm("R","U",n,q,1.d0,xtx,q,x,n,0.d0,xinvxtx,n)
c     alpha=(x'x)^{-1}x'y, beta=0
      call dgemv("T",n,q,1.d0,xinvxtx,n,y,1,0.d0,alpha,1)
      beta(:)=0.d0
c     calculate residual r= ( I - x(x'x)^{-1}x')) Y 
c     px=x(x'x)^{-1}x'
      call dgemm("N","T",n,n,q,1.d0,xinvxtx,n,x,n,0.d0,px,n)
c     qx=I-px
      do 1200 j=1,n
         do 1201 i=1,j
            if (i .eq. j) then
               qx(i,j)=1.d0-px(i,j)
            else
               qx(i,j)=-0.5*(px(i,j)+px(j,i))
            endif
 1201    continue
 1200 continue
c     r=qx*y
      call dsymv("U",n,1.d0,qx,n,y,1,0.d0,r,1)
c     tmp=z'_{ij}r
      do 1202 i=1,p
         tmp(i)=((dot_product(z(:,i),r))/dble(n))**2
 1202 continue
      do 1203 j=1,ngrp
         tmpgrp(j)=sqrt( sum(tmp(sinx(j):einx(j)))/dble(dfgrp(j)) )
 1203 continue
      lmdamax=maxval(tmpgrp)
      end

c     Lasso solution as initial for MCP
c***********************************************************************
      subroutine l2iniga(alpha,beta,r,xinvxtx,lmda,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),r(n),xinvxtx(n,q),lmda,
     +     y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),newlmda,
     +     tmp(p),tmpsq(p),m
c     begin of iteration
      count=0
 1204 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1205 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)= beta(idx) + dot_product(z(:,idx),r)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .le. newlmda) then 
            beta(sinx(j):einx(j))=0.d0
         else 
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=(1.d0-newlmda/m)*tmp(idx)
 1207       continue
         end if
c     .... update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1208 idx=sinx(j),einx(j)
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1208    continue
 1205 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 1209 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
 1209 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1210
      if (count .ge. maxit) call rexit("Lasso solution diverges! \n")
      if (count .lt. maxit) go to 1204
 1210 continue
c     end of iteration loop
      end            

c     MCP solution given lambda, kappa
c***********************************************************************
      subroutine l2lmga(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),newlmda,
     +     tmp(p),tmpsq(p),m
c     begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)=beta(idx) + dot_product(z(:,idx),r)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .le. newlmda) then 
            beta(sinx(j):einx(j))=0.d0
         else if (m .ge. newlmda/ka) then 
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=tmp(idx)
 1207       continue            
         else 
            do 1208 idx=sinx(j),einx(j)
               beta(idx)=(1.d0-newlmda/m)*tmp(idx)/(1.d0-ka)
 1208       continue
         end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1209 idx=sinx(j),einx(j)
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1209    continue
 1001 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end

c     SCAD solution given lambda, kappa
c***********************************************************************
      subroutine l2scadga(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,idx,tag
      double precision alphaold(q),betaold(p),newlmda,gamma,
     +     tmp(p),tmpsq(p),m
      gamma=1.d0/ka
c     begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)=beta(idx) + dot_product(z(:,idx),r)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .le. 2.d0*newlmda) then 
            if (m .le. newlmda) then
               beta(sinx(j):einx(j))=0.d0
            else
               do 9000 idx=sinx(j),einx(j)
                  beta(idx)=(1.d0-newlmda/m)*tmp(idx)
 9000          continue
            end if     
         else if (m .ge. newlmda*gamma) then 
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=tmp(idx)
 1207       continue            
         else 
!     if ((m .gt. 2.d0*newlmda) .and. (m .lt. newlmda*ka) ) then
            if (m .le. (gamma*newlmda/(gamma-1.d0)) ) then
               beta(sinx(j):einx(j))=0.d0
            else 
               do 90001 idx=sinx(j),einx(j)
                  beta(idx)=(gamma-1.d0)/(gamma-2.d0)*
     +                 (1.d0-gamma*newlmda/m/(gamma-1))*tmp(idx)
90001          continue
            end if 
         end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1209 idx=sinx(j),einx(j)
            r(:)=r(:)+z(:,idx)*( betaold(idx)-beta(idx) )
 1209    continue
 1001 continue
c     .. update of alpha, alpha=alpha+(x'x)^{-1}x'r
      call dgemv("T",n,q,1.d0,xinvxtx,n,r,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:)+x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end






cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL III FUNCTIONS: solution surface
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     L1+GMCP
c***********************************************************************
      subroutine l1gmcpga(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1lmga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call f1msevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c***********************************************************************
      subroutine bfl1gmcpga(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),y(n),x(n,q),z(n,p),kas(nka),
     +     minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1lmga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call f1msev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c     L1+GSCAD
c***********************************************************************
      subroutine l1gscadga(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1scadga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call f1msevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c***********************************************************************
      subroutine bfl1gscadga(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),y(n),x(n,q),z(n,p),kas(nka),
     +     minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
c     call intpr("nka",-1,j,1)
               call l1scadga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call f1msev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c     L2+GMCP
c***********************************************************************
      subroutine l2gmcpga(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2lmga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmsevabic(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c     get solution path, df, convex index only
c***********************************************************************
      subroutine bfl2gmcpga(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),y(n),x(n,q),z(n,p),kas(nka),
     +     minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2lmga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue
      do 10 i=1,(nka*nlmda)
         call fmsev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c     L2+GSCAD
c***********************************************************************
      subroutine l2gscadga(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2scadga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call fmsevabic(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c     get solution path, df, convex index only
c***********************************************************************
      subroutine bfl2gscadga(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),y(n),x(n,q),z(n,p),kas(nka),
     +     minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2scadga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue
      do 10 i=1,(nka*nlmda)
         call fmsev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end


c     LEVEL IV function
c***********************************************************************
c     Tuning parameter selection part
c***********************************************************************

c     L1+GMCP
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l1cvga(pmse,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pmse(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpr(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call l1maxga(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l1iniga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l1lmga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call f1msev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the pmses for the prediction set
      do 1008 i=1,(nka*nlmda)
         do 1009 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1009    continue
         cvpr(:)=cvpy(:)-cvpeta(:)
         tmp=0.d0
         do 1010 j=1,cvpn
            tmp=cvpr(j)**2+tmp
 1010    continue
         pmse(i)=tmp/dble(cvpn)
 1008 continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl1gmcpga(out,pmse,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pmse(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpmse(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompmse(nka*nlmda)
c     solution path of original dataset
      call bfl1gmcpga(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pmse for cv datasets
         call l1cvga(cvpmse(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pmse, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pmse(i)=sum(cvpmse(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompmse(use) =pmse(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get minimum
      loc=minloc( ompmse(1:use), use )
      out(1)=ompmse(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L1+GSCAD
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l1cvscadga(pmse,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pmse(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpr(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call l1maxga(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l1iniga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l1scadga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call f1msev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the pmses for the prediction set
      do 1008 i=1,(nka*nlmda)
         do 1009 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1009    continue
         cvpr(:)=cvpy(:)-cvpeta(:)
         tmp=0.d0
         do 1010 j=1,cvpn
            tmp=cvpr(j)**2+tmp
 1010    continue
         pmse(i)=tmp/dble(cvpn)
 1008 continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl1gscadga(out,pmse,olmda,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pmse(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpmse(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompmse(nka*nlmda)
c     solution path of original dataset
      call bfl1gscadga(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pmse for cv datasets
         call l1cvscadga(cvpmse(:,cv),cvfullm(:,cv),cvcvxm(:,cv),
     +        lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pmse, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pmse(i)=sum(cvpmse(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompmse(use) =pmse(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get minimum
      loc=minloc( ompmse(1:use), use )
      out(1)=ompmse(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L2+MCP
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l2cvga(pmse,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pmse(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvcholmat(p,p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpr(cvpn),tmp
!     groupwise standardization
      do 1000 j=1,ngrp
         call grpstd(cvcholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        cvtsz(:,sinx(j):einx(j)),cvtz(:,sinx(j):einx(j)),cvtn,
     +        dfgrp(j))
 1000 continue
c     solution path along lambdas, initial values for along kas
      call l2maxga(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l2iniga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l2lmga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fmsev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cvcholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        cvbs(sinx(j):einx(j),:),dfgrp(j))
10005 continue
c     get the pmses for the prediction set
      do 1008 i=1,(nka*nlmda)
         do 1009 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1009    continue
         cvpr(:)=cvpy(:)-cvpeta(:)
         tmp=0.d0
         do 1010 j=1,cvpn
            tmp=cvpr(j)**2+tmp
 1010    continue
         pmse(i)=tmp/dble(cvpn)
 1008 continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl2gmcpga(out,pmse,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pmse(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpmse(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompmse(nka*nlmda)
c     solution path of original dataset
      call bfl2gmcpga(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pmse for cv datasets
         call l2cvga(cvpmse(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pmse, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pmse(i)=sum(cvpmse(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompmse(use) =pmse(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get minimum
      loc=minloc( ompmse(1:use), use )
      out(1)=ompmse(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L2+SCAD
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l2cvscadga(pmse,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pmse(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvcholmat(p,p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpr(cvpn),tmp
!     groupwise standardization
      do 1000 j=1,ngrp
         call grpstd(cvcholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        cvtsz(:,sinx(j):einx(j)),cvtz(:,sinx(j):einx(j)),cvtn,
     +        dfgrp(j))
 1000 continue
c     solution path along lambdas, initial values for along kas
      call l2maxga(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l2iniga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l2scadga(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call fmsev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cvcholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        cvbs(sinx(j):einx(j),:),dfgrp(j))
10005 continue
c     get the pmses for the prediction set
      do 1008 i=1,(nka*nlmda)
         do 1009 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1009    continue
         cvpr(:)=cvpy(:)-cvpeta(:)
         tmp=0.d0
         do 1010 j=1,cvpn
            tmp=cvpr(j)**2+tmp
 1010    continue
         pmse(i)=tmp/dble(cvpn)
 1008 continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl2gscadga(out,pmse,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pmse(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpmse(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompmse(nka*nlmda)
c     solution path of original dataset
      call bfl2gscadga(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pmse for cv datasets
         call l2cvscadga(cvpmse(:,cv),cvfullm(:,cv),cvcvxm(:,cv),
     +        lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pmse, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pmse(i)=sum(cvpmse(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompmse(use) =pmse(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get minimum
      loc=minloc( ompmse(1:use), use )
      out(1)=ompmse(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This section is using the solution surface along lambda 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     L2+ GMCP
c     get solution path, df, convex index only
c***********************************************************************
      subroutine bfl2gmcpgalm(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),y(n),x(n,q),z(n,p),kas(nka),
     +     minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j,idx
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia0(q),inib0(p),rmat0(n),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max, compute lambdas
      call l2maxga(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      inia0(:)=alpha(:)
      inib0(:)=beta(:)
      rmat0(:)=r(:)
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
10001 continue
c     compute the solution path along lambdas, for a given kappa
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         alpha(:)=inia0(:)
         beta(:) =inib0(:)
         r(:)    =rmat0(:)
         ocoef(1:q,1)=alpha(:)
         ocoef((q+1):(q+p),1)=beta(:)
         do 1002 i=2,nlmda
            call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +           ngrp,dfgrp,sinx,einx,epsilon,maxit)
            ocoef(1:q,i)=alpha(:)
            ocoef((q+1):(q+p),i)=beta(:)
 1002    continue
      else
         j=1
         i=1
         idx=(i-1)*nka+j
         olmdas(idx)=lmdas(i)
         okas(idx)  =kas(j)
         alpha(:)=inia0(:)
         beta(:) =inib0(:)
         r(:)    =rmat0(:)
         ocoef(1:q,idx)=alpha(:)
         ocoef((q+1):(q+p),idx)=beta(:)
         do 1003 i=2,nlmda
            idx=(i-1)*nka+j
            olmdas(idx)=lmdas(i)
            okas(idx)  =kas(j)
            call l2iniga(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +           ngrp,dfgrp,sinx,einx,epsilon,maxit)
            ocoef(1:q,idx)=alpha(:)
            ocoef((q+1):(q+p),idx)=beta(:)
1003     continue   
         do 1004 j=2,nka
            i=1
            idx=(i-1)*nka+j
            olmdas(idx)=lmdas(i)
            okas(idx)  =kas(j)
            alpha(:)=inia0(:)
            beta(:) =inib0(:)
            r(:)    =rmat0(:)
            ocoef(1:q,idx)=alpha(:)
            ocoef((q+1):(q+p),idx)=beta(:)
            do 1005 i=2,nlmda
               idx=(i-1)*nka+j
               olmdas(idx)=lmdas(i)
               okas(idx)  =kas(j)
               call l2lmga(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q,idx)=alpha(:)
               ocoef((q+1):(q+p),idx)=beta(:)
1005        continue
1004     continue
      endif
c     compute model size,eigenvalue
      do 10 i=1,(nka*nlmda)
         call fmsev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     fuctions used to compute AIC, BIC etc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     For L2 group penalty
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Given solution, kappa,lmda,
c     determine model size,aic, bic, objective value
C     Note: convex dianosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fmsevabic(ms,mszg,evgx,aic,bic,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,mszg,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision aic,bic,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer j,idx,i,info,lwork
      double precision alpha(q),beta(p),xnz(n,q+p),anb(q+p),eta(n),
     +     sumsq,logl2
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups, individuals
      xnz(:,1:q)=x(:,:)
      anb(1:q)=alpha(:)      
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
               xnz(:,ms)=z(:,idx)
               anb(ms)=beta(idx)
 1001       continue
         endif
 1000 continue
!     .. sum of square, aic, bic
      call dgemv("N",n,ms,1.d0,xnz(:,1:ms),n,anb(1:ms),1,0.d0,eta,1)
      sumsq=0.d0
      do 11001 i=1,n
         sumsq=sumsq+(y(i)-eta(i))**2
11001 continue
      logl2=dble(n)*(log(2* 3.141593)+log(sumsq)-log(dble(n))+1)
      aic=logl2 + ms*2.d0
      bic=logl2 + ms*log(dble(n))
      end

c     Given solution, kappa,lmda,
c     determine model size
C     convex diagnosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fmsev(ms,evgx,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer mszg,j,idx
      double precision alpha(q),beta(p)
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups, individuals
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
 1001       continue
         endif
 1000 continue
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fmsev2(full,evgx,alpha,beta,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer full,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision alpha(q),beta(p),y(n),x(n,q),z(n,p)
!     .. local arguments
      integer mszg,ms,j,idx
c     model size
c     .. non-zero groups, individuals
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
 1001       continue
         endif
 1000 continue
      full=1
      if (ms .le. n) full=0
      end

c     For L1 group penalty
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Given solution, kappa,lmda,
c     determine model size,aic, bic, objective value
C     Note: convex dianosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f1msevab(ms,mszg,evgx,aic,bic,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,mszg,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision aic,bic,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer j,idx,i,info,lwork
      double precision alpha(q),beta(p),xnz(n,q+p),anb(q+p),eta(n),
     +     sumsq,logl2
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups,
      mszg=0
      do 11000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
         endif
11000 continue
c     .. non-zero individuals
      xnz(:,1:q)=x(:,:)
      anb(1:q)=alpha(:)      
      ms=q
      do 1000 j=1,p
         if ( beta(j) .ne. 0.d0 ) then
            ms=ms+1
            xnz(:,ms)=z(:,j)
            anb(ms)=beta(j)
         endif
 1000 continue
!     .. sum of square, aic, bic
      call dgemv("N",n,ms,1.d0,xnz(:,1:ms),n,anb(1:ms),1,0.d0,eta,1)
      sumsq=0.d0
      do 11001 i=1,n
         sumsq=sumsq+(y(i)-eta(i))**2
11001 continue
      logl2=dble(n)*(log(2* 3.141593)+log(sumsq)-log(dble(n))+1)
      aic=logl2 + ms*2.d0
      bic=logl2 + ms*log(dble(n))
      end

c     Given solution, kappa,lmda,
c     determine model size
C     convex diagnosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f1msev(ms,evgx,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer j
      double precision alpha(q),beta(p)
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero individuals
      ms=q
      do 1000 j=1,p
         if (  beta(j) .ne. 0.d0) then
            ms=ms+1
         endif
 1000 continue
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f1msev2(full,evgx,alpha,beta,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer full,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision alpha(q),beta(p),y(n),x(n,q),z(n,p)
!     .. local arguments
      integer ms,j
c     model size
c     .. non-zero individuals
      ms=q
      do 1000 j=1,p
         if ( beta(j) .ne. 0.d0) then
            ms=ms+1
         endif
 1000 continue
      full=1
      if (ms .le. n) full=0
c     call intpr("n",-1,n,1)
c     call intpr("ms",-1,ms,1)
c     call intpr("full",-1,full,1)
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Logistic regression
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     LEVEL II FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     L1 group
c     Calculate lambda_max, byproduct:alpha,beta,xinvxtx,r
c***********************************************************************
      subroutine l1maxbi(lmdamax,alpha,beta,xinvxtx,r,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision lmdamax,alpha(q),beta(p),xinvxtx(n,q),
     +     r(n),y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer info,count,i,j,tag
      double precision xtx(q,q),alphaold(q),pi(n),rpi(n),
     +     tmp(p),tmpgrp(ngrp)
c     .... xinvxtx=x(x'x)^{-1}
      call dsyrk("U","T",q,n,1.d0,x,n,0.d0,xtx,q)
      call dpotrf("U",q,xtx,q,info)
      call dpotri("U",q,xtx,q,info)
      call dsymm("R","U",n,q,1.d0,xtx,q,x,n,0.d0,xinvxtx,n)
!     Calculate initial value of alpha, beta
      beta(:)=0.d0
      alpha(:)=0.d0
      count=0
c     .. start of iteration
 1204 continue
      alphaold(:)=alpha(:)
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,r,1)
      do 1205 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         rpi(i)=y(i)-pi(i)
 1205 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
      count=count+1
      call converge1(tag,alpha,alphaold,q,epsilon)
      if (tag .eq. 1) go to 1206
      if (count .ge. maxit) call rexit("Diverge of alpha initials!\n")
      if (count .lt. maxit) go to 1204
c     .. end of iteration loop
 1206 continue
c     start of lambda_max
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,r,1)
      do 1207 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         rpi(i)=y(i)-pi(i)
 1207 continue
      do 1208 i=1,p
         tmp(i)=abs((dot_product(z(:,i),rpi))/dble(n))
 1208 continue
      do 1209 j=1,ngrp
         tmpgrp(j)= maxval(tmp(sinx(j):einx(j)))/dble(dfgrp(j))
 1209 continue
      lmdamax=maxval(tmpgrp)
      end


c     Lasso solution as initial for MCP
c***********************************************************************
      subroutine l1inibi(alpha,beta,r,xinvxtx,lmda,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),r(n),xinvxtx(n,q),lmda,
     +     y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),newlmda,pi(n),rpi(n),
     +     m,numor
c     begin of iteration
      count=0
 1204 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1205 j=1,ngrp
         newlmda= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            do 1213 i=1,n
               pi(i)=1/(1+exp(-r(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               rpi(i)=y(i)-pi(i)
 1213       continue
            m= beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            call soft(numor,m,newlmda)
            beta(idx)=4.d0*numor
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1205 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 1209 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
 1209 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1210
      if (count .ge. maxit) call rexit("Lasso solution diverges! \n")
      if (count .lt. maxit) go to 1204
 1210 continue
c     end of iteration loop
      end            


c     MCP solution given lambda, kappa
c***********************************************************************
      subroutine l1lmbi(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),thresh,pi(n),rpi(n),
     +     m,sumabs,newlmda,numor
c     begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         thresh= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            do 1213 i=1,n
               pi(i)=1/(1+exp(-r(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               rpi(i)=y(i)-pi(i)
 1213       continue
            m=beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            sumabs=sum(abs(beta(sinx(j):einx(j))))
            newlmda=thresh - (sumabs-abs(beta(idx)))*ka
            if (abs(m) .lt. 0.25*newlmda/ka) then
               call soft(numor,m,newlmda)
               beta(idx)=numor/(0.25-ka)
            else
               beta(idx)=4.d0*m
            end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1001 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end


c     SCAD solution given lambda, kappa
c***********************************************************************
      subroutine l1scadbi(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),thresh,pi(n),rpi(n),
     +     m,sumabs,newlmda,numor,gamma,tuna
c     begin of iteration
c     call dblepr("lmda",-1,lmda,1)
c     call dblepr("ka",-1,ka,1)
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         thresh= dble(dfgrp(j))*lmda
         do 1206 idx=sinx(j),einx(j)
            do 1213 i=1,n
               pi(i)=1/(1+exp(-r(i)))
               if (pi(i) .lt. 0.0001) pi(i)=0.d0
               if (pi(i) .gt. 0.9999) pi(i)=1.d0
               rpi(i)=y(i)-pi(i)
 1213       continue
            m=beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            sumabs=sum(abs(beta(sinx(j):einx(j))))
            newlmda=sumabs-abs(beta(idx))
            if (abs(m) .le. (1.25*thresh-0.25*newlmda) ) then
               call soft(numor,m,thresh)
               beta(idx)=4.d0*numor
            else if (abs(m) .ge. 0.25*(thresh*gamma-newlmda) ) then
               beta(idx)=4.d0*m
            else 
               tuna=(thresh*gamma-newlmda)/(gamma-1.d0)
               call soft(numor,m,tuna)
               beta(idx)=4.d0*(gamma-1.d0)*numor/(gamma-5.d0)
            end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1206    continue
 1001 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end




c     L2 group
c     Calculate lambda_max, byproduct:alpha,beta,xinvxtx,r
c***********************************************************************
      subroutine l2maxbi(lmdamax,alpha,beta,xinvxtx,r,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision lmdamax,alpha(q),beta(p),xinvxtx(n,q),
     +     r(n),y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer info,count,i,j,tag
      double precision xtx(q,q),alphaold(q),pi(n),rpi(n),
     +     tmp(p),tmpgrp(ngrp)
c     .... xinvxtx=x(x'x)^{-1}
      call dsyrk("U","T",q,n,1.d0,x,n,0.d0,xtx,q)
      call dpotrf("U",q,xtx,q,info)
      call dpotri("U",q,xtx,q,info)
      call dsymm("R","U",n,q,1.d0,xtx,q,x,n,0.d0,xinvxtx,n)
!     Calculate initial value of alpha, beta
      beta(:)=0.d0
      alpha(:)=0.d0
      count=0
c     .. start of iteration
 1204 continue
      alphaold(:)=alpha(:)
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,r,1)
      do 1205 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         rpi(i)=y(i)-pi(i)
 1205 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
      count=count+1
      call converge1(tag,alpha,alphaold,q,epsilon)
      if (tag .eq. 1) go to 1206
      if (count .ge. maxit) call rexit("Diverge of alpha initials!\n")
      if (count .lt. maxit) go to 1204
c     .. end of iteration loop
 1206 continue
c     start of lambda_max
      call dgemv("N",n,q,1.d0,x,n,alpha,1,0.d0,r,1)
      do 1207 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         rpi(i)=y(i)-pi(i)
 1207 continue
      do 1208 i=1,p
         tmp(i)=((dot_product(z(:,i),rpi))/dble(n))**2
 1208 continue
      do 1209 j=1,ngrp
         tmpgrp(j)=sqrt( sum(tmp(sinx(j):einx(j)))/dble(dfgrp(j)) )
 1209 continue
      lmdamax=maxval(tmpgrp)
      end


c     Lasso solution as initial for MCP
c***********************************************************************
      subroutine l2inibi(alpha,beta,r,xinvxtx,lmda,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),r(n),xinvxtx(n,q),lmda,
     +     y(n),x(n,q),z(n,p),epsilon
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),newlmda,pi(n),rpi(n),
     +     tmp(p),tmpsq(p),m
c     begin of iteration
      count=0
 1204 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1205 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1213 i=1,n
            pi(i)=1/(1+exp(-r(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            rpi(i)=y(i)-pi(i)
 1213    continue
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)= beta(idx)/4.d0+dot_product(z(:,idx),rpi)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .le. newlmda) then 
            beta(sinx(j):einx(j))=0.d0
         else 
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=4.d0*(1.d0-newlmda/m)*tmp(idx)
 1207       continue
         end if
c     .... update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1208 idx=sinx(j),einx(j)
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1208    continue
 1205 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 1209 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
 1209 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1210
      if (count .ge. maxit) call rexit("Lasso solution diverges! \n")
      if (count .lt. maxit) go to 1204
 1210 continue
c     end of iteration loop
      end            


c     MCP solution given lambda, kappa
c***********************************************************************
      subroutine l2lmbi(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),newlmda,pi(n),rpi(n),
     +     tmp(p),tmpsq(p),m
c     begin of iteration
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1213 i=1,n
            pi(i)=1/(1+exp(-r(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            rpi(i)=y(i)-pi(i)
 1213    continue
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)=beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .ge. 0.25*newlmda/ka) then
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=4.d0*tmp(idx)
 1207       continue            
         else 
            if (m .le. newlmda) then 
               beta(sinx(j):einx(j))=0.d0            
            else 
               do 1208 idx=sinx(j),einx(j)
                  beta(idx)=(1.d0-newlmda/m)*tmp(idx)/(0.25-ka)
 1208          continue
            end if
         end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1209 idx=sinx(j),einx(j)
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1209    continue
 1001 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end




c     MCP solution given lambda, kappa
c***********************************************************************
      subroutine l2scadbi(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),newlmda,pi(n),rpi(n),
     +     tmp(p),tmpsq(p),m,gamma,tuna,fish
c     begin of iteration
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1213 i=1,n
            pi(i)=1/(1+exp(-r(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            rpi(i)=y(i)-pi(i)
 1213    continue
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)=beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .ge. 0.25*newlmda/ka) then
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=4.d0*tmp(idx)
 1207       continue            
         else if (m .le. 1.25*newlmda) then
            if (m .le. newlmda) then 
               beta(sinx(j):einx(j))=0.d0            
            else 
               do 1208 idx=sinx(j),einx(j)
                  beta(idx)=(1.d0-newlmda/m)*tmp(idx)*4.d0
 1208          continue
            end if
         else 
            tuna=gamma*newlmda/(gamma-1.d0)
            if (m .le. tuna) then
               beta(sinx(j):einx(j))=0.d0    
            else
               do 1308 idx=sinx(j),einx(j)
                  fish=4.d0*(gamma-1.d0)/(gamma-5.d0)
                  beta(idx)=(1.d0-tuna/m)*tmp(idx)*fish
 1308          continue
            end if
         end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1209 idx=sinx(j),einx(j)
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1209    continue
 1001 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end



c     SCAD solution given lambda, kappa
c***********************************************************************
      subroutine l2scadb(alpha,beta,r,xinvxtx,lmda,ka,y,x,z,n,q,p,ngrp,
     +     dfgrp,sinx,einx,epsilon,maxit)
      integer n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp),maxit
      double precision alpha(q),beta(p),lmda,ka,r(n),y(n),x(n,q),z(n,p),
     +     epsilon,xinvxtx(n,q)
c     local vars
      integer count,j,i,idx,tag
      double precision alphaold(q),betaold(p),newlmda,pi(n),rpi(n),
     +     tmp(p),tmpsq(p),m,gamma,tuna
c     begin of iteration
      call dblepr("lmda",-1,lmda,1)
      call dblepr("ka",-1,ka,1)
      gamma=1.d0/ka
      count=0
 1000 continue
      alphaold(:)=alpha(:)
      betaold(:)=beta(:)
c     .. update of beta
      do 1001 j=1,ngrp
         newlmda=sqrt( dble(dfgrp(j)) )*lmda
         do 1213 i=1,n
            pi(i)=1/(1+exp(-r(i)))
            if (pi(i) .lt. 0.0001) pi(i)=0.d0
            if (pi(i) .gt. 0.9999) pi(i)=1.d0
            rpi(i)=y(i)-pi(i)
 1213    continue
         do 1206 idx=sinx(j),einx(j)
            tmp(idx)=beta(idx)/4.d0 + dot_product(z(:,idx),rpi)/dble(n)
            tmpsq(idx)=tmp(idx)**2
 1206    continue
         m=sqrt( sum(tmpsq(sinx(j):einx(j))) )
         if (m .ge. 0.25*newlmda*gamma) then
            do 1207 idx=sinx(j),einx(j)
               beta(idx)=4.d0*tmp(idx)
 1207       continue            
         else if (m .le. 1.25*newlmga) then
            if (m .le. newlmda) then 
               beta(sinx(j):einx(j))=0.d0            
            else 
               do 1208 idx=sinx(j),einx(j)
                  beta(idx)=4.d0*(1.d0-newlmda/m)*tmp(idx)
 1208          continue
            end if
         else
            tuna=gamma*newlmda/(gamma-1.d0)
            if (m .le. tuna) then
               beta(sinx(j):einx(j))=0.d0       
            else
               do 9000 idx=sinx(j),einx(j)
                  beta(idx)=4.d0*(gamma-1.d0)*(1.d0-tuna/m)*tmp(idx)
     +                 /(gamma-5.0)
 9000          continue
            end if
         end if
c     ....  update residuals, r(s+1)=r(s)+z(b(s)-b(s+1))
         do 1209 idx=sinx(j),einx(j)
            r(:)=r(:) - z(:,idx)*( betaold(idx)-beta(idx) )
 1209    continue
 1001 continue
!     .. update of alpha
      do 1214 i=1,n
         pi(i)=1/(1+exp(-r(i)))
         if (pi(i) .lt. 0.0001) pi(i)=0.d0
         if (pi(i) .gt. 0.9999) pi(i)=1.d0
         rpi(i)=y(i)-pi(i)
 1214 continue
c     .... alpha^(s+1)=alpha^{s} + 4 (X'X)^{-1}X'(Y-pi)
      call dgemv("T",n,q,4.d0,xinvxtx,n,rpi,1,1.d0,alpha,1)
c     .... update  r(s+1)=r(s)+x(a(s)-a(s+1))
      do 12109 j=1,q
         r(:)=r(:) - x(:,j)*( alphaold(j) - alpha(j) )
12109 continue
c     check convergence
      count=count+1
      call converge2(tag,alpha,alphaold,q,beta,betaold,p,epsilon)
      if (tag .eq. 1) go to 1003
      if (count .ge. maxit) call rexit("MCP solution diverges!")
      if (count .lt. maxit) go to 1000
 1003 continue
c     end of iteration loop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LEVEL III FUNCTIONS: solution surface
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     L1+GMCP
c***********************************************************************
      subroutine l1gmcpbi(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1lmbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bi1msevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c***********************************************************************
      subroutine bfl1gmcpbi(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1lmbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bi1msev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c     L1+GSCAD
c***********************************************************************
      subroutine l1gscadbi(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1scadbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bi1msevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c***********************************************************************
      subroutine bfl1gscadbi(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,as(p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
c     Column-wise standardization of Z
      call standard(as,sz,z,n,p)
c     calculate lambda_max
      call l1maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l1inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l1scadbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bi1msev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 i=1,p
         ocoef(q+i,:)=as(i)*ocoef(q+i,:)
10005 continue
      end

c     L2+GMCP
c***********************************************************************
      subroutine l2gmcpbi(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2lmbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bimsevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c***********************************************************************
      subroutine bfl2gmcpbi(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2lmbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,sz,
     +              n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bimsev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c     L2+GSCAD
c***********************************************************************
      subroutine l2gscadbi(olmdas,okas,ocoef,oaic,obic,odf,odfg,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     odfg(nka*nlmda),oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),oaic(nka*nlmda),obic(nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2scadbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bimsevab(odf(i),odfg(i),oevidx(i),oaic(i),obic(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end

c***********************************************************************
      subroutine bfl2gscadbi(olmdas,okas,ocoef,odf,oevidx,
     +     y,x,z,n,q,p,ngrp,dfgrp,
     +     nka,kas,nlmda,minlmda,epsilon,maxit)
      integer n,q,p,ngrp,nka,nlmda,maxit,dfgrp(ngrp),odf(nka*nlmda),
     +     oevidx(nka*nlmda)
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),
     +     y(n),x(n,q),z(n,p),kas(nka),minlmda,epsilon
c     local vars
      integer sinx(ngrp),einx(ngrp),i,j
      double precision unitka,cholmat(p,p),sz(n,p),
     +     lmdas(nlmda),alpha(q),beta(p),xinvxtx(n,q),r(n),
     +     inia(q,nlmda),inib(p,nlmda),rmat(n,nlmda),unitlmda
c     calculate sinx, einx  for later use
      call fseinx(sinx,einx,dfgrp,ngrp)
c     calculate kappas
!     if (nka .eq. 1) then
!     kas(1)=0.d0
!     else
!     unitka=maxka/dble(nka-1)
!     do 10000 i=1,nka
!     kas(i)=dble(i-1)*unitka
!     10000    continue
!     endif
!     perform groupwise standardization
      do 00005 j=1,ngrp
         call grpstd(cholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        sz(:,sinx(j):einx(j)),z(:,sinx(j):einx(j)),n,dfgrp(j))
00005 continue
c     calculate lambda_max
      call l2maxbi(lmdas(1),alpha,beta,xinvxtx,r,y,x,sz,n,q,p,
     +     ngrp,dfgrp,sinx,einx,epsilon,maxit)
      inia(:,1)=alpha(:)
      inib(:,1)=beta(:)
      rmat(:,1)=r(:)
c     lasso solution as initial for MCP
      unitlmda=log(minlmda)/dble(nlmda-1)
      do 10001 i=2,nlmda
         lmdas(i)=lmdas(1)*exp(unitlmda*dble(i-1))
         call l2inibi(alpha,beta,r,xinvxtx,lmdas(i),y,x,sz,n,q,p,
     +        ngrp,dfgrp,sinx,einx,epsilon,maxit)
         inia(:,i)=alpha(:)
         inib(:,i)=beta(:)
         rmat(:,i)=r(:)
10001 continue
c     MCP
      if (nka .eq. 1) then
         olmdas(:)=lmdas(:)
         okas(1:nlmda)=kas(1)
         ocoef(1:q,:)=inia(:,:)
         ocoef((q+1):(q+p),:)=inib(:,:)
      else
         do 1003 i=1,nlmda
            olmdas(((i-1)*nka+1):(i*nka))=lmdas(i)
            okas(((i-1)*nka+1):(i*nka))=kas(:)
 1003    continue
!     ... i=1 case
         i=1
         alpha(:)=inia(:,i)
         beta(:) =inib(:,i)
         r(:)  =rmat(:,i)
         do 1005 j=1,nka
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1005    continue
         do 1004 i=2,nlmda
            alpha(:)=inia(:,i)
            beta(:) =inib(:,i)
            r(:)  =rmat(:,i)
c     ...... j=1 case
            j=1
            ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
            ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
            do 1006 j=2,nka
               call l2scadbi(alpha,beta,r,xinvxtx,lmdas(i),kas(j),y,x,
     +              sz,n,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
               ocoef(1:q        ,(i-1)*nka+j)=alpha(:)
               ocoef((q+1):(q+p),(i-1)*nka+j)=beta(:)
 1006       continue
 1004    continue
      endif      
c     compute model size,eigenvalue,aic,bic, objective function
      do 10 i=1,(nka*nlmda)
         call bimsev(odf(i),oevidx(i),
     +        ocoef(:,i),y,x,sz,n,q,p,ngrp,dfgrp,sinx,einx,
     +        okas(i),olmdas(i))
 10   continue
!     change the penalized coefficients back 
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        ocoef((q+sinx(j)):(q+einx(j)),:),dfgrp(j))
10005 continue
      end


c***********************************************************************
c     LEVEL IV function: tuning parameter selection
c***********************************************************************

c     L1+GMCP
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l1cvbi(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call l1maxbi(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l1inibi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l1lmbi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call bi1msev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end
      
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl1gmcpbi(out,pauc,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call bfl1gmcpbi(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call l1cvbi(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L1+GSCAD
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l1cvscadbi(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvtas(p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     column-wise standardization
      call standard(cvtas,cvtsz,cvtz,cvtn,p)
c     solution path along lambdas, initial values for along kas
      call l1maxbi(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l1inibi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l1scadbi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call bi1msev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 i=1,p
         cvbs(i,:)=cvtas(i)*cvbs(i,:)
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end
      
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl1gscadbi(out,pauc,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call bfl1gscadbi(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call l1cvscadbi(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),
     +        lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L2+GMCP
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l2cvbi(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,
     +     ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvcholmat(p,p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     groupwise standardization
      do 1000 j=1,ngrp
         call grpstd(cvcholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        cvtsz(:,sinx(j):einx(j)),cvtz(:,sinx(j):einx(j)),cvtn,
     +        dfgrp(j))
 1000 continue
c     solution path along lambdas, initial values for along kas
      call l2maxbi(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l2inibi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l2lmbi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call bimsev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cvcholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        cvbs(sinx(j):einx(j),:),dfgrp(j))
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl2gmcpbi(out,pauc,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call bfl2gmcpbi(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call l2cvbi(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),lmdas,kas,
     +        nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end

c     L2+GSCAD
c***********************************************************************
c     CV comonent:  get alpha, beta for the CV dataset
c***********************************************************************
      subroutine l2cvscadbi(pauc,full,cvx,lmdas,kas,nlmda,nka,cvtn,cvpn,
     +     q,p,ngrp,dfgrp,sinx,einx,cvty,cvtx,cvtz,epsilon,maxit,cvpy,
     +     cvpx,cvpz)
      integer nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),
     +     einx(ngrp),maxit,
     +     full(nka*nlmda),cvx(nka*nlmda)
      double precision pauc(nlmda*nka),lmdas(nlmda),kas(nka),
     +     cvty(cvtn),cvtx(cvtn,q),cvtz(cvtn,p),epsilon,cvpy(cvpn),
     +     cvpx(cvpn,q),cvpz(cvpn,p)
c     local vars
      integer j,i
      double precision cvcholmat(p,p),cvtsz(cvtn,p),
     +     nulllmda,cvalpha(q),cvbeta(p),cvxinvxtx(cvtn,q),
     +     cvtr(cvtn),cvrmat(cvtn,nlmda),inicva(q,nlmda),
     +     inicvb(p,nlmda),cvas(q,nlmda*nka),cvbs(p,nlmda*nka),
     +     cvpeta(cvpn),cvpi(cvpn),tmp
!     groupwise standardization
      do 1000 j=1,ngrp
         call grpstd(cvcholmat(1:dfgrp(j),sinx(j):einx(j)),
     +        cvtsz(:,sinx(j):einx(j)),cvtz(:,sinx(j):einx(j)),cvtn,
     +        dfgrp(j))
 1000 continue
c     solution path along lambdas, initial values for along kas
      call l2maxbi(nulllmda,cvalpha,cvbeta,cvxinvxtx,cvtr,cvty,cvtx,
     +     cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
      do 1004 i=1,nlmda
         if (lmdas(i) .ge. nulllmda) then
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         else
            call l2inibi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),cvty,
     +           cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx,epsilon,maxit)
            inicva(:,i)=cvalpha(:)
            inicvb(:,i)=cvbeta(:)
            cvrmat(:,i)=cvtr(:)
         endif
1004  continue
c     solution path along kappa
      if (nka .eq. 1) then
         cvas(:,:)=inicva(:,:)
         cvbs(:,:)=inicvb(:,:)
      else
         do 1005 i=1,nlmda
            if (lmdas(i) .ge. nulllmda) then
               do 1006 j=1,nka
                  cvas(:,(i-1)*nka+j)=inicva(:,i)
                  cvbs(:,(i-1)*nka+j)=inicvb(:,i)
 1006          continue
            else
               cvalpha(:)=inicva(:,i)
               cvbeta(:) =inicvb(:,i)
               cvtr(:)=cvrmat(:,i)
c     ....     j=1 case 
               j=1
               cvas(:,(i-1)*nka+j)=cvalpha(:)
               cvbs(:,(i-1)*nka+j)=cvbeta(:)
               do 1007 j=2,nka
                  call l2scadbi(cvalpha,cvbeta,cvtr,cvxinvxtx,lmdas(i),
     +                 kas(j),cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,
     +                 sinx,einx,epsilon,maxit)
                  cvas(:,(i-1)*nka+j)=cvalpha(:)
                  cvbs(:,(i-1)*nka+j)=cvbeta(:)
 1007          continue
            endif
 1005    continue
      endif
!     determin whether satured or not 0: Not, 1: satured
      do 10006 i=1,nka*nlmda
         call bimsev2(full(i),cvx(i),cvas(:,i),cvbs(:,i),
     +        cvty,cvtx,cvtsz,cvtn,q,p,ngrp,dfgrp,sinx,einx)
10006 continue
!     change the penalized coefficients back
      do 10005 j=1,ngrp
         call dtrmm("L","U","N","N",dfgrp(j),(nka*nlmda),1.d0,
     +        cvcholmat(1:dfgrp(j),sinx(j):einx(j)),dfgrp(j),
     +        cvbs(sinx(j):einx(j),:),dfgrp(j))
10005 continue
c     get the paucs for the prediction set
      do 1011 i=1,(nka*nlmda)
         do 1012 j=1,cvpn
            cvpeta(j)=dot_product(cvpx(j,:),cvas(:,i))+dot_product(
     +           cvpz(j,:),cvbs(:,i))
 1012    continue
         do 1013 j=1,cvpn
            cvpi(j)=1.d0/(1.d0+exp(-cvpeta(j)))
1013     continue
         call aucroc(tmp,cvpy,cvpi,cvpn)
         pauc(i)=tmp
1011  continue
      end
	  
c***********************************************************************
c     Algorithm I  CV for lambda and kappa
c***********************************************************************
      subroutine cvl2gscadbi(out,pauc,olmdas,okas,ocoef,
     +     ofull,odf,ocvx,cvfull,cvcvx,
     +     nindex,cvk,y,x,z,n,q,p,ngrp,dfgrp,nka,kas,nlmda,minlmda,
     +     epsilon,maxit)
      integer cvk,n,q,p,ngrp,nka,nlmda,nindex(n),dfgrp(ngrp),maxit,
     +     ofull(nka*nlmda),odf(nka*nlmda),ocvx(nka*nlmda),
     +     cvfull(nka*nlmda),cvcvx(nka*nlmda)
      double precision out(3+q+p),pauc(nka*nlmda),y(n),x(n,q),z(n,p),
     +     kas(nka),minlmda,epsilon
c     local vars
      integer i,sinx(ngrp),einx(ngrp),cv,cvtn,cvpn,k,
     +     cvfullm(nka*nlmda,cvk),cvcvxm(nka*nlmda,cvk),
     +     use,loc
      double precision olmdas(nka*nlmda),okas(nka*nlmda),
     +     ocoef(q+p,nka*nlmda),lmdas(nlmda),kkas(nka),
     +     cvpy(n),cvpx(n,q),cvpz(n,p),cvty(n),cvtx(n,q),cvtz(n,p),
     +     cvpauc(nka*nlmda,cvk),
     +     omlmdas(nka*nlmda),omkas(nka*nlmda),omcoef(q+p,nka*nlmda),
     +     ompauc(nka*nlmda)
c     solution path of original dataset
      call bfl2gscadbi(olmdas,okas,ocoef,odf,ocvx,y,x,z,n,q,p,
     +     ngrp,dfgrp,nka,kas,nlmda,minlmda,epsilon,maxit)
c     Determine sature model, 0:not, 1: satured
      do 10000 i=1,nka*nlmda
         if ( odf(i) .le. n) then
            ofull(i)=0
         else
            ofull(i)=1
         endif
10000 continue
c     get lambdas,kappas
      do 10005 i=1,nlmda
         lmdas(i)= olmdas((i-1)*nka+1)
10005 continue
      kas(:)=okas(1:nka)
c     get sinx einx
      call fseinx(sinx,einx,dfgrp,ngrp)
c     start cross validation
      do 10006 cv=1,cvk
         cvtn=0
         cvpn=0
         do 10007 k=1,n
            if (nindex(k) .eq. cv) then
               cvpn=1+cvpn
               cvpy(cvpn)=y(k)
               cvpx(cvpn,:)=x(k,:)
               cvpz(cvpn,:)=z(k,:)
            else
               cvtn=1+cvtn
               cvty(cvtn)=y(k)
               cvtx(cvtn,:)=x(k,:)
               cvtz(cvtn,:)=z(k,:)
            end if
10007    continue
c     .... solve alpha, beta, and pauc for cv datasets
         call l2cvscadbi(cvpauc(:,cv),cvfullm(:,cv),cvcvxm(:,cv),
     +        lmdas,kas,nlmda,nka,cvtn,cvpn,q,p,ngrp,dfgrp,sinx,einx,
     +        cvty(1:cvtn),cvtx(1:cvtn,:),cvtz(1:cvtn,:),
     +        epsilon,maxit,
     +        cvpy(1:cvpn),cvpx(1:cvpn,:),cvpz(1:cvpn,:))
10006 continue
!     get average pauc, satured, convex info  in cross-validation
      do 10008 i=1,(nlmda*nka)
         pauc(i)=sum(cvpauc(i,:))/dble(cvk)
         if ( sum(cvfullm(i,:)) .eq. 0) then
            cvfull(i)=0
         else
            cvfull(i)=1
         endif
!     .....if ( sum(cvcvx(i,:)) .eq. 5) cvcvx(i)=1
10008 continue
!     assign useful output
      use=0
      do 10009 i=1,nka*nlmda
         if ( (ofull(i) .eq. 0) .and. (cvfull(i) .eq. 0) ) then
            use=use+1
            ompauc(use) =pauc(i)
            omlmdas(use)=olmdas(i)
            omkas(use)  =okas(i)
            omcoef(:,use)=ocoef(:,i)
         endif
10009 continue
!     get maximum
      loc=maxloc( ompauc(1:use), use )
      out(1)=ompauc(loc)
      out(2)=omlmdas(loc)
      out(3)=omkas(loc)
      out(4:(3+q+p))=omcoef(:,loc)
      end


c     L1 group
c     Given solution, kappa,lmda,
c     determine model size,aic, bic, objective value
C     Note: convex dianosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bi1msevab(ms,mszg,evgx,aic,bic,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,mszg,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision aic,bic,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer j,idx,i,info,lwork
      double precision alpha(q),beta(p),xnz(n,q+p),anb(q+p),eta(n),
     +     pi(n),logl
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups
      mszg=0
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
         endif
 1000 continue
c     .. non-zero individuals
      xnz(:,1:q)=x(:,:)
      anb(1:q)=alpha(:)      
      ms=q
      do 1001 idx=1,p
         if ( beta(idx) .ne. 0.d0 ) then
            ms=ms+1
            xnz(:,ms)=z(:,idx)
            anb(ms)=beta(idx)
         endif
 1001 continue
!     .. sum of square, aic, bic
      call dgemv("N",n,ms,1.d0,xnz(:,1:ms),n,anb(1:ms),1,0.d0,eta,1)
      logl=0.d0
      do 11001 i=1,n
         pi(i)=1/(1+exp(-eta(i)))
         logl=logl+ y(i)*eta(i) - log( 1.d0+exp(eta(i)) )
11001 continue
      aic= -2.d0*logl+ ms*2.d0
      bic= -2.d0*logl+ ms*log(dble(n))
      end
      
      

c     Given solution, kappa,lmda,
c     determine model size
C     convex diagnosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bi1msev(ms,evgx,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer idx
      double precision alpha(q),beta(p)
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero individuals
      ms=q
      do 1001 idx=1,p
         if ( beta(idx) .ne. 0.d0 ) then
            ms=ms+1
         endif
 1001 continue
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bi1msev2(full,evgx,alpha,beta,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer full,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision alpha(q),beta(p),y(n),x(n,q),z(n,p)
!     .. local arguments
      integer mszg,ms,j,idx
c     model size
c     .. non-zero individuals
      ms=q
      do 1001 idx=1,p
         if ( beta(idx) .ne. 0.d0 ) then
            ms=ms+1
         endif
 1001 continue
      full=1
      if (ms .le. n) full=0
c     call intpr("ms",-1,ms,1)
c     call intpr("full",-1,full,1)
      end



c     L2 group
c     Given solution, kappa,lmda,
c     determine model size,aic, bic, objective value
C     Note: convex dianosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bimsevab(ms,mszg,evgx,aic,bic,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,mszg,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision aic,bic,ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer j,idx,i,info,lwork
      double precision alpha(q),beta(p),xnz(n,q+p),anb(q+p),eta(n),
     +     pi(n),logl
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups, individuals
      xnz(:,1:q)=x(:,:)
      anb(1:q)=alpha(:)      
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
               xnz(:,ms)=z(:,idx)
               anb(ms)=beta(idx)
 1001       continue
         endif
 1000 continue
!     .. sum of square, aic, bic
      call dgemv("N",n,ms,1.d0,xnz(:,1:ms),n,anb(1:ms),1,0.d0,eta,1)
      logl=0.d0
      do 11001 i=1,n
         pi(i)=1/(1+exp(-eta(i)))
         logl=logl+ y(i)*eta(i) - log( 1.d0+exp(eta(i)) )
11001 continue
      aic= -2.d0*logl+ ms*2.d0
      bic= -2.d0*logl+ ms*log(dble(n))
      end


c     Given solution, kappa,lmda,
c     determine model size
C     convex diagnosis is not performed, though space is left
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bimsev(ms,evgx,ab,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx,ka,lmda)
      integer ms,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision ab(q+p),y(n),x(n,q),z(n,p),ka,lmda
!     .. local arguments
      integer mszg,j,idx
      double precision alpha(q),beta(p)
      alpha(:)=ab(1:q)
      beta(:)=ab((q+1):(q+p))
c     model size
c     .. non-zero groups, individuals
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
 1001       continue
         endif
 1000 continue
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bimsev2(full,evgx,alpha,beta,y,x,z,n,q,p,
     +     ngrp,dfgrp,sinx,einx)
      integer full,evgx,n,q,p,ngrp,dfgrp(ngrp),sinx(ngrp),einx(ngrp)
      double precision alpha(q),beta(p),y(n),x(n,q),z(n,p)
!     .. local arguments
      integer mszg,ms,j,idx
c     model size
c     .. non-zero groups, individuals
      mszg=0
      ms=q
      do 1000 j=1,ngrp
         if ( sum(abs( beta(sinx(j):einx(j)) )) .ne. 0.d0) then
            mszg=mszg+1
            do 1001 idx=sinx(j),einx(j)
               ms=ms+1
 1001       continue
         endif
 1000 continue
      full=1
      if (ms .le. n) full=0
c     call intpr("ms",-1,ms,1)
c     call intpr("full",-1,full,1)
      end
