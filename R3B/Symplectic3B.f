c An example of data file: (see comments in the code)
c   1  30  1  1000.           Iwr, NS Ist  tmax (iwr & ist affect writing, NS accuracy , tmax time interval) 
c   1.0   0 0 0   0 0 0          'SUN'    mass_1  coords_1  vels_1
c   1.d-6  1 0 .1   0  1 0        'EARTH' mass_2  .............
c   1.e-3  5  .1 .2  0 .27 0      'JUPITER' (very eccentric)
C try and see what the code writes 
c-------------------------------------------------------------------------
c       CODE: G3.f : reads data from 'standard input'
        IMPLICIT REAL*8 (A-H,M,O-Z)
        PARAMETER(NMX=3,NMX3=3*NMX)
        parameter(PI=3.141592653589793D0,
     &  GaussK=0.01720209895D0,twopi=pi+pi)
        real*8 xw(9),vw(9)
        real*8 M(NMX),X(NMX3),V(NMX3)
        real*8 XJ(NMX3),VJ(NMX3)
        real*8 cmx(9),cmv(9)
        real*8 inc12,inc3 
        Common /AAA/A1,A2,Ene0,ax1,prd1
        common/Diagno/Energy,iwr
        common/paskel/hs
c       Read data from file 'input'
          open(55,file='input')
c        IWR=writing index, NS = (average) # steps/binary period,
c        Ist = Output after every Ist steps, tmax= maximum time
1        READ(55,*,end=999)Iwr,NS,Ist,yearmax,a1,a2
          TMAX=twopi*yearmax
c          these a's [ a1, a2 ] are parameters in the time transformation. Optimal values not known!
         !a1=1  ! this may be odopted 'correct'
         !a2=15 ! but nobody knows what this really should be!
         time=0.0
c         Read mass, x y z  vx vy vz in cm-coordinates 
c         (or heliocentric or in any other system for which the 
c         reduction to Jacobian coordinates (below) works).
         do i=1,3
          k1=3*i-2
          k2=k1+2
        read(55,*,end=999)m(i),(cmx(k),k=k1,k2),(cmv(k),k=k1,k2)
          end do
c       transform to Jacobian coordinates
         do i=1,3
         x(i)=cmx(3+i)-cmx(i)
         x(i+3)=cmx(6+i)-(cmx(3+i)*m(2)+cmx(i)*m(1))/(m(1)+m(2))
         v(i)=cmv(3+i)-cmv(i)
         v(i+3)=cmv(6+i)-(cmv(3+i)*m(2)+cmv(i)*m(1))/(m(1)+m(2))
         end do

c         open output file
         open(66,file='output')
            do k=1,6
            xj(k)=x(k)
            vj(k)=v(k)
            end do

           do  k=1,6
           xw(k)=xj(k)
           vw(k)=vj(k)
           end do
        call EnergyTest(xw,vw,m,Ene0,g,A1,A2,NS,Hs)
        ene00=ene0
        prd0=prd1
c         find the parametric 'symplectic space coordinates'
c         by application of Wisdom-Holman-Touma corrector
c         backwards. (iteration )
         do ite=1,10
         do k=1,6
         xw(k)=xj(k)
         vw(k)=vj(k)
         end do
         call cor(hs,xw,vw,m)
         do k=1,6
         xj(k)=xj(k)+x(k)-xw(k)
         vj(k)=vj(k)+v(k)-vw(k)
         end do
         end do
         do k=1,6
         xw(k)=xj(k)
         vw(k)=vj(k)
         end do
        call cor(hs,xw,vw,m)
        call EnergyTest(xw,vw,m,Ene1,g,dum,dum,NS,dum)
        call EnergyTest(xj,vj,m,Ene0,g,dum,dum,ns,dum)
        prd0=prd1
        ene0=ene1
c     Here begins the integration
          call SI3(Hs/2,Xj,Vj,M,Time,0)
c        Integration stops if t>tmax or more that a billion steps 
         do is=0,1000 000 000
            call SI3(Hs,Xj,Vj,M,Time,1)
            ytime=time/twopi
            if(is.eq.is/Ist*Ist)then
c           write the JACOBIAN coordinates. 
         do k=1,6
         xw(k)=xj(k)
         vw(k)=vj(k)
         end do
          TY=Time  !added
            call SI3(Hs/2,Xw,Vw,M,Ty,1) !added
            ytime=Ty/twopi ! added
        call cor(hs,xw,vw,m)
        call EnergyTest(xw,vw,m,Enet,g,dum,dum,NS,dum)
         dene=enet-ene00
            write(6,106)ytime,dene*g,enet/ene00-1!,enet
106         format(1x,f12.3,1p,2g9.1)!,g16.8)
            if(iwr.gt.0)then ! WRITE IF IWR>0 !!!!!!!!!!!!!!!
            write(66,123)ytime,(xw(j),j=1,6)!,(vw(j),j=1,6) ! write coords
            end if 
123        format(1x,f12.2,6f12.6,' 0 0 0')
                m12=m(1)+m(2)
                m3=m12+m(3)
        call elmnts
     & (xw(1),vw(1),m12,a12,e12,mo12,inc12,Om12,oo12,alfa12,q12,tq12)
        call elmnts
     & (xw(4),vw(4),m3,a3,e3,mo3,inc3,Om3,oo3,alfa3,q3,tq3)
       write(91,191)ytime,a12,e12,inc12,a3,e3,inc3 ! write some elemets
c                   t/prd1  a ecc  inc  a  ecc  inc
191    format(1x,f12.3,2f9.5,f7.2,2f9.5,f7.2)
          if(time.gt.tmax)goto 1
          end if
         end  do! is
99           continue 
        GOTO 1
999        END
        SUBROUTINE SI3(h,x,v,m,time,Ifi)
            IMPLICIT REAL*8 (A-H,M,O-Z)
        real*8 x(9),v(9),m(3),ga(0:5),gb(0:5),GN(5),dv1(3),dv2(3)
        real*8 xa(3),xb(3),dg1(3),dg2(3),dr1(3),dr2(3),dp1(3),dp2(3)
          real*8 oy2(5)!,ot(5)
        common/AAA/A1,A2,ener0,ax1,prd1
        common/diagno/Energy,iwr
        common/paskel/hs
        common/mass/mp1,mp2,my1,my2,af0,af1,af2,ak1,ak2
        SAVE
        data td,y1,y2,rest/4*0.0/,dv1,dv2/6*0.0/
c       Auxiliary quantities
        twopi=8.d0*atan(1.d0)
         AR=a2/a1
        mj1=m(1)+m(2)
        mj2=mj1+m(3)
c        mp1=m(1)*m(2)
c        mp2=mj1*m(3)
c        my1=mp1/mj1
c        my2=mp2/mj2
c        af0=(m(1)+m(2))*m(3)
c        af1=-m(1)*m(3)
c        af2=-m(2)*m(3)
c        ak1=m(2)/(m(1)+m(2))
c        ak2=ak1-1.
        r10=sqrt(cdot(x,x))
        r20=sqrt(cdot(x(4),x(4)))
c        Perturbations
c-----------------------------------------------------------V
         do k=1,3
         xa(k)=x(3+k)+ak1*x(k)
         xb(k)=x(3+k)+ak2*x(k)
         end do
         ra=sqrt(cdot(xa,xa))
         rb=sqrt(cdot(xb,xb))
         Rest=af0/r20+af1/ra+af2/rb
         g=1./(A1/r10+A2/r20)
c************************************************************
         yr13=1/r10**3
         yr23=1/r20**3
         yra3=1/ra**3
         yrb3=1/rb**3
         if(ifi.gt.0)then
          do k=1,3
          dg1(k)=yr13*g**2*A1*x(k)
          dg2(k)=yr23*g**2*A2*x(k+3)
          dr1(k)=-(af1*ak1*yra3*xa(k)+af2*ak2*yrb3*xb(k))
          dr2(k)=-(af1*yra3*xa(k)+af2*yrb3*xb(k)+af0*x(3+k)*yr23)
          dp1(k)=-dg1(k)*Rest-g*dr1(k)
          dp2(k)=-dg2(k)*Rest-g*dr2(k)
          dv1(k)=hs*dp1(k)/my1
          dv2(k)=hs*dp2(k)/my2
          v(k  )=v(k  )+dv1(k)
          v(3+k)=v(3+k)+dv2(k)
          end do
c-----------------------------------------------------------A
          end if

c       Two-body motions
        W12=CDOT(v(1),v(1))
        w22=cdot(v(4),v(4))
        beta1=2.D0*mj1/R10-W12
          if(beta1.lt.0.0)NEGAX=1
        beta2=2.d0*mj2/r20-w22
        Ener2=(-.5d0*my1*beta1-.5d0*my2*beta2)
c          ENERGY EVALUATION
        Energy=Ener2+Rest
c*************************************
        eps=g*(Ener2-Ener0)
        eps1=eps*A1/mp1
        eps2=eps*A2/mp2
c       Modified masses
        mj1=mj1*(1+eps1)
        mj2=mj2*(1+eps2)
c       Modified two-body motions
        beta1=2*mj1/r10-w12
        beta2=2*mj2/r20-w22
        ETA1=CDOT(X(1),v(1))
        eta2=cdot(x(4),v(4))
        ZETA1=mj1-beta1*R10
        zeta2=mj2-beta2*r20

         if(h.eq.hs)then
         kie=kie+1 
         do i=5,2,-1
         oy2(i)=oy2(i-1)
         end do
         oy2(1)=y2van
        if(kie.gt.5)then
        y2=4*oy2(1)-6*oy2(2)+4*oy2(3)-oy2(4)
         y2van=y2
        end if
         end if

          vi1=0.0
          vi2=h/A2
          p1=-1
          p2= 1
           if(h.ne.hs)y2=.5*y2
           if(y2.lt.0.0.or.y2.gt.vi2)y2=.5*(vi2)
                                      
         do 10 ite=1,100
         y1=(h-A2*y2)/A1
         z1=beta1*y1*y1
         z2=beta2*y2*y2
         call gfun(beta1,y1,ga)
         call gfun(beta2,y2,gb)
          R1=R10+eta1*ga(1)+zeta1*ga(2)
          R2=R20+eta2*gb(1)+zeta2*gb(2)
          t1=r10*y1+eta1*ga(2)+zeta1*ga(3)
          t2=r20*y2+eta2*gb(2)+zeta2*gb(3)
          ft=t2-t1
          f1=R2+AR*R1
          if(ft.lt.0.0)then
          ak=vi2
          Pk=p2
          vi1=y2
          p1=ft
          else
          ak=vi1
          Pk=p1
          vi2=y2
          p2=ft
          end if

          et1=eta1*ga(0)+zeta1*ga(1)
          et2=eta2*gb(0)+zeta2*gb(1)
          ze1=zeta1*ga(0)-beta1*eta1*ga(1)
          ze2=zeta2*gb(0)-beta2*eta2*gb(1)
          f2=et1-(ar)**2*et2
          f3=ze1+(ar)**3*ze2
          dy2=-ft/f1
          dy2=-ft/(f1+dy2*f2/2)
          dy2=-ft/(f1+dy2*(f2/2+f3*dy2/6))
          dy1=-AR*dy2
          test=a1*abs(dy1)+a2*abs(dy2)
         if(test.lt.1.e-9*h)goto 11
          y2=y2+dy2
          if( (y2-vi1).lt.0.0)y2=.5d0*(vi1+vi2)
          if( (vi2-y2).lt.0.0)y2=.5d0*(vi1+vi2)
          if(ite.eq.ite/7*7)y2=.5d0*(vi1+vi2)
          y1=(h-A2*y2)/A1
1231    format(1x,2i10,1p,9e12.4)
10       continue
             write(6,*)' no convergence'
             STOP
11         continue
         if(h.eq.hs)then
         y2van=y2
         end if

         ga(5)=ga(5)+ga(4)*dy1+ga(3)*dy1**2/2
         ga(4)=ga(4)+ga(3)*dy1+ga(2)*dy1**2/2
         gb(5)=gb(5)+gb(4)*dy2+gb(3)*dy2**2/2
         gb(4)=gb(4)+gb(3)*dy2+gb(2)*dy2**2/2
         y1=y1+dy1
         y2=y2+dy2
         ga(3)=y1**3/6-beta1*ga(5)
         gb(3)=y2**3/6-beta2*gb(5)
         ga(2)=y1**2/2-beta1*ga(4)
         gb(2)=y2**2/2-beta2*gb(4)
         ga(1)=y1-beta1*ga(3)
         gb(1)=y2-beta2*gb(3)
         r1=r10+eta1*ga(1)+zeta1*ga(2)
         r2=r20+eta2*gb(1)+zeta2*gb(2)
         td=r20*y2+eta2*gb(2)+zeta2*gb(3)
         time=time+td
         do 22 k=1,2
         if(k.eq.1)then
         do i=1,5
         gn(i)=ga(i)
         end do
         r0=r10
         rx=r1
         m12=mj1
         else
         do i=1,5
         gn(i)=gb(i)
         end do
         r0=r20
         rx=r2
         m12=mj2
         end if

        F=1.D0-m12*GN(2)/R0
        G=(td-m12*GN(3))
        dF=-m12*GN(1)/(R0*RX)
        dG=1.D0-m12*GN(2)/RX
        i1=1+(k-1)*3
        i3=i1+2
        DO 1 I=i1,i3
        W0I=v(I)
        v(I)=dF*X(I)+dG*v(I)
        X(I)= F*X(I)+ G*W0I
1       CONTINUE
22      continue
        RETURN
        END
       FUNCTION CDOT(A,B)
       IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 A(3),B(3),CDOT
       CDOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
       RETURN
       END

C
       subroutine gfun(beta,y,g)
       implicit real*8 (a-h,o-z)
        real*8 g(0:5),c(5)
        z=beta*y*y
        call cfun(z,c)
        g(0)=1-z*c(2)
        s=y
       do i=1,5
       g(i)=s*c(i)
       s=s*y
       end do
       return
       end
       SUBROUTINE CFUN(Z,C)
       IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(C6=1.D0/6.D0,C132=1.D0/132.D0,C56=1.D0/56.D0,
     &  C30=1.D0/30D0,C24=1.D0/24.D0,C156=1.D0/156.D0,
     &  C90=1.D0/90.D0,C110=1.D0/110.D0,C16=1.D0/16.D0,C8=1.D0/8.D0,
     &  C72=1.D0/72.D0,C42=1.D0/42.D0,C120=1.D0/120.D0,U=1.d0)
       DIMENSION C(5)
       data large/0/
       h=z
       DO 1 k=0,10
       IF(ABS(H).lt. 0.1)goto 2
       h=0.25d0*h
 1     CONTINUE
         if(large.ge.10)stop
         if(large.lt.10)then
         large=large+1
         write(6,*)' too large Z in CFUN!',z,char(7)
         end if
2       CONTINUE
        C(4)=(U-H*(U-H*(U-H*C90/(U+H*C132))*C56)*C30)*C24
        C(5)=(U-H*(U-H*(U-H*C110/(U+H*C156))*C72)*C42)*C120
       DO 4 I=1,K
       C3=c6-h*C(5)
       C2=.5D0-h*C(4)
       C(5)=(C(5)+C(4)+C2*C3)*c16
       C(4)=C3*(2.d0-h*C3)*c8
 4     h=4.d0*h
       C(3)=c6-Z*C(5)
        C(2)=.5D0-Z*C(4)
        C(1)=1.d0-Z*C(3)
       RETURN
       END
        SUBROUTINE EnergyTest(x,v,m,Energy,g,an1,an2,NS,Hs)
        IMPLICIT REAL*8 (A-H,M,O-Z)
        real*8 x(9),v(9),m(3)
        real*8 xa(3),xb(3)
        common/AAA/A1,A2,dummy,ax1,prd1
        common/mass/mp1,mp2,my1,my2,af0,af1,af2,ak1,ak2
        data i0/0/
        save
        twopi=8.d0*atan(1.d0)
c       Auxiliary quantities
        mj1=m(1)+m(2)
        mj2=mj1+m(3)
        mp1=m(1)*m(2)
        mp2=mj1*m(3)
        my1=mp1/mj1
        my2=mp2/mj2
        af0=(m(1)+m(2))*m(3)
        af1=-m(1)*m(3)
        af2=-m(2)*m(3)
        ak1=m(2)/(m(1)+m(2))
        ak2=ak1-1.
        r10=sqrt(cdot(x,x))
        r20=sqrt(cdot(x(4),x(4)))
c        Perturbations
c-----------------------------------------------------------V
         do k=1,3
         xa(k)=x(3+k)+ak1*x(k)
         xb(k)=x(3+k)+ak2*x(k)
         end do
         ra=sqrt(cdot(xa,xa))
         rb=sqrt(cdot(xb,xb))
         Rest=af0/r20+af1/ra+af2/rb
        W12=CDOT(v(1),v(1))
        w22=cdot(v(4),v(4))
        beta1=2.D0*mj1/R10-W12
        beta2=2.d0*mj2/r20-w22
        alf1=beta1/mj1
        alf2=beta2/mj2
        ax1=1/alf1
        Prd1=twopi*ax1*sqrt(ax1/mj1)
        hs=Prd1/NS
        Ener2=(-.5d0*my1*beta1-.5d0*my2*beta2)
c          ENERGY EVALUATION
        Energy=Ener2+Rest
        if(i0.eq.0)then
        ga=1/(A1*alf1+A2*alf2)
        an1=A1*ga
        an2=A2*ga
        i0=1
        end if
        g=1/(a1/r10+a2/r20)
       RETURN
       END
        SUBROUTINE Cor(h,xj,vj,m)
        IMPLICIT REAL*8 (A-H,M,O-Z)
        real*8 x(6),m(3),dv1(3,3),dv2(3,3),xj(6),vj(6),gam(3)
        real*8 xa(3),xb(3),dg1(3),dg2(3),dr1(3),dr2(3),dp1(3),dp2(3)
        common/AAA/A1,A2,ener0,ax1,prd1
        common/diagno/Energy,iwr
        common/mass/mp1,mp2,my1,my2,af0,af1,af2,ak1,ak2
        save
c       Auxiliary quantities
        mj1=m(1)+m(2)
        mj2=mj1+m(3)
c        mp1=m(1)*m(2)
c        mp2=mj1*m(3)
c        my1=mp1/mj1
c        my2=mp2/mj2
c        af0=(m(1)+m(2))*m(3)
c        af1=-m(1)*m(3)
c        af2=-m(2)*m(3)
c        ak1=m(2)/(m(1)+m(2))
c       ak2=ak1-1.

         r10=sqrt(cdot(xj(1),xj(1)))
         r20=sqrt(cdot(xj(4),xj(4)))
         g0=1/(A1/r10+A2/r20)
         ceps=.1*g0
         s=ceps*h
         twos=s+s
        Do kh=1,3
        s=-s
        if(kh.eq.3)s=0.
        do j=1,6
        x(j)=xj(j)+s*vj(j)
        end do

        r10=sqrt(cdot(x,x))
        r20=sqrt(cdot(x(4),x(4)))
         do k=1,3
         xa(k)=x(3+k)+ak1*x(k)
         xb(k)=x(3+k)+ak2*x(k)
         end do
         ra=sqrt(cdot(xa,xa))
         rb=sqrt(cdot(xb,xb))
         Rest=af0/r20+af1/ra+af2/rb
         g=1./(A1/r10+A2/r20)
         gam(kh)=g*Rest
         yr13=1/r10**3
         yr23=1/r20**3
         yra3=1/ra**3
         yrb3=1/rb**3
          do k=1,3
          dg1(k)=yr13*g**2*A1*x(k)
          dg2(k)=yr23*g**2*A2*x(k+3)
          dr1(k)=-(af1*ak1*yra3*xa(k)+af2*ak2*yrb3*xb(k))
          dr2(k)=-(af1*yra3*xa(k)+af2*yrb3*xb(k)+af0*x(3+k)*yr23)
          dp1(k)=-dg1(k)*Rest-g*dr1(k)
          dp2(k)=-dg2(k)*Rest-g*dr2(k)
          dv1(k,kh)=dp1(k)/my1
          dv2(k,kh)=dp2(k)/my2
          end do
           end do
                      hc=-h**2/24
             rdot=(gam(2)-gam(1))/(twos)

           do k=1,3
           xj(k  )=xj(k  )+g0*hc*dv1(k,3)
        vj(k  )=vj(k  )-g0*hc*(dv1(k,2)-dv1(k,1))/(twos)
     &     + dg1(k)*rdot*hc /my1
           xj(k+3)=xj(k+3)+g0*hc*dv2(k,3)
        vj(k+3)=vj(k+3)-g0*hc*(dv2(k,2)-dv2(k,1))/(twos)
     &     + dg2(k)*rdot*hc /my2
           end do
        RETURN
        END
       SUBROUTINE JtoCM(XJ,VJ,M,X,V)
       IMPLICIT REAL*8 (A-H,M,O-Z)
       REAL*8 X(9),M(3),XJ(6),cm(3),cv(3),v(9),vj(6)

C
       x(1)=0
       x(2)=0
       x(3)=0
       x(4)=xj(1)
       x(5)=xj(2)
       x(6)=xj(3)
       p=m(2)/(m(1)+m(2))
       x(7)=xj(1)*p+xj(4)
       x(8)=xj(2)*p+xj(5)
       x(9)=xj(3)*p+xj(6)
       v(1)=0
       v(2)=0
       v(3)=0
       v(4)=vj(1)
       v(5)=vj(2)
       v(6)=vj(3)
       v(7)=vj(1)*p+vj(4)
       v(8)=vj(2)*p+vj(5)
       v(9)=vj(3)*p+vj(6)
               do ite=1,2
        do k=1,3
        cm(k)=x(k)*m(1)
        cv(k)=v(k)*m(1)
        end do

        sm=m(1)
        do i=2,3
        sm=sm+m(i)
        do k=1,3
        cm(k)=cm(k)+m(i)*x(3*i-3+k)
        cv(k)=cv(k)+m(i)*v(3*i-3+k)
        end do ! k
        end do ! i
c        cm ok

        do i=1,3
        do k=1,3
        x(3*i-3+k)=x(3*i-3+k)-cm(k)/sm
        v(3*i-3+k)=v(3*i-3+k)-cv(k)/sm
        end do! k
        end do! i
         end do !ite
        return
        end

       SUBROUTINE JtoCMoB(XJ,VJ,M,X,V)
       IMPLICIT REAL*8 (A-H,M,O-Z)
       REAL*8 X(9),M(3),XJ(6),v(9),vj(6)
C
       p=m(2)/(m(1)+m(2))
       x(7)=+xj(4)
       x(8)=+xj(5)
       x(9)=+xj(6)
       x(1)=-p*xj(1)
       x(2)=-p*xj(2)
       x(3)=-p*xj(3)
       x(4)=(1-p)*xj(1)
       x(5)=(1-p)*xj(2)
       x(6)=(1-p)*xj(3)
        return
        end

        subroutine elmnts
     & (x,v,m,a,e,mo,inc,Om,oo,alfa,q,tq)
c       NOTE: wrong results can be produced in exeptional situations
c       where some angles are undefined in terms of the expressions used.
c       This may happen in exactly planar, rectilinear .. orbits
c       Troubles can often be avoided by a very small 'perturbation' of x and/or v.
        implicit real*8 (a-h,m,o-z)
        parameter(rad=180.d0/3.141592653589793d0 )
        real*8 x(3),w(3),v(3),inc,jx,jy,jz
        mu=sqrt(m)
        do k=1,3
        w(k)=v(k)/mu
        end do
        r=sqrt(x(1)**2+x(2)**2+x(3)**2)
        w2=w(1)**2+w(2)**2+w(3)**2
        eta=x(1)*w(1)+x(2)*w(2)+x(3)*w(3)
        alfa=2/r-w2
        zeta=1-alfa*r

c       areal velocity vector (jx,jy,jz)
        jx=x(2)*w(3)-x(3)*w(2)
        jy=x(3)*w(1)-x(1)*w(3)
        jz=x(1)*w(2)-x(2)*w(1)
        d=sqrt(jx*jx+jy*jy+jz*jz)

c       eccentricity vector (ex,ey,ez)
        ex=w(2)*jz-w(3)*jy-x(1)/r
        ey=w(3)*jx-w(1)*jz-x(2)/r
        ez=w(1)*jy-w(2)*jx-x(3)/r

        e=sqrt(ex*ex+ey*ey+ez*ez)
        b=sqrt(jx*jx+jy*jy)
        inc=atn2(b,jz)*rad
        Om=atn2(jx,-jy)*rad
        oo=atn2(ez*d,ey*jx-ex*jy)*rad
        a=1/alfa
        sqaf=sqrt(abs(alfa))
        q=d*d/(1+e) 
        too=oot(alfa,eta,zeta,q,e,sqaf)
        tq=too/mu
        mo=too*sqaf**3*rad
        return
        end
        function atn2(s,c)
        implicit real*8 (a-h,o-z)
        parameter(twopi=2*3.141592653589793d0)
        atn2=atan2(s,c)
        if(atn2.lt.0.0)atn2=atn2+twopi
        return
        end
        function oot(alfa,eta,zeta,q,e,sqaf) ! oot=pericentre time
c       alfa=1/a; eta=sqrt(a) e sin(E); zeta=e Cos(E),
c       q=a(1-e), e=ecc, sqaf=sqrt(|a|)
        implicit real*8 (a-h,o-z)
        parameter(tiny=1.d-18)
        save
        if(zeta.gt.0.0)then 
c        ellipse (near peri), parabola or hyperbola.
         ecc=max(e,tiny)
         X=eta/ecc
         Z=alfa*X*X
         oot=X*(q+X*X*g3(Z))
        else
c       upper half of an elliptic orbit.
        oot=(atan2(eta*sqaf,zeta)/sqaf-eta)/alfa
        end if
         return
         end
        function g3(z)
        implicit real*8 (a-h,o-z)
        common/mita/zero
        save
        if(z.gt.0.025d0)then ! elliptic
        x=sqrt(z)
        g3 = (asin(x)-x)/x**3
        elseif(z.lt.-0.025d0)then ! hyperbolic
        x = sqrt(-z)
        g3 = (log(x+sqrt(1+x*x))-x )/x/z
        else ! Pade approximant for small  |z|
c       g3 = (1/6.d0-19177*z/170280 + 939109*z*z/214552800)/
c     &  (1-7987*z/7095 + 54145*z*z/204336)
       g3 = (1+6*(-19177*z/170280 + 939109*z*z/214552800))/
     &  (6*(1-7987*z/7095 + 54145*z*z/204336))
        zero=0
        end if
        return
        end
