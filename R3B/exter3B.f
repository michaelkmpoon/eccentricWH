
!        A Three-Body Code, with the possibility of external
!        accelerations. To program your own version modify
!        the routines 'external acceleration' 
!        and  'external coordinate  accelerations' 
         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9)
         real*8 g0(3),gt(3)
         real*8 xcm(3),vcm(3) 
         common/diagno/neval
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         eval=0.0
         told=0
         open(66,file='xyz') ! output file for t x1 y1 z1 x2 y2 z2 x3 y3 z3
         read(5,*)tmx,tincr,eps,icm ! read some parameters 
c                   tmx=max time
c                      tincr=output interval (at least)
c                            eps= error tolerance
c                                icm= cm-index (>0=>reduce2cm)
         do i=1,3
         k1=i*3-2
         k2=k1+2
         read(5,*)m(i),(xr(k),k=k1,k2),(vr(k),k=k1,k2) ! read m, x & v (inertial coords)
         end do
         nb=3
         newcase=1
         time=0
         if(icm.gt.0)then ! cm-only?
         call reduce 2 cm(xr,m,xcm) 
         call reduce 2 cm(vr,m,vcm) 
         end if
         call Constants Of Motion(xr,vr,m,nb,ene0,g0,alag)
1        continue         
         write(66,123)time,XR!write (inertial system) coords
         call ThreeB(m,time,xr,vr,b,tincr,newcase,eps) 
c        begin{diagnostics}
         call Constants Of Motion(xr,vr,m,nb,enet,gt,alag)
         vrh=abs(enet+b)/alag
         glag=mass**2.5/abs(ene0)*sqrt(alag)
         dgg=(gt(1)-g0(1))**2+(gt(2)-g0(2))**2+(gt(3)-g0(3))**2
         dgg=sqrt(dgg)/glag
         write(6,106)time,vrh,dgg,neval,enet,time-told 
!         vrh=error in the internal 3B energy (actually = difference
!         of the integrated and computed values).          
!         dgg=three-body angular momentum, not necessarily conserved if you
!         have external forces
!         neval = number of acceleration evaluations
!         enet= 3B energy at the moment
!         time-told = true time difference between outputs 
         told=time
106      format(1x,' T: ',f12.4,1p,2g11.2,i9,1p,2g20.12)
123      format(1x,f12.4,1p,9g14.6)
         call flush(6)
         call flush(66)
         if(time.lt.tmx)goto 1
         end

        subroutine Reduce2cm(x,m,cm)
        real*8 x(*),m(*),cm(*)
        save
        cm(1)=0
        cm(2)=0
        cm(3)=0
        sm=0
        do i=1,3
        sm=sm+m(i)
        do k=1,3
        cm(k)=cm(k)+m(i)*x(k+3*(i-1))
        end do
        end do
        do k=1,3
        cm(k)=cm(k)/sm
        end do
        do i=1,3
        do k=1,3
        x(k+3*(i-1))=x(k+3*(i-1))-cm(k)
        end do
        end do
        
        return
        end

         subroutine ThreeB(mr,time,xr,vr,b,tincr,newcase,eps)
         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9),x(0:12),v(0:12),y(26),SY(26),mr(3)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         equivalence(y(1),x(0)),(y(14),v(0))
         save
c        initializations         
         if(newcase.eq.1)then
         newcase=0
         m(1)=mr(1)
         m(2)=mr(2)
         m(3)=mr(3)
         mass=m(1)+m(2)+m(3)
         mu(1)=m(1)/mass
         mu(2)=m(2)/mass
         mu(3)=m(3)/mass
         mm(1)=m(2)*m(3)
         mm(2)=m(1)*m(3)
         mm(3)=m(1)*m(2)
         mmm=mm(1)+mm(2)+mm(3)
         call reduce 2 cm(xr,m,x(10)) ! obtain cm-x [x(10),x(11),x(12)]
         call reduce 2 cm(vr,m,v(10)) ! obtain cm-v [v(10),v(11),v(12)]
         call Compute Triangle Differences(xr(1),x(1))
         call Compute Triangle Differences(vr(1),v(1))
         call Evaluate CM Energy(v,x,E)
         eta=0.01
         if(E.lt.0.0)then
         h=eta*mmm*sqrt(mass/abs(E))!this scales correctly
         else 
         h=eta*mmm*sqrt(mass/(tkin+E))
         end if
         Energy=E ! to common
         v(0)=-E ! binding energy
         x(0)=0  ! time
         end if ! newcase
c         
         told=x(0)
         nsteps=0
         NVX=20+6
1        continue
c        order of magnitude of variables
         rsmall=mmm/UPOT
         vsmall=sqrt(TKIN/mass)
         tsmall=rsmall*sqrt(rsmall/mass)
         esmall=abs(v(0))  
         sy(1)=tsmall
         sy(14)=esmall
         do k=2,10+3
         sy(k)=rsmall
         sy(k+13)=vsmall
         end do
            if(tincr.gt.0.0)h=min(h,1.1*tincr*Upot)
            JMAX=10
         call  DIFSYAB(NVX,EPS,SY,h,dummy,Y,JMAX) ! integrate 
         call Triangle Correction(v,x)
                 if(h.eq.0.0)then
                 write(6,*)' STEP=0!',char(7)
                 STOP
                 end if
         nsteps=nsteps+1
         if(x(0).lt.told+tincr.and.nsteps.lt.10000000)goto 1
         call evaluate cm energy(v,x,enow)
         call Transform 2 inertial(x,v,x(10),v(10),xr,vr,time)
c         write(67,*)x(0),x(10),x(11),x(12)
         b=v(0)
         return
         end

         subroutine Transform 2 inertial(x,v,xcm,vcm,xr,vr,time)
         implicit real*8 (a-h,o-z)
         real*8 x(0:9),v(0:9),xcm(3),vcm(3),xr(9),vr(9)
         time=x(0)
         call Transform to CM system(X(1),XR(1))
         call Transform to CM system(V(1),VR(1))
         do i=1,3
         i0=3*i-3
         do k=1,3 ! add cm-contributions
         xr(i0+k)=xr(i0+k)+xcm(k)!+x(0)*vcm(k) 
         vr(i0+k)=vr(i0+k)+vcm(k)
         end do !k
         end do !i 
         return
         end
         subroutine Triangle Correction(v,x)
         implicit real*8 (a-h,o-z)
         real*8 v(0:9),x(0:9)
         save
c        find longest distance
         rrmx=0
         do i=1,3
         i1=3*i-2
         rr=sqa(x(i1))
         if(rrmx.lt.rr)then
         rrmx=rr
         imx=i
         end if
         end do ! k

         ix0=3*imx-3
c        take -sum of two shorter ones for the longest
         do k=1,3
         x(ix0+k)=0
         v(ix0+k)=0
         end do

         do i=1,3
         i0=3*i-3
         if(i.ne.imx)then
         do k=1,3
         x(ix0+k)=x(ix0+k)-x(i0+k) ! re-evaluate the longest vector
         v(ix0+k)=v(ix0+k)-v(i0+k) ! and the velocity in the same way
         end do
         end if
         end do
         return
         end 
         
         subroutine SubSteps(Y,H,Leaps) ! Leaps substeps of size H/Leaps.
         implicit real*8 (a-h,o-z)
         real*8 y(*)
         save
         call LeapFrog(y(14),y(1),H,Leaps)
         return
         end
         subroutine LeapFrog(v,x,h,Leaps)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:12),v(0:12)!,xx(0:9),vv(0:9)
         save
         call Move Coordinates(v,x,h/2)
         call Move Velocities(v,x,h)
         do L=1,Leaps-1
         call Move Coordinates(v,x,h)
         call Move Velocities(v,x,h)
         end do
         call Move Coordinates(v,x,h/2)
         return
         end
         subroutine Evaluate CM Energy(v,x,E)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/T,U,Energy
         save
         T=0
         U=0
         do i=1,3
         i1=3*i-2
         vv=sqa(v(i1))
         rr=sqa(x(i1))
         R=sqrt(rr)
         T=T+mm(i)*vv/(2*mass)
         U=U+mm(i)/R
         end do
         E=T-U
         return
         end
         subroutine Transform to CM system(X,Y)
         implicit real*8 (a-h,m,o-z)
         real*8 Y(9),X(9)
         parameter(I1=0,I2=3,I3=6)
         common/masses/mass,m(3),mu(3),mm(3)
         save
         do k=1,3
         Y(I1+k)=X(I2+k)*mu(3)-X(I3+k)*mu(2)
         Y(I2+k)=X(I3+k)*mu(1)-X(I1+k)*mu(3)
         Y(I3+k)=X(I1+k)*mu(2)-X(I2+k)*mu(1)
         end do 
         return
         end
         subroutine Compute Triangle Differences(X,Y)
         real*8 Y(9),x(9)
         parameter(I1=0,I2=3,I3=6)
         save
         do k=1,3
         Y(I1+k)=X(I3+k)-X(I2+k)
         Y(I2+k)=X(I1+k)-X(I3+k)
         Y(I3+k)=X(I2+k)-X(I1+k)
         end do
         return
         end 
         subroutine Move Coordinates(v,x,s)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/Tk,Up,En
         save
         Tk=0
         do k=1,3
         Tk=Tk+mm(k)*sqa(v(3*k-2))/(2*mass)
         end do
         Te=Tk+v(0)
         dt=s/Te                 ! time increment
         x(0)=x(0)+dt            ! add it 
         do k=1,9+3
         x(k)=x(k)+dt*v(k)       ! coords move linearly
         end do
         return
         end
         function sqa(a)
         real*8 a(3),sqa
         sqa=a(1)**2+a(2)**2+a(3)**2
         return
         end
         subroutine Move Velocities(v,x,s)
         implicit real*8 (a-h,m,o-z)
         parameter(I1=0,I2=3,I3=6)
         real*8 x(0:12),v(0:12),A(9),SumA(3),exacc(0:12),va(9)
         common/masses/mass,m(3),mu(3),mm(3)
         common/diagno/neval
         common/kinpot/Tk,U,En
         save
         neval=neval+1
         R1=sqrt(sqa(x(1)))
         R2=sqrt(sqa(x(4)))
         R3=sqrt(sqa(x(7)))   
         U=mm(1)/R1+mm(2)/R2+mm(3)/R3
         do k=1,3
         A(I1+k)=x(I1+k)/R1**3
         A(I2+k)=x(I2+k)/R2**3
         A(I3+k)=x(I3+k)/R3**3
         SumA(k)=A(I1+k)+A(I2+k)+A(I3+k)
         end do
         call external acceleration(x,exacc) 
         dt=s/U

         do k=1,9
         va(k)=v(k)
         end do ! k

         do k=1,3
         v(I1+k)=v(I1+k)+dt*(-mass*A(I1+k)+m(1)*SumA(k)+exacc(i1+k)) 
         v(I2+k)=v(I2+k)+dt*(-mass*A(I2+k)+m(2)*SumA(k)+exacc(i2+k))  
         v(I3+k)=v(I3+k)+dt*(-mass*A(I3+k)+m(3)*SumA(k)+exacc(i3+k))  
         end do
         do k=1,9
         va(k)=(va(k)+v(k))/2 ! average velocities
         end do
         dTk=0
         do k=1,3
         dTk=dTk+mm(k)*cdot(va(3*k-2),exacc(3*k-2))/(mass)
         end do ! k
          v(0)=v(0)-dt*dTk !binding energy change
         do k=10,12
         v(k)=v(k)+dt*exacc(k) ! cm v jump
         end do
         return
         end
         subroutine external acceleration(x,exacc)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:12),exacc(0:12),xr(9),a(9)
         common/masses/mass,m(3),mu(3),mm(3)
         call  X 2 inertial(x,x(10),xr,time)
         exacc(10)=0
         exacc(11)=0
         exacc(12)=0
         call external coordinate accelerations(xr,a) ! here the ex-acc (=a) is really evaluated
         do i=1,3
         i0=3*i-3
         do k=1,3
         exacc(9+k)=exacc(9+k)+m(i)*a(i0+k)/mass !cm acceleration
         end do! i
         end do! k
         call Compute Triangle Differences(a(1),exacc(1))! side acc
         return
         end
         subroutine external coordinate accelerations(x,a)
         implicit real*8 (a-h,m,o-z)
         real*8 x(9),a(9) 
         common/masses/mass,m(3),mu(3),mm(3) 
c        THIS IS JUST AN EXAMPLE. 
c        a(3*(i-1)+k)= acceration of coordinate x_k of body number I 
c        To remove  external acceleration, set a(k)=0, k=1,9
c        
         do i=1,3
         i0=3*i-3
         i1=i0+1
         rr=cdot(x(i1),x(i1))
         do k=1,3
         a(i0+k)=-x(i0+k)/sqrt(10.d0+rr)**3 -1.d-3*x(i0+k) !  Plummer+(weak) harmonic
         end do
         end do
         return
         end

         subroutine X 2 inertial(x,xcm,xr,time)
         implicit real*8 (a-h,o-z)
         real*8 x(0:9),xcm(3),xr(9)
         time=x(0)
         call Transform to CM system(X(1),XR(1))
         do i=1,3
         i0=3*i-3
         do k=1,3 ! add cm-contributions
         xr(i0+k)=xr(i0+k)+xcm(k) 
         end do !k
         end do !i 
         return
         end
         function cdot(a,b)
         implicit real*8 (a-h,m,o-z)
         real*8 a(3),b(3)
         cdot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
         return
         end



        SUBROUTINE DIFSYAB(N,EPS,S,h,t,Y,Jmax)
        implicit real*8 (a-h,o-z)
c       N=coordin. määrä (=3*NB)
c       F=funktion nimi (FORCE)
        parameter (NMX=125,NMX2=2*NMX,nmx27=nmx2*7) ! NMX=MAX(N),N=3*NB
        REAL*8 Y(N),YR(NMX2),YS(NMX2)
     +  ,DT(NMX2,7),D(7),S(N),EP(4)
        LOGICAL KONV,BO,KL,GR
        DATA EP/.4D-1,.16D-2,.64D-4,.256D-5/
        data dt/nmx27*0.0d0/
        save
        IF(EPS.LT.1.D-14)EPS=1.D-14
        IF(N.gt.NMX)write(6,*) ' too many variables!', char(7)
        if(jmax.lt.4)write(6,*)' too small Jmax (=',jmax,')'
        JTI=0
        FY=1
        redu=0.8d0
        Odot7=0.7
10      tN=t+H
        BO=.FALSE.
C
        M=1
        JR=2
        JS=3
        DO  J=1,Jmax! 10

        do i=1,N
        ys(i)=y(i)
        s(i)=max(abs(ys(i)),s(i)) 
        end do
C

        IF(BO)then
        D(2)=1.777777777777778D0
        D(4)=7.111111111111111D0
        D(6)=2.844444444444444D1
        else
        D(2)=2.25D0
        D(4)=9.D0
        D(6)=36.0D0
        end if

        IF(J.gt.7)then
        L=7
        D(7)=6.4D1
        else
        L=J
        D(L)=M*M
        end if

        KONV=L.GT.3
           subH=H/M
           call SubSteps(YS,subH,M) ! M substeps of size H/M.
        KL=L.LT.2
        GR=L.GT.5
        FS=0.



        DO  I=1,N 
        V=DT(I,1)
        C=YS(I)
        DT(I,1)=C
        TA=C

        IF(.NOT.KL)THEN
        DO  K=2,L
        B1=D(K)*V
        B=B1-C
        W=C-V
        U=V
        if(B.ne.0.0)then
        B=W/B
        U=C*B
        C=B1*B
        end if
        V=DT(I,K)
        DT(I,K)=U
        TA=U+TA
        END DO ! K=2,L
        SI=max(S(I),abs(TA))
        IF(DABS(YR(I)-TA).GT.SI*EPS) KONV=.FALSE.
        IF(.NOT.(GR.OR.SI.EQ.0.D0))THEN
        FV=DABS(W)/SI
        IF(FS.LT.FV)FS=FV
        END IF
        END IF ! .NOT.KL.
        YR(I)=TA
        END DO ! I=1,N

c       end of I-loop
        IF(FS.NE.0.D0)THEN
        FA=FY
        K=L-1
        FY=(EP(K)/FS)**(1.d0/FLOAT(L+K))
        IF(.NOT.((L.NE.2.AND.FY.LT.Odot7*FA).OR.FY.GT.Odot7))THEN
        H=H*FY
               JTI=JTI+1
               IF(JTI.GT.25)THEN
               H=0.0
               RETURN
               END IF
        GO TO 10 ! Try again with a smaller step.
        END IF
        END IF

        IF(KONV)THEN
        t=tN
        H=H*FY
        DO  I=1,N
        Y(I)=YR(I)
        END DO
        RETURN
        END IF

        D(3)=4.D0
        D(5)=1.6D1
        BO=.NOT.BO
        M=JR
        JR=JS
        JS=M+M
        END DO ! J=1,Jmax
        redu=redu*redu+.001d0 ! square the reduction factor (but minimum near 0.001)
        H=H*redu 
        GO TO 10 ! Try again with smaller step.
        END


        SUBROUTINE CONSTANTS OF MOTION(X,V,M,NB,ENERGY,G,Alag)
        IMPLICIT real*8 (A-H,m,O-Z)
        DIMENSION X(*),V(*),G(3),M(*)
        common/kinpot/Tkin,Upot,En
        save
        ee=0
c       Contants of motion in the centre-of-mass system.        
        T=0.0
        U=0.0
        G(1)=0.
        G(2)=0.
        G(3)=0.
        RMIN=1.D30
        mass=M(NB)
        DO I=1,NB-1
        mass=mass+M(I)
        DO J=I+1,Nb
        MIJ=M(I)*M(J)
        KI=(I-1)*3
        KJ=(J-1)*3
        xx=X(KI+1)-X(KJ+1)
        yy=X(KI+2)-X(KJ+2)
        zz=X(KI+3)-X(KJ+3)
        R2=xx*xx+yy*yy+zz*zz+ee
        vx=V(KI+1)-V(KJ+1)
        vy=V(KI+2)-V(KJ+2)
        vz=V(KI+3)-V(KJ+3)
        U=U+MIJ/SQRT(R2)
        T=T+MIJ*(vx*vx+vy*vy+vz*vz)
        G(1)=G(1)+MIJ*(yy*vz-zz*vy)
        G(2)=G(2)+MIJ*(zz*vx-xx*vz)
        G(3)=G(3)+MIJ*(xx*vy-yy*vx)
        END DO
        END DO
        T=T/(2*mass)
        G(1)=G(1)/mass
        G(2)=G(2)/mass
        G(3)=G(3)/mass
        ENERGY=T-U
        Alag=T+U
        Tkin=T ! 
        Upot=U ! 
        En=ENERGY
        RETURN
        END
