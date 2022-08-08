         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9)
         real*8 g0(3),gt(3)
         real*8 xcm(3),vcm(3) 
         common/diagno/neval
         common/soft/soft2
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         common/DIFSYorder/Jmax
         eval=0.0
         open(66,file='xyz')! outputfile for t,x,y,z,vx,vy,vz
         read(5,*)tmx,tincr,eps,icm,Jmax ,soft1
         soft2=soft1**2
c                   tmx=max time
c                 tincr=output interval in physical time (approximate only)
c            eps= error tolerance
c                  icm= cm-index (>0=>reduce2cm)
c        Jmax determines integration order, use 4 for short steps, 10 for long steps
c         soft is softening parameter. Usually =0, it is here obnly for testing the effect.
         do i=1,3
         k1=i*3-2
         k2=k1+2
         read(5,*)m(i),(xr(k),k=k1,k2),(vr(k),k=k1,k2) ! read mass & initial coord/vels
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

         write(66,123)time,XR,VR!write (inertial system) coords and vels

         call ThreeB(m,time,xr,vr,tincr,newcase,eps) 
c        begin{diagnostics}
         call Constants Of Motion(xr,vr,m,nb,enet,gt,alag)
         vrh=abs(enet-ene0)/upot ! ~ LogH=the primary integral (should be =0) 
         glag=mass**2.5/abs(enet)**2*sqrt(alag)**3
         dgg=(gt(1)-g0(1))**2+(gt(2)-g0(2))**2+(gt(3)-g0(3))**2
         dgg=sqrt(dgg)/glag ! ang.mom. 'relative' error
         write(6,106)time,vrh,dgg,neval,enet
106      format(1x,' T: ',f12.4,1p,2g11.2,i9,1p,g20.12)
123      format(1x,f12.4,1p,19g14.6)
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

         subroutine ThreeB(mr,time,xr,vr,tincr,newcase,eps)
         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9),x(0:9),v(0:9),y(20),SY(20),mr(3)
         real*8 xcm(3),vcm(3)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         equivalence(y(1),x(0)),(y(11),v(0))
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
         call reduce 2 cm(xr,m,xcm)
         call reduce 2 cm(vr,m,vcm)
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
         NVX=20
1        continue
c        order of magnitude of variables
         rsmall=mmm/UPOT
         vsmall=sqrt(TKIN/mass)
         tsmall=rsmall*sqrt(rsmall/mass)
         esmall=abs(v(0))  
         sy(1)=tsmall/10
         sy(11)=esmall/10
         do k=2,10
         sy(k)=rsmall/10
         sy(k+10)=vsmall/10
         end do
c
         h=min(h,tincr*Upot) !timestep at most = tincr (approx)
         call  DIFSYAB(NVX,EPS,SY,h,dummy,Y) ! best??
         call Triangle Correction(v,x)
                 if(h.eq.0.0)then
                 write(6,*)' STEP=0!',char(7)
                 STOP
                 end if
         nsteps=nsteps+1
         if(x(0).lt.told+tincr.and.nsteps.lt.1000000)goto 1

         time=x(0)
         call Transform to CM system(X(1),XR(1))
         call Transform to CM system(V(1),VR(1))
         do i=1,3
         i0=3*i-3
         do k=1,3 ! add cm-contributions
         xr(i0+k)=xr(i0+k)+xcm(k)+x(0)*vcm(k) 
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
         real*8 y(20)
         save
         call LeapFrog(y(11),y(1),H,Leaps)
         return
         end
         subroutine LeapFrog(v,x,h,Leaps)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9)!,xx(0:9),vv(0:9)
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
         common/soft/soft2
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/T,U,Energy
         save
         T=0
         U=0
         do i=1,3
         i1=3*i-2
         vv=sqa(v(i1))
         rr=sqa(x(i1)) + soft2
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
         do k=1,9
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
         real*8 x(0:9),v(0:9),A(9),SumA(3)
         common/masses/mass,m(3),mu(3),mm(3)
         common/diagno/neval
         common/soft/soft2
         common/kinpot/Tk,U,En
         save
         neval=neval+1
         R1=sqrt(sqa(x(1))+soft2)
         R2=sqrt(sqa(x(4))+soft2)
         R3=sqrt(sqa(x(7))+soft2)   
         U=mm(1)/R1+mm(2)/R2+mm(3)/R3
         do k=1,3
         A(I1+k)=x(I1+k)/R1**3
         A(I2+k)=x(I2+k)/R2**3
         A(I3+k)=x(I3+k)/R3**3
         SumA(k)=A(I1+k)+A(I2+k)+A(I3+k)
         end do
         
         dt=s/U
         do k=1,3
         v(I1+k)=v(I1+k)+dt*(-mass*A(I1+k)+m(1)*SumA(k) ) 
         v(I2+k)=v(I2+k)+dt*(-mass*A(I2+k)+m(2)*SumA(k) )  
         v(I3+k)=v(I3+k)+dt*(-mass*A(I3+k)+m(3)*SumA(k) )  
         end do
         return
         end
         function cdot(a,b)
         implicit real*8 (a-h,m,o-z)
         real*8 a(3),b(3)
         cdot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
         return
         end


        SUBROUTINE DIFSYAB(N,EPS,S,h,t,Y)! Burlirch-Stoer extrapolation of the results of subroutine SubSteps
        implicit real*8 (a-h,o-z)
c       N=number of variables  (=3*NB)
c       EPS error tolerance (recommendation =1.e-13 for high accuracy)
c       S = vector of variable size estimate for the test |dY|<EPS*S
c       T= independent variable (e.g. time or whateve U prefer)
c       Y=the vctor of integrated variables. Lenth at least N.        
        parameter (NMX=25,NMX2=2*NMX,nmx27=nmx2*7) ! NMX=MAX(N),N=3*NB
        REAL*8 Y(N),YR(NMX2),YS(NMX2)
     +  ,DT(NMX2,7),D(7),S(N),EP(4)
        LOGICAL KONV,BO,KL,GR
         common/DIFSYorder/Jmax
        DATA EP/.4D-1,.16D-2,.64D-4,.256D-5/
        data dt/nmx27*0.0d0/
        save
        IF(EPS.LT.1.D-14)EPS=1.D-14 ! max accuracy, may be better to use at most 1.e-13
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

        subroutine constants of motion(x,xd,m,nb,energy,g,alag)
        implicit real*8 (a-h,m,o-z)
        dimension x(*),xd(*),m(*),g(3)
        common/kinpot/TKIN,UPOT,Enrg
        common/soft/soft2
        save
        t=0.0
        ug=0.0
        u=ug
        g(1)=0.
        g(2)=0.
        g(3)=0.
        rmin=1.d30
        do 10 i=1,nb
        k1=(i-1)*3+1
        k2=k1+1
        k3=k2+1
        t=t+.5d0*m(i)*(xd(k1)**2+xd(k2)**2+xd(k3)**2)
        g(1)=g(1)+m(i)*(x(k2)*xd(k3)-x(k3)*xd(k2))
        g(2)=g(2)-m(i)*(x(k1)*xd(k3)-x(k3)*xd(k1))
        g(3)=g(3)+m(i)*(x(k1)*xd(k2)-x(k2)*xd(k1))
        if(i.eq.nb)go to 10
        j1=i+1
        do 9 j=j1,nb
        ki=(i-1)*3
        kj=(j-1)*3
        r2=0.+soft2
        do 8 k=1,3
        ki=ki+1
        kj=kj+1
8       r2=r2+(x(ki)-x(kj))**2
        u=u+m(i)*m(j)/sqrt(r2)
9       continue
10      continue
        energy=t-u
        alag=t+u
               tkin=t
               upot=u
               enrg=energy
        return
        end
