         implicit real*8 (a-h,m,o-z)
         real*8 xr(0:9),vr(0:9),dxr(0:9),dvr(0:9)
         real*8 g0(3),gt(3),dg0(3),dgt(3)
         common/diagno/neval
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         open(67,file='xyz.L')
         eval=0.0
         read(5,*)tmx,tincr,eps ! read step & max time
         do i=1,3
         k1=i*3-2
         k2=k1+2
         read(5,*)m(i),(xr(k),k=k1,k2),(vr(k),k=k1,k2) ! read mass & initial coord/vels
         end do
         do i=0,9
         dxr(i)=0
         dvr(i)=0
         end do
         dxr(1)=1
         drlog=0
         call Reduce2cm(xr,m,3)
         call Reduce2cm(vr,m,3)
         call Reduce2cm(dxr,m,3)
         call Reduce2cm(dvr,m,3)
         nb=3
         newcase=1
         xr(0)=0 !==time
         call Constants Of Motion
     &   (xr,vr,dxr,dvr,m,nb,ene0,dene0,g0,dg0,alag0)
1        continue         
         call ThreeB(m,xr,vr,dxr,dvr,drlog,tincr,newcase,eps) 
         call Constants Of Motion
     &   (xr,vr,dxr,dvr,m,nb,enet,denet,gt,dgt,alag)
            ss=0
         do i=0,9
         ss=ss+abs(dxr(i))+abs(dvr(i))
         end do
           sl=log(ss)
         vrh=abs(enet-ene0)/alag ! measure of energy error
         glag=mass**2.5/abs(alag0)**0.5!1.5*(alag)
         dgg=(gt(1)-g0(1))**2+(gt(2)-g0(2))**2+(gt(3)-g0(3))**2
         dgg=sqrt(dgg)/glag      !measure of ang.mom. error
         time=xr(0)
         drw=drlog+sl
         alam=drw/(time+1.e-6)
         write(6,106)time,vrh,dgg,drw,alam


106      format(1x,' T: ',f12.4,1p,2g11.2,6g12.4)
         write(67,123)XR
         write(68,123)xr(0),dxr(0)-1,dxr
123      format(1x,f12.4,1p,100g12.4)
         call flush(6)
         call flush(66)
         if(time.lt.tmx)goto 1
         end
         subroutine  re_scale(dx,dv,drlog)
         implicit real*8 (a-h,o-z)
         real*8 dx(0:9),dv(0:9)
         s=0
         do k=0,9
         s=s+abs(dx(k))+abs(dv(k))
         end do 
         drlog=drlog+log(s)
         do k=0,9
         dx(k)=dx(k)/s
         dv(k)=dv(k)/s
         end do         
         return
         end
         subroutine ThreeB(mr,xr,vr,dxr,dvr,drlog,tincr,newcase,eps)
         implicit real*8 (a-h,m,o-z)
         real*8 xr(0:9),vr(0:9),x(0:9),v(0:9),y(40),SY(40),mr(3)
         real*8 dxr(0:9),dvr(0:9),dx(0:9),dv(0:9)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         equivalence(y(1),x(0)),(y(11),v(0)),(y(21),dx(0)),(y(31),dv(0))
         data sy/40*0.0/
         save
c        initializations         
         time=xr(0)! ???
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
         call Compute Triangle Differences(xr,x)
         call Compute Triangle Differences(vr,v)
         call Compute Triangle Differences(dxr,dx)
         call Compute Triangle Differences(dvr,dv)

         call Transform to CM system(X,XR)
         call Transform to CM system(V,VR)
         call Transform to CM system(dx,dxr)
         call Transform to CM system(dv,dvr)
         call Evaluate Energy(v,x,dv,dx,E,dE)
         h=.1/(TKIN+UPOT) 
         Energy=E ! to common
         v(0)=-E ! binding energy
         x(0)=0  ! time
         dv(0)=-dE
         dx(0)=dxr(0)!0
         end if ! newcase
c         
         told=x(0)
         nsteps=0
         NVX=20*2
         nzero=0
1        continue
c        order of magnitude of variables
         rsmall=mmm/UPOT
         vsmall=sqrt(TKIN/mass)
         tsmall=rsmall*sqrt(rsmall/mass)
         esmall=abs(v(0))  
         sy(1)=tsmall*10 +abs(time)
         sy(11)=esmall
         call orderY(Y,SY)
         hold=h
         if(tincr.ne.0.0)h=min(h,tincr*Upot)
         JMAX=10 ! =DIFSY order= h^(2*JMAX) if necessary
         call  DIFSYAB(NVX,EPS,SY,h,dummy,Y,JMAX) ! best??
         do kmet=0,1
         if(h.eq.0.0)then
         h=hold
         call gdfsy(nvx,h,eps,SY,dummy,y,kmet) ! if h=0, try two other methods
         end if
         end do ! kmet
         call Triangle Correction(v,x)
         call Triangle Correction(dv,dx)
c-------------------------------------try to change the step (if zero)
            if(h.eq.0.0)then
                 nzero=nzero+1
                 write(6,*)' STEP=0!',char(7)
              if(nzero.gt.10)then
                 STOP
                 else
                  if(nzero.ne.nzero/2*2)then
                  h=hold*1.1d0  ! try increase
                  else
                  h=hold*.1d0   ! try decrease
                  end if
              end if
                 else
                 nzero=0
            end if
c--------------------------------------------------------
         nsteps=nsteps+1
         if(x(0).lt.told+tincr.and.nsteps.lt.10000000)goto 1
         time=x(0)
         xr(0)=time
         call Transform to CM system(X,XR)
         call Transform to CM system(V,VR)
         call re_scale(dx,dv,drlog)
         call Transform to CM system(dX,dXR)
         call Transform to CM system(dV,dVR)
         return
         end
         subroutine orderY(Y,SY)
         real*8 Y(*),SY(*)
         sy(2)=abs(y(2))+abs(y(3))+abs(y(4))
         sy(3)=sy(2)
         sy(4)=sy(3)
         sy(5)=abs(y(5))+abs(y(6))+abs(y(7))
         sy(6)=sy(5)
         sy(7)=sy(6)
         sy(8)=abs(y(8))+abs(y(9))+abs(y(10))
         sy(9)=sy(8)
         sy(10)=sy(9)
         sy(12)=abs(y(12))+abs(y(13))+abs(y(14))
         sy(13)=sy(12)
         sy(14)=sy(13)
         sy(15)=abs(y(15))+abs(y(16))+abs(y(17))
         sy(16)=sy(15)
         sy(17)=sy(16)
         sy(18)=abs(y(18))+abs(y(19))+abs(y(20))
         sy(19)=sy(18)
         sy(20)=sy(19)
         do i=21,40
         sy(i)=1+abs(y(i))*100 ! for variations
         end do
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
         call LeapFrog(y(11),y(1),y(31),y(21),H,Leaps)
         return
         end
         subroutine LeapFrog(v,x,dv,dx,h,Leaps)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9),dx(0:9),dv(0:9)
         save
         call Move Coordinates(v,x,dv,dx,h/2)
         call Move Velocities(v,x,dv,dx,h)
         do L=1,Leaps-1
         call Move Coordinates(v,x,dv,dx,h)
         call Move Velocities(v,x,dv,dx,h)
         end do
         call Move Coordinates(v,x,dv,dx,h/2)
         return
         end
         subroutine Evaluate Energy(v,x,dv,dx,E,dE)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9),dx(0:9),dv(0:9)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/T,U,Energy
         save
         T=0
         U=0
         dT=0
         dU=0
         do i=1,3
         i1=3*i-2
         vv=sqa(v(i1))
         dvv=2*cdot(v(i1),dv(i1))
         rr=sqa(x(i1))
         drr=2*cdot(x(i1),dx(i1))
         R=sqrt(rr)
         dR=1/(2*R)*drr
         T=T+mm(i)/(2*mass)*vv
         dT=dT+mm(i)/(2*mass)*dvv
         U=U+mm(i)/R
         dU=dU-mm(i)/rr*dR
         end do
         E=T-U
         dE=dT-dU
         return
         end
         subroutine Transform to CM system(X,Y)
         implicit real*8 (a-h,m,o-z)
         real*8 Y(0:9),X(0:9)
         parameter(I1=0,I2=3,I3=6)
         common/masses/mass,m(3),mu(3),mm(3)
         save
         Y(0)=X(0)
         do k=1,3
         Y(I1+k)=X(I2+k)*mu(3)-X(I3+k)*mu(2)
         Y(I2+k)=X(I3+k)*mu(1)-X(I1+k)*mu(3)
         Y(I3+k)=X(I1+k)*mu(2)-X(I2+k)*mu(1)
         end do 
         return
         end
         subroutine Compute Triangle Differences(X,Y)
         real*8 Y(0:9),x(0:9)
         parameter(I1=0,I2=3,I3=6)
         save
         Y(0)=X(0)
         do k=1,3
         Y(I1+k)=X(I3+k)-X(I2+k)
         Y(I2+k)=X(I1+k)-X(I3+k)
         Y(I3+k)=X(I2+k)-X(I1+k)
         end do
         return
         end 
         subroutine Move Coordinates(v,x,dv,dx,s)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9),dx(0:3),dv(0:3)
         common/masses/mass,m(3),mu(3),mm(3)
         save
         Te=v(0)                            !binding energy
         dTe=dv(0)
         do k=1,3
         Te=Te+mm(k)/(2*mass)*sqa(v(3*k-2)) ! kinetic energy
         dTe=dTe+mm(k)/(2*mass)*2*cdot(v(3*k-2),dv(3*k-2))
         end do
         dt=s/Te                 ! time increment
         ddt=-s/Te**2*dTe
         x(0)=x(0)+dt            ! add it 
         dx(0)=dx(0)+ddt
         do k=1,9
         x(k)=x(k)+dt*v(k)       ! coords move linearly
         dx(k)=dx(k)+ddt*v(k)+dt*dv(k)
         end do
         return
         end
         subroutine Move Velocities(v,x,dv,dx,s)
         implicit real*8 (a-h,m,o-z)
         parameter(I1=0,I2=3,I3=6)
         real*8 x(0:9),v(0:9),A(9),SumA(3)
         real*8 dx(0:9),dv(0:9),dA(9),dSumA(3)
         common/masses/mass,m(3),mu(3),mm(3)
         common/diagno/neval
         save
         neval=neval+1
         R1=sqrt(sqa(x(1))) ! triangle side1
         dR1=cdot(x(1),dx(1))/R1
         R2=sqrt(sqa(x(4))) !          side2
         dR2=cdot(x(4),dx(4))/R2
         R3=sqrt(sqa(x(7))) !          side3  
         dR3=cdot(x(7),dx(7))/R3
         U=mm(1)/R1+mm(2)/R2+mm(3)/R3 ! potential
         dU=-mm(1)/R1**2*dR1-mm(2)/R2**2*dR2-mm(3)/R3**2*dR3
         do k=1,3  ! auxiliary quantities for acceleration computation
         A(I1+k)=x(I1+k)/R1**3 
         dA(i1+k)=dx(I1+k)/R1**3-3*x(i1+k)/R1**4*dR1
         A(I2+k)=x(I2+k)/R2**3
         dA(I2+k)=dx(i2+k)/R2**3-3*x(i2+k)/R2**4*dR2
         A(I3+k)=x(I3+k)/R3**3
         dA(I3+k)=dx(i3+k)/R3**3-3*x(i3+k)/R3**4*dR3
         SumA(k)=A(I1+k)+A(I2+k)+A(I3+k)
         dSumA(k)=dA(I1+k)+dA(I2+k)+dA(I3+k)
         end do
         
         dt=s/U     ! time increment (only for moving v's)
         ddt=-s/U**2*dU
         do k=1,3           !  accelerations of triangle sides
         v(I1+k)=v(I1+k)+dt*(-mass*A(I1+k)+m(1)*SumA(k) ) 
         dv(I1+k)=dv(I1+k)+ddt*(-mass*A(I1+k)+m(1)*SumA(k) )
     &   +dt*(-mass*dA(I1+k)+m(1)*dSumA(k) )
         v(I2+k)=v(I2+k)+dt*(-mass*A(I2+k)+m(2)*SumA(k) )  
         dv(I2+k)=dv(I2+k)+ddt*(-mass*A(I2+k)+m(2)*SumA(k) )
     &   +dt*(-mass*dA(I2+k)+m(2)*dSumA(k) )
         v(I3+k)=v(I3+k)+dt*(-mass*A(I3+k)+m(3)*SumA(k) )  
         dv(I3+k)=dv(I3+k)+ddt*(-mass*A(I3+k)+m(3)*SumA(k) )
     &   +dt*(-mass*dA(I3+k)+m(3)*dSumA(k) )
         end do
         return
         end
         function cdot(a,b)
         implicit real*8 (a-h,m,o-z)
         real*8 a(3),b(3)
         cdot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
         return
         end
         function sqa(a)
         real*8 a(3),sqa
         sqa=a(1)**2+a(2)**2+a(3)**2
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


        subroutine constants of motion
     &  (x,v,dx,dv,m,nb,energy,denergy,g,dg,alag)
        implicit real*8 (a-h,m,o-z)
        dimension x(0:9),v(0:9),m(3),g(3),dx(0:9),dv(0:9),dg(3)
        save
        T=0
        U=0
        dT=0
        dU=0
        do k=1,3
        g(k)=0
        dg(k)=0
        end do
        rmin=1.d30

        do 10 i=1,3
        k1=(i-1)*3+1
        k2=k1+1
        k3=k2+1
        T=T+.5d0*m(i)*cdot(v(k1),v(k1))
        dT=dT+m(i)*cdot(v(k1),dv(k1))
        g(1)=g(1)+m(i)*(x(k2)*v(k3)-x(k3)*v(k2))
        dg(1)=dg(1)+m(i)*
     &  (dx(k2)*v(k3)-dx(k3)*v(k2)+x(k2)*dv(k3)-x(k3)*dv(k2))
        g(2)=g(2)-m(i)*(x(k1)*v(k3)-x(k3)*v(k1))
        dg(2)=dg(2)-m(i)*
     &  (dx(k1)*v(k3)-dx(k3)*v(k1)+x(k1)*dv(k3)-x(k3)*dv(k1))
        g(3)=g(3)+m(i)*(x(k1)*v(k2)-x(k2)*v(k1))
        dg(3)=dg(3)-m(i)*
     &  (dx(k1)*v(k2)-dx(k2)*v(k1)+x(k1)*dv(k2)-x(k2)*dv(k1))
        if(i.eq.nb)go to 10
        j1=i+1
        do  j=j1,3
        ki=(i-1)*3
        kj=(j-1)*3
        r2=0.
        dr2=0
        do  k=1,3
        ki=ki+1
        kj=kj+1
        r2=r2+(x(ki)-x(kj))**2
        dr2=dr2+2*(x(ki)-x(kj))*(dx(ki)-dx(kj))
        end do
        r1=sqrt(r2)
        u=u+m(i)*m(j)/r1
        dU=dU-m(i)*m(j)*dr2/(r1*r2)/2
        end do
10      continue
        energy=T-U
        denergy=dT-dU
        alag=t+u
        return
        end

        subroutine gdfsy(nv,h,eps,yscal,t,y,imet)
c	y dependent variable
c	nv =legth of y
c	t =independent variable
c	eps =accuracy
c	yscal(i) =order of y(i)
       implicit real*8 (a-h,o-z)
       parameter(nmax=100,imax=11,nuse=7,one=1.d0,
     -  shrink=.95d0,grow=1.1d0)
        dimension y(nv),yscal(nv),yerr(nmax)
     -  ,yseq(nmax),nseq(imax),yextr(nmax)
         data nseq/2,4,6,8,10,11,12,13,14,15,16/
c          data nseq/1,2,3,4,5,6,7,8,9,10,11/
         save
       tsav=t
1       continue
        do 10 i=1,imax
       do  j=1,nv
       yseq(j)=y(j)
        end do

        call substeps(yseq,h/nseq(i),nseq(i))

        test=(h/nseq(i))**2
         if(imet.eq.0)then
       call rzextr(i,test,yseq,yextr,yerr,nv,nuse)
         else
       call pzextr(i,test,yseq,yextr,yerr,nv,nuse)
         end if
         conv=0
        do 12 j=1,nv
         SI=max(abs(yextr(j)),yscal(j))
        if(abs(yerr(j)).gt.eps*SI)then
          conv=2.
          goto 13
          end if
        IF(i.lt.nuse-3)conv=2.
12       continue
13       if(conv.lt.one)then
       t=t+h
       if(i.gt.nuse)then
               h=h*shrink
       elseif(i.eq.nuse-1)then
               h=h*grow
       else
               h=(h*nseq(nuse-1))/nseq(i)
        end if
        do j=1,nv
        y(j)=yextr(j)! return the extrapolated results
        end do       
       return
       end if
10       continue
c       if it gets here, then step failed. Reduce step and try again.
       h=.25*h/2**((imax-nuse)/2)
       goto 1
       end
      subroutine pzextr(iest,xest,yest,yz,dy,nv,nuse)
      implicit real*8 (a-h,o-z)
      parameter(imax=11,ncol=7,nmax=100)
      real*8 x(imax),yest(nv),yz(nv),dy(nv),qcol(nmax,ncol),d(nmax)
      save
      x(iest)=xest ! store current dependent value
      do 11 j=1,nv
      dy(j)=yest(j)
      yz(j)=yest(j)
11    continue
      if(iest.eq.1)then
      do 12 j=1,nv
       qcol(j,1)=yest(j)
12    continue
      else
      m1=min(iest,nuse)
      do 13 j=1,nv
      d(j)=yest(j)
13    continue
      do 15 k1=1,m1-1
      delta=1./(x(iest-k1)-xest)
      f1=xest*delta
      f2=x(iest-k1)*delta
      do 14 j=1,nv
      q=qcol(j,k1)
      qcol(j,k1)=dy(j)
      delta=d(j)-q
      dy(j)=f1*delta
      d(j)=f2*delta
      yz(j)=yz(j)+dy(j)
14    continue
15    continue
      do 16 j=1,nv
      qcol(j,m1)=dy(j)
16    continue
      end if
      return
      end

       subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse)
       implicit real*8 (a-h,o-z)
       parameter(imax=11,nmax=100,ncol=9)
       dimension x(imax),yest(nv),yz(nv),dy(nv),d(nmax,ncol),fx(ncol)
       save
       x(iest)=xest
       if(iest.eq.1) then
         do 11 j=1,nv
         yz(j)=yest(j)
         d(j,1)=yest(j)
         dy(j)=yest(j)
11         continue
       else
        m1=min(iest,nuse)
        do 12 k=1,m1-1
        fx(k+1)=x(iest-k)/xest
12        continue
       do 14 j=1,nv
       yy=yest(j)
       v=d(j,1)
       c=yy
       d(j,1)=yy
       do 13 k=2,m1
       b1=fx(k)*v
       b=b1-c
          if(b.ne.0.)then
       b=(c-v)/b
       ddy=c*b
       c=b1*b
              else
       ddy=v
              endif
       if(k.ne.m1)v=d(j,k)
       d(j,k)=ddy
       yy=yy+ddy
13       continue
       dy(j)=ddy
       yz(j)=yy
14       continue
       end if
       return
       end

        subroutine Reduce2cm(x,m,nb)
        real*8 x(0:9),m(3),cm(3)
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
        do i=1,3
        do k=1,3
        x(k+3*(i-1))=x(k+3*(i-1))-cm(k)/sm
        end do
        end do
        return
        end
