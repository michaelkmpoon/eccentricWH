         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9)
         real*8 g0(3),gt(3)
         common/diagno/neval
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
         eval=0.0
         read(5,*)eta,tmx,tincr ! read accu, max-time, output-interval
         open(66,file='xyz.s')
         do i=1,3
         k1=i*3-2
         k2=k1+2
         read(5,*)m(i),(xr(k),k=k1,k2),(vr(k),k=k1,k2) ! read mass & initial coord/vels
         end do
         call reduce2cm(xr,m) !only cm-system possible
         call reduce2cm(vr,m) !  -"-
         nb=3
         newcase=1
         time=0
         call Constants Of Motion(xr,vr,m,nb,ene0,g0,alag)
1        continue         
         call ThreeB(eta,m,time,xr,vr,tincr,newcase,eps) 
         call Constants Of Motion(xr,vr,m,nb,enet,gt,alag)
         upot=(alag-enet)/2
         vrh=abs(enet-ene0)/upot
         glag=mass**2.5/abs(ene0)*sqrt(alag)
         dgg=(gt(1)-g0(1))**2+(gt(2)-g0(2))**2+(gt(3)-g0(3))**2
         dgg=sqrt(dgg)/upot
         xscale=(m(1)*(m(2)+m(3))+m(2)*m(3))/abs(ene0)
         write(6,106)time,vrh,dgg,enet-ene0,neval
106      format(1x,' T: ',f12.4,1p,3g10.1,I12,6g12.4)
         write(66,123)time,XR
123      format(1x,f12.4,1p,9g12.4)
         call flush(6)
         call flush(66)
         if(time.lt.tmx)goto 1
         write(6,*)' x_scale ',xscale
         end

         subroutine ThreeB(eta,mr,time,xr,vr,tincr,newcase,eps)
         implicit real*8 (a-h,m,o-z)
         real*8 xr(9),vr(9),x(0:9),v(0:9),mr(3)
         common/masses/mass,m(3),mu(3),mm(3)
         common/kinpot/TKIN,UPOT,Energy
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
         call Compute Triangle Differences(xr(1),x(1))
         call Compute Triangle Differences(vr(1),v(1))
         call Transform to CM system(X(1),XR(1))
         call Evaluate Energy(v,x,E)
         Energy=E ! to common
         v(0)=-E ! binding energy
         x(0)=0  ! time
         if(energy.lt.0.0)then
         h=eta*mmm*sqrt(mass/abs(Energy))!this scales correctly
         else
         h=eta*mmm*sqrt(mass/(tkin+energy))
         end if !E<0.
         end if ! newcase
c         
         told=x(0)
         nsteps=0
1        continue
c        order of magnitude of variables
         tnext=told+tincr
         call LeapFrog(v,x,h,tnext)
c
         call Triangle Correction(v,x)
                 if(h.eq.0.0)then
                 write(6,*)' STEP=0!',char(7)
                 STOP
                 end if
         nsteps=nsteps+1
         if(x(0).lt.told+tincr.and.nsteps.lt.10000000)goto 1
         time=x(0)
         call Transform to CM system(X(1),XR(1))
         call Transform to CM system(V(1),VR(1))

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
         
         subroutine LeapFrog(v,x,h,tnext)
         implicit real*8 (a-h,m,o-z)
         real*8 x(0:9),v(0:9)!,xx(0:9),vv(0:9)
         save
         Leaps=1000! just something to prevent an infinite loop
         call Move Coordinates(v,x,h/2) ! first halfstep
         call Move Velocities(v,x,h)
         do L=1,Leaps-1
         call Move Coordinates(v,x,h)
         call Move Velocities(v,x,h)
         if(x(0).gt.tnext)goto 22 ! 'output' when next time reached
         end do
22       call Move Coordinates(v,x,h/2) !last halfstep
         return
         end
         subroutine Evaluate Energy(v,x,E)
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
         save
         Te=v(0)
         do k=1,3
         Te=Te+mm(k)*sqa(v(3*k-2))/(2*mass)
         end do
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


        subroutine constants of motion(x,xd,m,nb,energy,g,alag)
        implicit real*8 (a-h,m,o-z)
        dimension x(*),xd(*),m(*),g(3)
        common/kinpot/t,u,enedummy
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
        r2=0.
        do 8 k=1,3
        ki=ki+1
        kj=kj+1
8       r2=r2+(x(ki)-x(kj))**2
        u=u+m(i)*m(j)/sqrt(r2)
9       continue
10      continue
        energy=t-u
        alag=t+u
        return
        end

        subroutine reduce2cm(x,m)
        implicit real*8 (a-h,m,o-z)
        real*8 cm(3),x(9),m(3)
        mass=m(1)+m(2)+m(3)
        do k=1,3
        cm(k)=0
        do i=1,3
        i0=3*i-3
        cm(k)=cm(k)+m(i)*x(i0+k) 
        end do
        end do

        do i=1,3
        i0=3*i-3
        do k=1,3
        x(i0+k)=x(i0+k)-cm(k)/mass
        end do
        end do
        return
        end

