module lptsub
    use demoflow
    implicit none
    real(WP), parameter :: dp=0.05_WP
    real(WP), parameter :: rhop=10.0_WP
    real(WP), parameter :: taup=rhop*dp**2.0_WP/mu/18.0_WP
    real(WP), parameter :: e=0.8 ! Coefficient of restitution
    integer :: k,ip,jp
    real(WP) :: fsn,Rep
    real(WP) :: ufp,vfp,xphalf,yphalf,dup,dvp,uphalf,vphalf,duphalf,dvphalf
 ! Spring constant, damping coefficient, force range, distance from/into particle/wall, collision force (x,y), effective mass
    real(WP) :: ksp,eta,lam,dab,delab,fcolx, fcoly, mab
    integer :: ii,jj
    real(WP), dimension(Np,Np*2) :: norm
end module lptsub

subroutine lpt_init
    use lptsub
    implicit none
    integer :: i
    ! Create initial positions/velocites of particles TBD
    do i=1,Np
        xp(i)=.51_WP/i
        yp(i)=0.45_WP/i
        mp(i)=4.0_WP/3.0_WP*pi*(dp/2)**3*rhop
    end do

    up=0.0_WP
    vp=0.4_WP
    liter= 1
    dtp=0.1*dp/sqrt(vp(1)**2+up(1)**2)
    if (dtp.ge.maxdt) dtp=maxdt

end subroutine

subroutine lpt_solve
    use lptsub
    implicit none
    
    ! This time step needs to be matched up with overall flow solver timestep still
    dtp=0.1*dp/sqrt(vp(1)**2+up(1)**2)
    if (dtp.ge.maxdt) dtp=maxdt
    liter=1

    do k=1,Np
        ! ufp, vfp, fsn
        call lpt_fvel(xp(k),yp(k))
        xphalf=xp(k)+dtp/2.0_WP*up(k)
        yphalf=yp(k)+dtp/2.0_WP*vp(k)
        ! Periodic BC
        if (xphalf.ge.Lx-Lx/nx) then
            xp(k)=xp(k)+dtp*up(k)-Lx
            yp(k)=yp(k)+dtp*vp(k)
            ! call lpt_fvel(xp(k),yp(k))
            ! up(k)=ufp
            ! vp(k)=vfp
            up(k)=0.0_WP
            vp(k)=0.0_WP
            cycle
        end if
        call lpt_collisions(xp(k),yp(k),up(k),vp(k))

        ! Calculate dup/dt
        dup=fsn*(ufp-up(k))/taup+gravity(1)+fcolx/mp(k)
        dvp=fsn*(vfp-vp(k))/taup+gravity(2)+fcoly/mp(k)

        
        uphalf=up(k)+dtp/2.0_WP*dup
        vphalf=vp(k)+dtp/2.0_WP*dvp
        xp(k)=xp(k)+dtp*uphalf
        yp(k)=yp(k)+dtp*vphalf
        
        ! Periodic BC
        if (xp(k).ge.Lx-Lx/nx) then
            xp(k)=xp(k)-Lx
            ! call lpt_fvel(xp(k),yp(k))
            ! up(k)=ufp
            ! vp(k)=vfp
            up(k)=0.0_WP
            vp(k)=0.0_WP
            cycle
        end if

        ! new ufp,vfp -> magud
        ! new Rep -> fsn
        call lpt_fvel(xphalf,yphalf)
        ! Calculate dup/dt half
        call lpt_collisions(xphalf,yphalf,uphalf,vphalf)
        duphalf=fsn*(ufp-uphalf)/taup+gravity(1)+fcolx/mp(k)
        dvphalf=fsn*(vfp-vphalf)/taup+gravity(2)+fcoly/mp(k)

        up(k)=up(k)+dtp*duphalf
        vp(k)=vp(k)+dtp*dvphalf
    end do
    
contains
    subroutine lpt_fvel(xpo,ypo)
        implicit none
        integer :: i,j
        real(WP) :: xpo,ypo,stppt,magud,dx2,dy2,xvp,yup,xup,yvp
        stppt=-1
        dx2=(x(2)-x(1))/2
        dy2=(y(2)-y(1))/2
        ! Get surrounding indices
        i=1
        do while (stppt.lt.0)
        stppt=x(i)-xpo
        i=i+1
        end do
        i=i-1
        if ((xpo+dx2).lt.x(i)) then
            xvp=x(i-1)
        else
            xvp=x(i)
        endif

        xup=x(i)-dx2

        stppt=-1
        j=1
        do while (stppt.lt.0)
        stppt=y(j)-ypo
        j=j+1
        end do
        j=j-1

        if ((ypo+dy2).lt.y(j)) then
            yup=y(j-1)
        else
            yup=y(j)
        endif

        yvp=y(j)-dy2
        
        ! interpolate velocities
        ! i=cell right of cell particle is indices
        ufp=1/(x(i)-x(i-1))/(2*dy2)*( &
        &   u(i-1,j-1)*(xup-xpo)*(yup-ypo)+ &
        &   u(i,j-1)*(-(xup-2*dx2)+xpo)*(yup-ypo)+ &
        &   u(i-1,j)*(xup-xpo)*(-(yup-2*dy2)+ypo)+ &
        &   u(i,j)*(-(xup-2*dx2)+xpo)*(-(yup-2*dy2)+ypo))
        
        vfp=1/(2*dx2)/(y(j)-y(j-1))*( &
        &   v(i-1,j-1)*(xvp-xpo)*(yvp-ypo)+ &
        &   v(i,j-1)*(-(xvp-2*dx2)+xpo)*(yvp-ypo)+ &
        &   v(i-1,j)*(xvp-xpo)*(-(yvp-dy2*2)+ypo)+ &
        &   v(i,j)*(-(xvp-2*dx2)+xpo)*(-(yvp-dy2*2)+ypo))

        magud=((ufp-up(k))**2+(vfp-vp(k))**2)**0.5
        Rep=dp*magud/knu
        fsn=1+0.15*Rep**0.687
        ip=i-1
        jp=j-1
        return
    end subroutine lpt_fvel

subroutine lpt_collisions(xpc,ypc,upc,vpc)
    use lptsub
    implicit none
    real(WP) :: xpc,ypc,upc,vpc
    fcolx=0.0_WP
    fcoly=0.0_WP
    lam=0.05_WP*dp
    ksp=mp(k)/((pi**2+log(e)**2)*(100.0_WP*dtp**2))
    eta=-2.0_WP*log(e)*sqrt(mp(k)*ksp/(pi**2+log(e)**2))

    !x force at right wall
    if (mask(ip+1,jp).eq.1) then
        dab=abs(x(ip+1)-xpc)
        delab=abs(dab-dp/2)
        if (dab.lt.(dp/2.0_WP+lam)) then            
            fcolx=-ksp*delab-eta*upc
        end if
    end if

    !x force at left wall
    if (mask(ip-1,jp).eq.1) then
        dab=abs(x(ip)-xpc)
        delab=abs(dab-dp/2)
        if (dab.lt.(dp/2.0_WP+lam)) then
            fcolx=ksp*delab+eta*upc
        end if
    end if

    !y force at top wall
    if (mask(ip,jp+1).eq.1) then
        dab=abs(y(jp+1)-ypc)
        delab=abs((dp/2.0_WP)-dab)
        if (dab.lt.(dp/2.0_WP+lam)) then
            fcoly=-ksp*delab-eta*abs(vpc)
        end if
    end if

    !y force at bottom wall
    if (mask(ip,jp-1).eq.1) then
        dab=abs(y(jp)-ypc)
        delab=abs(dp/2.0_WP-dab)
        if (dab.lt.(dp/2.0_WP+lam)) then
            fcoly=ksp*delab+eta*abs(vpc)
        end if
    end if
    !print *, delab,fcoly/mp(k),y(jp+1),dtp*vpc/dp
    print *, vpc, ypc, y(jp),fcoly/mp(k)*dtp, dab, dp/2


    !particle-particle collisions
    


    !particle-particle collisions

    !get normal vector
    do ii=1,Np
        ! loop over all points
         do jj=1,Np
            norm(ii,2*jj-1) = (xp(jj)-xpc)/(sqrt((xpc-xp(jj))**2+(ypc-yp(k))**2)) ! x-normal of particle ii with particle jj
            norm(ii,2*jj  ) = (yp(jj)-ypc)/(sqrt((xpc-xp(jj))**2+(ypc-yp(k))**2)) ! y-normal of particle ii with particle jj
        end do
    end do

!print *, norm(2,1),xpc
end subroutine lpt_collisions

end subroutine lpt_solve