module lptsub
    use demoflow
    implicit none
    real(WP), parameter :: dp=0.1_WP
    real(WP), parameter :: rhop=30.0_WP
    real(WP), parameter :: taup=rhop*dp**2.0_WP/mu/18.0_WP
    real(WP), parameter :: e=0.8 ! Coefficient of restitution
    integer :: k,ip,jp
    real(WP) :: fsn,Rep,dtp
    real(WP) :: ufp,vfp,xphalf,yphalf,dup,dvp,uphalf,vphalf,duphalf,dvphalf
 ! Spring constant, damping coefficient, force range, distance into particle/wall, collision force (x,y), effective mass
    real(WP) :: ksp,eta,lam,dab,fcolx, fcoly, mab
end module lptsub

subroutine lpt_init
    use lptsub
    implicit none
    integer :: i
    ! Create initial positions/velocites of particles TBD
    do i=1,Np
        xp(i)=.5_WP/i
        yp(i)=0.3_WP/i
        mp(i)=4.0_WP/3.0_WP*pi*(dp/2)**3*rhop
    end do

    up=0.0_WP
    vp=0.0_WP
    liter= 10

end subroutine

subroutine lpt_solve
    use lptsub
    implicit none
    
    ! This time step needs to be matched up with overall flow solver timestep still
    dtp=0.1*dp/maxval(sqrt(vp**2+up**2))

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

        ! Calculate dup/dt
        dup=fsn*(ufp-up(k))/taup+gravity(1)
        dvp=fsn*(vfp-vp(k))/taup+gravity(2)

        
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
        duphalf=fsn*(ufp-uphalf)/taup+gravity(1)
        dvphalf=fsn*(vfp-vphalf)/taup+gravity(2)

        up(k)=up(k)+dtp*duphalf
        vp(k)=vp(k)+dtp*dvphalf
        
    end do
    
contains
    subroutine lpt_fvel(xpo,ypo)
        implicit none
        integer :: i,j
        real(WP) :: xpo,ypo,stppt,magud,dx2,dy2,xvp,yup,xup,yvp
        stppt=-1
        dx2=(xm(2)-xm(1))/2
        dy2=(ym(2)-ym(1))/2
        ! Get surrounding indices
        i=1
        do while (stppt.lt.0)
        stppt=xm(i)-xpo
        i=i+1
        end do

        if ((xpo+dx2).lt.xm(i)) then
            xvp=xm(i-1)
        else
            xvp=xm(i)
        endif

        xup=xm(i)-dx2

        stppt=-1
        j=1
        do while (stppt.lt.0)
        stppt=ym(j)-ypo
        j=j+1
        end do

        if ((ypo+dy2).lt.ym(j)) then
            yup=ym(j-1)
        else
            yup=ym(j)
        endif

        yvp=ym(j)-dy2
        
        ! interpolate velocities
        ! i=cell right of cell particle is indices
        ufp=1/(xm(i)-xm(i-1))/(2*dy2)*( &
        &   u(i-1,j-1)*(xup-xpo)*(yup-ypo)+ &
        &   u(i,j-1)*(-(xup-2*dx2)+xpo)*(yup-ypo)+ &
        &   u(i-1,j)*(xup-xpo)*(-(yup-2*dy2)+ypo)+ &
        &   u(i,j)*(-(xup-2*dx2)+xpo)*(-(yup-2*dy2)+ypo))
        
        vfp=1/(2*dx2)/(ym(j)-ym(j-1))*( &
        &   v(i-1,j-1)*(xvp-xpo)*(yvp-ypo)+ &
        &   v(i,j-1)*(-(xvp-2*dx2)+xpo)*(yvp-ypo)+ &
        &   v(i-1,j)*(xvp-xpo)*(-(yvp-dy2*2)+ypo)+ &
        &   v(i,j)*(-(xvp-2*dx2)+xpo)*(-(yvp-dy2*2)+ypo))

        magud=((ufp-up(k))**2+(vfp-vp(k))**2)**0.5
        Rep=dp*magud/knu
        fsn=1+0.15*Rep**0.687
        ip=i
        jp=j
        return
    end subroutine lpt_fvel
end subroutine lpt_solve

subroutine lpt_collisions
    use lptsub
    implicit none
    fcolx=0
    fcoly=0
    lam=0.05*dp
    ksp=mp(k)/(pi**2+log(e)**2)/(100*dtp**2)
    eta=-2*log(e)*sqrt(mp(k)*ksp)/(pi**2+log(e)**2)

    !x force at right wall
    if (mask(ip+1,jp).eq.1) then
        dab=abs(x(ip+1)-xp(k))
        if (dab.lt.(dp/2+lam)) then
            fcolx=-ksp*dab-eta*up(k)
        end if
    end if

    !x force at left wall
    if (mask(ip,jp).eq.1) then
        dab=abs(x(ip)-xp(k))
        if (dab.lt.(dp/2+lam)) then
            fcolx=ksp*dab-eta*up(k)
        end if
    end if

    !y force at top wall
    if (mask(ip,jp+1).eq.1) then
        dab=abs(y(jp+1)-yp(k))
        if (dab.lt.(dp/2+lam)) then
            fcoly=-ksp*dab-eta*vp(k)
        end if
    end if

    !y force at bottom wall
    if (mask(ip,jp).eq.1) then
        dab=abs(y(jp)-yp(k))
        if (dab.lt.(dp/2+lam)) then
            fcoly=ksp*dab-eta*vp(k)
        end if
    end if


end subroutine lpt_collisions