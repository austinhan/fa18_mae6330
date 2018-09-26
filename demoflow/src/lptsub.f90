module lptsub
    use demoflow
    implicit none
    real(WP), parameter :: dp=0.01_WP
    real(WP), parameter :: taup=dp**2.0_WP/knu/18.0_WP
end module lptsub

subroutine lpt_init
    use lptsub
    implicit none
    integer :: i
    ! Create initial positions/velocites of particles TBD
    do i=1,Np
        xp(i)=.5_WP/i
        yp(i)=0.3_WP/i
    end do

    up=1.0_WP
    vp=0.0_WP
    liter= int(dt/taup)

end subroutine

subroutine lpt_solve
    use lptsub
    implicit none
    integer :: k
    ! real(WP) :: taup
    real(WP) :: fsn,Rep
    real(WP) :: ufp,vfp,xphalf,yphalf,dup,dvp,uphalf,vphalf,duphalf,dvphalf

    liter= int(dt/taup)

    ! taup=dp**2.0_WP/knu/18.0_WP (if not defined in lptsub)
    do k=1,Np
        ! ufp, vfp, fsn
        call lpt_fvel(xp(k),yp(k))
        xphalf=xp(k)+taup/2.0_WP*up(k)
        yphalf=yp(k)+taup/2.0_WP*vp(k)
        ! Periodic BC
        if (xphalf.ge.Lx-Lx/nx) then
            xp(k)=xp(k)+taup*up(k)-Lx
            yp(k)=yp(k)+taup*vp(k)
            ! call lpt_fvel(xp(k),yp(k))
            ! up(k)=ufp
            ! vp(k)=vfp
            up(k)=1.0_WP
            vp(k)=0.0_WP
            cycle
        end if

        ! Calculate dup/dt
        dup=fsn*(ufp-up(k))/taup+gravity(1)
        dvp=fsn*(vfp-vp(k))/taup+gravity(2)
        
        uphalf=up(k)+taup/2.0_WP*dup
        vphalf=vp(k)+taup/2.0_WP*dvp
        xp(k)=xp(k)+taup*uphalf
        yp(k)=yp(k)+taup*vphalf
        
        ! Periodic BC
        if (xp(k).ge.Lx-Lx/nx) then
            xp(k)=xp(k)-Lx
            ! call lpt_fvel(xp(k),yp(k))
            ! up(k)=ufp
            ! vp(k)=vfp
            up(k)=1.0_WP
            vp(k)=0.0_WP
            cycle
        end if

        ! new ufp,vfp -> magud
        ! new Rep -> fsn
        call lpt_fvel(xphalf,yphalf)
        ! Calculate dup/dt half
        duphalf=fsn*(ufp-uphalf)/taup+gravity(1)
        dvphalf=fsn*(vfp-vphalf)/taup+gravity(2)

        up(k)=up(k)+taup*duphalf
        vp(k)=vp(k)+taup*dvphalf
        
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
        return
    end subroutine lpt_fvel
end subroutine lpt_solve

