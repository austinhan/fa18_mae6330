module lptsub
    use demoflow
    implicit none
    integer, parameter :: N=1
    real(WP), parameter :: dp=0.001_WP
    real(WP), parameter :: taup=dp**2.0_WP/knu/18.0_WP
    real(WP), dimension(N) :: xp,yp,up,vp  
end module lptsub

subroutine lpt_init
    use lptsub
    implicit none
    integer :: i
    ! Create initial positions of particles TBD
    do i=1,N
        xp(i)=0.5_WP
        yp(i)=0.6_WP
    end do
    ! Initial velocities
    up=1.0_WP
    vp=0.0_WP
end subroutine

subroutine lpt_solve
    use lptsub
    implicit none
    integer :: k
    ! real(WP) :: taup
    real(WP) :: fsn,Rep
    real(WP) :: ufp,vfp,xphalf,yphalf,dup,dvp,uphalf,vphalf,duphalf,dvphalf

    ! taup=dp**2.0_WP/knu/18.0_WP (if not defined in lptsub)
    do k=1,N
        ! ufp, vfp, fsn
        call lpt_fvel(xp(k),yp(k))
        xphalf=xp(k)+dt/2.0_WP*up(k)
        yphalf=yp(k)+dt/2.0_WP*vp(k)
        
        ! Periodic BC
        if (xphalf.ge.Lx) then
            xp(k)=xp(k)+dt*up(k)-Lx
            yp(k)=yp(k)+dt*vp(k)
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
        
        uphalf=up(k)+dt/2.0_WP*dup
        vphalf=vp(k)+dt/2.0_WP*dvp
        xp(k)=xp(k)+dt*uphalf
        yp(k)=yp(k)+dt*vphalf
        
        ! Periodic BC
        if (xp(k).ge.Lx) then
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

        up(k)=up(k)+dt*duphalf
        vp(k)=vp(k)+dt*dvphalf
        
        
    end do
contains
    subroutine lpt_fvel(xpo,ypo)
        implicit none
        integer :: i,j
        real(WP) :: xpo,ypo,stppt,magud
        stppt=-1
        ! Get surrounding indices
        i=1
        do while (stppt.lt.0)
        stppt=x(i)-xpo
        i=i+1
        end do

        stppt=-1
        j=1
        do while (stppt.lt.0)
        stppt=y(j)-ypo
        j=j+1
        end do
        
        ! interpolate velocities
        ! i=cell right of cell particle is indices
        ufp=1/(x(i)-x(i-1))/(y(j)-y(j-1))*( &
        &   u(i-1,j-1)*(x(i)-xpo)*(y(j)-ypo)+ &
        &   u(i,j-1)*(-x(i-1)+xpo)*(y(j)-ypo)+ &
        &   u(i-1,j)*(x(i)-xpo)*(-y(j-1)+ypo)+ &
        &   u(i,j)*(x(i-1)-xpo)*(y(j-1)-ypo))
        
        vfp=1/(x(i)-x(i-1))/(y(j)-y(j-1))*( &
        &   v(i-1,j-1)*(x(i)-xpo)*(y(j)-ypo)+ &
        &   v(i,j-1)*(-x(i-1)+xpo)*(y(j)-ypo)+ &
        &   v(i-1,j)*(x(i)-xpo)*(-y(j-1)+ypo)+ &
        &   v(i,j)*(x(i-1)-xpo)*(y(j-1)-ypo))

        magud=((ufp-up(k))**2+(vfp-vp(k))**2)**0.5
        Rep=dp*magud/knu
        fsn=1+0.15*Rep**0.687
        return
    end subroutine lpt_fvel
end subroutine lpt_solve

