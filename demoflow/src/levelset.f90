module levelset
use demoflow
implicit none
real(WP) :: cx,cy
real(WP) :: phihxp,phihyp,phihxm,phihym,Vhp,Vhm,Uhp,Uhm
integer :: iphi,jphi,l
real(WP), dimension(0:nx+1,0:ny+1) :: Hphi1,Hphi2,diffH
real(WP), parameter :: eps=1.0e-9_WP
real(WP), dimension(-1:2) :: weno3_xp, weno3_yp
real(WP), dimension(-2:1) :: weno3_xm, weno3_ym

end module

subroutine levelsetinit
    use levelset

    ! Initialize residuals
    Hphi1=0.0_WP
    Hphi2=0.0_WP

    open(unit=88,file='lvlset.txt',action="write")
    do i=1,nx
        do j=1,ny

! Check negative and positive

        ! Nearest Point on Circle
        cx = 0.15*xm(i)/sqrt(xm(i)**2+(ym(j)-0.25)**2)
        cy = 0.25 + 0.15*(ym(j)-0.25)/sqrt(xm(i)**2+(ym(j)-0.25)**2)

        ! Outside circle and slot, but inside slot extended
        if ((xm(i).gt.-0.025_WP).and.(xm(i).lt.0.025_WP).and.(ym(j).lt.0.1)) then
            phi(i,j)= -min(sqrt((xm(i)+0.025)**2+(ym(j)-0.1)**2), sqrt((xm(i)-0.025)**2+(ym(j)-0.1)**2))

        ! Inside slot
        else if ((xm(i).gt.-0.025_WP).and.(xm(i).lt.0.025_WP).and.(ym(j).lt.0.35)) then
            phi(i,j)=-min(-xm(i)+0.025_WP,xm(i)+0.025_WP,0.35-ym(j))

        ! Inside circle, bottom left corner
        else if ((xm(i).gt.cx) .and. (xm(i).lt. -0.025_WP) .and. (ym(j).lt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), -xm(i)-0.025_WP)

        ! Inside circle, top left corner
        else if ((xm(i).gt.cx) .and. (xm(i).lt. -0.025_WP) .and. (ym(j).gt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), sqrt((-xm(i)-0.025_WP)**2+(ym(j)-0.35_WP)**2))

        ! Inside circle, bottom right corner
        else if ((xm(i).lt.cx) .and. (xm(i).gt. 0.025_WP) .and. (ym(j).lt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2),xm(i)-0.025_WP)

        ! Inside circle, top right corner
        else if ((xm(i).lt.cx) .and. (xm(i).gt. 0.025_WP) .and. (ym(j).gt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), sqrt((0.025_WP-xm(i))**2+(ym(j)-0.35_WP)**2))

        ! Inside circle, inside slot extended
        else if ((xm(i).gt.-0.025_WP).and.(xm(i).lt.0.025_WP).and.(ym(j).lt.cy).and.(ym(j).gt.0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), ym(j)-0.35_WP)

        ! Outside Circle and slot extended
        else if ((sqrt(xm(i)**2+(ym(j)-0.25)**2)).gt.0.15) then
            phi(i,j) = -sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2)
        end if

        ! Edges
        phi(0,:)=phi(1,:)
        phi(nx+1,:)=phi(nx,:)
        phi(:,0)=phi(:,1)
        phi(:,ny+1)=phi(:,ny)
        !phinew=phi
        end do
    end do
write(88,*) phi
end subroutine

subroutine levelset_step
    use levelset

    ! Remember previous residual
    Hphi2=Hphi1; Hphi1=0.0_WP
    if (time-dt.le.0.0_WP) Hphi2=0.0_WP
    do i=2,nx-1
        do j=2,ny-1

            if (U(i+1,j).lt.0.0_WP) then
                phihxp=   phi(i,j)/3.0_WP+5.0_WP*phi(i+1,j)/6.0_WP-phi(i+2,j)/6.0_WP
            else
                phihxp=  -phi(i-1,j)/6.0_WP+5.0_WP*phi(i,j)/6.0_WP+2.0_WP*phi(i+1,j)/6.0_WP
            end if

            if (U(i,j).lt.0.0_WP) then
                phihxm=  phi(i-1,j)/3.0_WP+5.0_WP*phi(i,j)/6.0_WP-phi(i+1,j)/6.0_WP
            else
                phihxm= -phi(i-2,j)/6.0_WP+5.0_WP*phi(i-1,j)/6.0_WP+phi(i,j)/3.0_WP
            end if

            if (V(i,j+1).lt.0.0_WP) then
                phihyp=  phi(i,j)/3.0_WP+5.0_WP*phi(i,j+1)/6.0_WP-phi(i,j+2)/6.0_WP
            else
                phihyp= -phi(i,j-1)/6.0_WP+5.0_WP*phi(i,j)/6.0_WP+phi(i,j+1)/3.0_WP
            end if

            if (V(i,j).lt.0.0_WP) then
                phihym=  phi(i,j-1)/3.0_WP+5.0_WP*phi(i,j)/6.0_WP-phi(i,j+1)/6.0_WP
            else
                phihym= -phi(i,j-2)/6.0_WP+5.0_WP*phi(i,j-1)/6.0_WP+phi(i,j)/3.0_WP
            end if

            Hphi1(i,j)=-dt/d*(phihxp*U(i+1,j)-phihxm*U(i,j)+phihyp*V(i,j+1)-phihym*V(i,j))
        enddo
    enddo

        diffH=3*Hphi1-Hphi2

phi = phi+ABcoeff*diffH

        !if (reinitcount.eq.10) call levelset_reinit

        phi(1,:)=phi(2,:)
        phi(nx,:)=phi(nx-1,:)
        phi(:,1)=phi(:,2)
        phi(:,ny)=phi(:,ny-1)

        phi(0,:)=phi(1,:)
        phi(nx+1,:)=phi(nx,:)
        phi(:,0)=phi(:,1)
        phi(:,ny+1)=phi(:,ny)    

contains

subroutine levelset_reinit
    use levelset
    implicit none
    real(WP) ::dtau, G, dphidx_m,dphidx_p,dphidy_m,dphidy_p
    real(WP), dimension(0:nx+1,0:ny+1) :: Se,phi0

    dtau=d/30.0_WP
    phi0=phi
    Se=phi/sqrt(phi**2+d**2)
    do l=1,10
        do i=1,nx
            do j=1,ny

              ! Calculate WENO3 coefficients

              call weno3_coeff

              ! Use WENO3 to calculate grad(G)

              dphidx_m=sum(weno3_xm(-2:+1)*phi(i-2:i+1,j))

              dphidx_p=sum(weno3_xp(-1:+2)*phi(i-1:i+2,j))

              dphidy_m=sum(weno3_ym(-2:+1)*phi(i,j-2:j+1))

              dphidy_p=sum(weno3_yp(-1:+2)*phi(i,j-1:j+2))

              if (phi0(i,j).gt.0) then
                G=sqrt(max(max(dphidx_m,0.0_WP)**2,min(dphidx_p,0.0_WP)**2) &
                &     +max(max(dphidy_m,0.0_WP)**2,min(dphidy_p,0.0_WP)**2))-1.0_WP
              else if (phi0(i,j).lt.0) then
                G=sqrt(max(min(dphidx_m,0.0_WP)**2,max(dphidx_p,0.0_WP)**2) &
                &     +max(min(dphidy_m,0.0_WP)**2,max(dphidy_p,0.0_WP)**2))-1.0_WP
              else
                G=0.0_WP
            
            end if

            !print *,dphidx_m,G

            phi(i,j)=phi(i,j)-Se(i,j)*G*dtau

            end do
        end do
    end do



end subroutine levelset_reinit


subroutine weno3_coeff

    implicit none

    real(WP) :: r_minus,r_plus,w_minus,w_plus


    ! X DIRECTION ============================================

    ! Direction x - left biased stencil

    r_minus=(eps+(phi(i  ,j)-2.0_WP*phi(i-1,j)+phi(i-2,j))**2.0_WP)&

         & /(eps+(phi(i+1,j)-2.0_WP*phi(i  ,j)+phi(i-1,j))**2.0_WP)

    w_minus=1.0_WP/(1.0_WP+2.0_WP*r_minus**2.0_WP)

    weno3_xm(-2)=0.5_WP/d*(       +1.0_WP*w_minus)

    weno3_xm(-1)=0.5_WP/d*(-1.0_WP-3.0_WP*w_minus)

    weno3_xm( 0)=0.5_WP/d*(       +3.0_WP*w_minus)

    weno3_xm(+1)=0.5_WP/d*(+1.0_WP-1.0_WP*w_minus)

    ! Direction x - right biased stencil

    r_plus=(eps+(phi(i+2,j)-2.0_WP*phi(i+1,j)+phi(i  ,j))**2.0_WP)&

         &/(eps+(phi(i+1,j)-2.0_WP*phi(i  ,j)+phi(i-1,j))**2.0_WP)

    w_plus=1.0_WP/(1.0_WP+2.0_WP*r_plus**2.0_WP)

    weno3_xp(-1)=0.5_WP/d*(-1.0_WP+1.0_WP*w_plus)

    weno3_xp( 0)=0.5_WP/d*(       -3.0_WP*w_plus)

    weno3_xp(+1)=0.5_WP/d*(+1.0_WP+3.0_WP*w_plus)

    weno3_xp(+2)=0.5_WP/d*(       -1.0_WP*w_plus)

    

    ! Y DIRECTION ============================================

    ! Direction y - left biased stencil

    r_minus=(eps+(phi(i,j  )-2.0_WP*phi(i,j-1)+phi(i,j-2))**2.0_WP)&

         & /(eps+(phi(i,j+1)-2.0_WP*phi(i,j  )+phi(i,j-1))**2.0_WP)

    w_minus=1.0_WP/(1.0_WP+2.0_WP*r_minus**2.0_WP)

    weno3_ym(-2)=0.5_WP/d*(       +1.0_WP*w_minus)

    weno3_ym(-1)=0.5_WP/d*(-1.0_WP-3.0_WP*w_minus)

    weno3_ym( 0)=0.5_WP/d*(       +3.0_WP*w_minus)

    weno3_ym(+1)=0.5_WP/d*(+1.0_WP-1.0_WP*w_minus)

    ! Direction y - right biased stencil

    r_plus=(eps+(phi(i,j+2)-2.0_WP*phi(i,j+1)+phi(i,j  ))**2.0_WP)&

         &/(eps+(phi(i,j+1)-2.0_WP*phi(i,j  )+phi(i,j-1))**2.0_WP)

    w_plus=1.0_WP/(1.0_WP+2.0_WP*r_plus**2.0_WP)

    weno3_yp(-1)=0.5_WP/d*(-1.0_WP+1.0_WP*w_plus)

    weno3_yp( 0)=0.5_WP/d*(       -3.0_WP*w_plus)

    weno3_yp(+1)=0.5_WP/d*(+1.0_WP+3.0_WP*w_plus)

    weno3_yp(+2)=0.5_WP/d*(       -1.0_WP*w_plus)

    return

  end subroutine weno3_coeff

end subroutine