module levelset
use demoflow
implicit none
real(WP) :: cx,cy
real(WP) :: phihxp,phihyp,phihxm,phihym,Vhp,Vhm,Uhp,Uhm
integer :: iphi,jphi
real(WP), dimension(0:nx+1,0:ny+1) :: Hphi1,Hphi2

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

    do i=1,nx
        do j=1,ny

            Uhp=U(i+1,j)/2.0_WP+U(i,j)/2.0_WP
            Uhm=U(i-1,j)/2.0_WP+U(i,j)/2.0_WP
            Vhp=V(i,j+1)/2.0_WP+V(i,j)/2.0_WP
            Vhm=V(i,j-1)/2.0_WP+V(i,j)/2.0_WP

            if (Uhp.lt.0.0_WP) then
                phihxp=   phi(i,j)/3.0_WP+5.0_WP*phi(i+1,j)/6.0_WP-phi(i+2,j)/6.0_WP
            else
                phihxp=  -phi(i-1,j)/6.0_WP+5.0_WP*phi(i,j)/6.0_WP+2.0_WP*phi(i+1,j)/6.0_WP
            end if

            if (Uhm.lt.0.0_WP) then
                phihxm=  phi(i-1,j)/3.0_WP+5.0_WP*phi(i,j)/6.0_WP-phi(i+1,j)/6.0_WP
            else
                phihxm= -phi(i-2,j)/6.0_WP+5.0_WP*phi(i-1,j)/6.0_WP+phi(i,j)/3.0_WP
            end if

            if (Vhp.lt.0.0_WP) then
                phihyp=  phi(i,j)/3.0_WP+5.0_WP*phi(i,j+1)/6.0_WP-phi(i,j+2)/6.0_WP
            else
                phihyp= -phi(i,j-1)/6.0_WP+5.0_WP*phi(i,j)/6.0_WP+phi(i,j+1)/3.0_WP
            end if

            if (Vhm.lt.0.0_WP) then
                phihym=  phi(i,j-1)/3.0_WP+5.0_WP*phi(i,j)/6.0_WP-phi(i,j+1)/6.0_WP
            else
                phihym= -phi(i,j-2)/6.0_WP+5.0_WP*phi(i,j-1)/6.0_WP+phi(i,j)/3.0_WP
            end if

            Hphi1(i,j)=-dt/d*(phihxp*U(i,j)-phihxm*U(i,j)+phihyp*V(i,j)-phihym*V(i,j))
        enddo
    enddo
        !phi(1,:)=phi(2,:)
        !phi(nx,:)=phi(nx-1,:)
        !phi(:,1)=phi(:,2)
        !phi(:,ny)=phi(:,ny-1)

        phi(0,:)=phi(1,:)
        phi(nx+1,:)=phi(nx,:)
        phi(:,0)=phi(:,1)
        phi(:,ny+1)=phi(:,ny)   


phi = phi+Hphi1!+ABcoeff*(Hphi1-Hphi2)

        !phi(1,:)=phi(2,:)
        !phi(nx,:)=phi(nx-1,:)
        !phi(:,1)=phi(:,2)
        !phi(:,ny)=phi(:,ny-1)

        !phi(0,:)=phi(1,:)
        !phi(nx+1,:)=phi(nx,:)
        !phi(:,0)=phi(:,1)
        !phi(:,ny+1)=phi(:,ny)    

end subroutine