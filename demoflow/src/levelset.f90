module levelset
use demoflow
implicit none
real(WP) :: cx,cy
real(WP) :: phihxp,phihyp,phihxm,phihym,Vhp,Vhm,Uhp,Uhm
integer :: iphi,jphi
real(WP), dimension(0:nx+1,0:ny+1) :: phinew

end module

subroutine levelsetinit
    use levelset
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
        phinew=phi
        end do
    end do
write(88,*) phi
end subroutine

subroutine levelset_step
    use levelset

    do i=1,nx
        do j=1,ny
            

            if (U(i,j).lt.0.0_WP) then
                phihxp= -phi(i,j)/6+5/6*phi(i+1,j)+phi(i+2,j)/3
                phihxm= -phi(i-1,j)/6+5/6*phi(i,j)+phi(i+1,j)/3
            else
                phihxp=  phi(i-1,j)/3+5*phi(i,j)/6-phi(i+1,j)/6
                phihxm=  phi(i-2,j)/3+5*phi(i-1,j)/6-phi(i,j)/6
            end if

            if (V(i,j).lt.0.0_WP) then
                phihyp=-phi(i,j)/6+5/6*phi(i,j+1)+phi(i,j+2)/3
                phihym=-phi(i,j-1)/6+5/6*phi(i,j)+phi(i,j+1)/3
            else
                phihyp=phi(i,j-1)/3+5*phi(i,j)/6-phi(i,j+1)/6
                phihym=phi(i,j-2)/3+5*phi(i,j-1)/6-phi(i,j)/6
            end if

            phinew(i,j)=phi(i,j)-dt/d* &
      &     ((phi(iphi,j   )*U(iphi,j   )-phi(iphi-1,j     )*U(iphi-1,j     ))+ &
      &      (phi(i   ,jphi)*V(i   ,jphi)-phi(i     ,jphi-1)*V(i     ,jphi-1)))
        enddo
    enddo
phi=phinew








end subroutine