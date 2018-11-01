module levelset
use demoflow
implicit none
real(WP), dimension(0:nx+1,0:ny+1) :: phi=0
real(WP) :: cx,cy




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
        
        ! Inside slot, not done
        else if ((xm(i).gt.-0.025_WP).and.(xm(i).lt.0.025_WP).and.(ym(j).lt.0.1)) then
            phi(i,j)=0

        ! Inside circle, bottom left corner
        else if ((xm(i).gt.cx) .and. (xm(i).lt. -0.025_WP) .and. (ym(j).lt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), -xm(i)+0.025_WP)

        ! Inside circle, top left corner
        else if ((xm(i).gt.cx) .and. (xm(i).lt. -0.025_WP) .and. (ym(j).gt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), sqrt((-xm(i)+0.025_WP)**2+(ym(j)-0.35_WP)**2))

        ! Inside circle, bottom right corner
        else if ((xm(i).lt.cx) .and. (xm(i).gt. 0.025_WP) .and. (ym(j).lt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2),xm(i)-0.025_WP)

        ! Inside circle, top right corner
        else if ((xm(i).lt.cx) .and. (xm(i).gt. 0.025_WP) .and. (ym(j).gt. 0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), sqrt((0.025_WP-xm(i))**2+(ym(j)-0.35_WP)**2))

        ! Inside circle, inside slot extended
        else if ((xm(i).gt.-0.025_WP).and.(xm(i).lt.0.025_WP).and.(ym(j).lt.cy).and.(ym(j).gt.0.35_WP)) then
            phi(i,j) = min(sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2), ym(j)-0.035_WP)

            
        ! Outside Circle and slot extended
        else
            phi(i,j) = -sqrt((xm(i)-cx)**2 + (ym(j)-cy)**2)
        end if

        ! Edges
        phi(0,:)=phi(1,:)
        phi(nx+1,:)=phi(nx,:)
        phi(:,0)=phi(:,1)
        phi(:,ny+1)=phi(:,ny)

        end do
    end do
write(88,*) phi


end subroutine