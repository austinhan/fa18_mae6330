module pressure
  use demoflow
  implicit none
  
  ! Work arrays for pressure solver
  real(WP), dimension( 0:nx+1, 0:ny+1) :: Z,R,diag    ! Arrays for conjugate gradient solver
  real(WP), dimension( 1:nx  , 1:ny  ) :: Q           ! Arrays for conjugate gradient solver
  real(WP), dimension( 1:nx  , 1:ny  ) :: Pres        ! Pressure residual
  
end module pressure


! ======================================== !
! Pressure solver using conjugate gradient !
! preconditioned by an incomplete Cholesky !
! ======================================== !
subroutine pressure_solve
  use pressure
  implicit none
  integer  :: i,j
  real(WP) :: alpha,beta,rho1,rho2,rho0
  
  ! Prepare preconditioner - Gauss-Seidel at this point
  diag=1.0_WP
  do j=1,ny
     do i=1,nx
        diag(i,j)=sum(plap(i,j,:,0))
     end do
  end do
  
  ! This now makes our preconditioner an incomplete Cholesky factorization
  do j=2,ny
     do i=2,nx
        diag(i,j)=diag(i,j)&
             -plap(i,j,1,-1)*plap(i-1,j,1,+1)/diag(i-1,j)&
             -plap(i,j,2,-1)*plap(i,j-1,2,+1)/diag(i,j-1)
     end do
  end do
  
  ! Check initial pressure residual
  !P=0.0_WP ! allows us to change the initial guess
  rho0=0.0_WP
  do j=1,ny
     do i=1,nx
        Pres(i,j)=div(i,j)/dt-(sum(plap(i,j,1,-1:+1)*P(i-1:i+1,j))+sum(plap(i,j,2,-1:+1)*P(i,j-1:j+1)))
        rho0=rho0+Pres(i,j)*Pres(i,j)
     end do
  end do
  
  ! Preconditioner
  call precond(Pres,Z)
  rho1=0.0_WP
  do j=1,ny
     do i=1,nx
        rho1=rho1+Z(i,j)*Pres(i,j)
     end do
  end do
  rho2=0.0_WP
  
  ! Solve for pressure using conjugate gradient
  pit=0
  pressureloop: do while (pit.lt.maxpit .and. sqrt(abs(rho1)).gt.abscvg .and. sqrt(abs(rho1)/rho0).gt.relcvg)
     pit=pit+1
     if (pit.eq.1) then
        R=0.0_WP
        R(1:nx,1:ny)=Z(1:nx,1:ny)
     else
        R(1:nx,1:ny)=Z(1:nx,1:ny)+rho1/rho2*R(1:nx,1:ny)
     end if
     beta=0.0_WP
     do j=1,ny
        do i=1,nx
           Q(i,j)=sum(plap(i,j,1,-1:+1)*R(i-1:i+1,j))+sum(plap(i,j,2,-1:+1)*R(i,j-1:j+1))
           beta=beta+Q(i,j)*R(i,j)
        end do
     end do
     alpha=rho1/beta
     P=P+alpha*R
     Pres=Pres-alpha*Q
     call precond(Pres,Z)
     rho2=rho1
     rho1=0.0_WP
     do j=1,ny
        do i=1,nx
           rho1=rho1+Z(i,j)*Pres(i,j)
        end do
     end do
  end do pressureloop
  
contains
  
  ! ======================= !
  ! Pressure preconditioner !
  ! ======================= !
  subroutine precond(A,B)
    implicit none
    real(WP), dimension(1:nx  ,1:ny  ) :: A
    real(WP), dimension(0:nx+1,0:ny+1) :: B
    real(WP), dimension(0:nx+1,0:ny+1) :: tmp
    integer :: i,j
    ! Copy residual in larger array
    tmp(1:nx,1:ny)=A(1:nx,1:ny)
    ! Perform symmetric Gauss-Seidel
    do j=1,ny
       do i=1,nx
          B(i,j)=(tmp(i,j)-plap(i,j,1,-1)*B(i-1,j)-plap(i,j,2,-1)*B(i,j-1))/diag(i,j)
       end do
    end do
    do j=1,ny
       do i=1,nx
          tmp(i,j)=B(i,j)*diag(i,j)
       end do
    end do
    do j=ny,1,-1
       do i=nx,1,-1
          B(i,j)=(tmp(i,j)-plap(i,j,1,+1)*B(i+1,j)-plap(i,j,2,+1)*B(i,j+1))/diag(i,j)
       end do
    end do
    return
  end subroutine precond
  
end subroutine pressure_solve
