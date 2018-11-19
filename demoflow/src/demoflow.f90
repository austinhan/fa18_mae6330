! ======================================================================== !
! Code for solving a 2D incompressible flow using a fractional step method !
! Written by O. Desjardins in preparation for GIAN course, IIT Kanpur 2017 !
! ======================================================================== !
module demoflow
  implicit none
  integer, parameter :: WP=8                          ! This is our working precision
  integer, parameter :: SP=4                          ! This is our visualization precision
  
  ! ==========================================
  ! ========= PARAMETERS TO MODIFY ===========
  ! ==========================================
  ! Time integration
  real(WP), parameter :: maxdt=1e-4_WP
  real(WP), parameter :: maxCFL=0.5_WP
  real(WP), parameter :: viztime=0.1_WP
  ! End of time integration
  real(WP), parameter :: maxtime=1.0_WP
  integer , parameter :: maxstep=10000
  ! Pressure convergence criterion
  real(WP), parameter :: relcvg=1.0e-4_WP
  real(WP), parameter :: abscvg=1.0e-4_WP
  integer , parameter :: maxpit=100
  ! Mesh size
  integer,  parameter :: nx=100
  integer,  parameter :: ny=100
  ! Domain size
  real(WP), parameter :: Lx=1.0_WP
  real(WP), parameter :: Ly=1.0_WP
  ! Kinematic viscosity
  real(WP), parameter :: knu=0.0001_WP
  ! Gravity
  real(WP), parameter, dimension(2) :: gravity=(/0.0_WP,0.0_WP/)
  ! Do we use levelset
  logical , parameter :: use_levelset=.true.
  real(WP), dimension(1:nx,1:ny) :: jcx, jcy
  real(WP), dimension(-1:nx+2,-1:ny+2) :: curv
  real(WP), parameter :: sigma=1.0_WP
    ! Levelset field
  real(WP), dimension(-1:nx+2,-1:ny+2) :: G
  real(WP) :: add, icurv, rhoij, nrhoij, rhog, rhol, rhori,rhole,rhoab,rhobe, thet
  ! ==========================================
  
  ! Mesh
  real(WP) :: d                                       ! Mesh size - cells are assumed square
  real(WP), dimension(0:nx+2) :: x                    ! Location of cell faces in x
  real(WP), dimension(0:ny+2) :: y                    ! Location of cell faces in y
  real(WP), dimension(0:nx+1) :: xm                   ! Location of cell center in x
  real(WP), dimension(0:ny+1) :: ym                   ! Location of cell center in y
  integer , dimension(0:nx+1,0:ny+1) :: mask          ! Masks are used to specify walls
  
  ! Solution vectors
  real(WP), dimension(0:nx+1,0:ny+1) :: U,V,P         ! Velocity and pressure variables
  
  ! Residual vectors
  real(WP), dimension(0:nx+1,0:ny+1) :: HU1,HU2       ! U velocity residuals
  real(WP), dimension(0:nx+1,0:ny+1) :: HV1,HV2       ! V velocity residuals
  real(WP), dimension(1:nx  ,1:ny  ) :: div           ! Velocity divergence
  
  ! Time info
  real(WP) :: time,CFL,dt,dt_old
  integer  :: ntime,pit
  
  ! Discretization coefficients
  real(WP), dimension(1:nx,1:ny,1:2,-1:+1) :: plap   ! Pressure Laplacian
  real(WP), dimension(2:nx,1:ny,1:2,-1:+1) :: ulap   ! U velocity Laplacian
  real(WP), dimension(1:nx,2:ny,1:2,-1:+1) :: vlap   ! V velocity Laplacian
  real(WP) :: ABcoeff                                ! Adams-Bashforth coefficient
  
  ! Named constant
  real(WP), parameter :: Pi=3.141592653589793_WP     ! Pi
  
end module demoflow


! ================ !
! Main driver code !
! ================ !
program main
  use demoflow
  implicit none
  real(WP) :: walltime_ref,walltime
  
  ! Initialize timer
  call cpu_time(walltime_ref)
  
  ! Initialize mesh, walls, and boundary conditions
  call demoflow_setup
  
  ! Initialize various aspects of the code
  call demoflow_init

  ! Initialize levelset method
  if (use_levelset) call levelset_init
  
  ! Initialize time
  time=0.0_WP
  dt=maxdt
  ntime=0
  
  ! Initialize visualization
  call visualize_init
  
  ! Main time loop
  timeloop: do while (time.lt.maxtime .and. ntime.lt.maxstep)
     
     ! Adjust timestep size
     call time_adjust
     
     ! Increment time
     time=time+dt; ntime=ntime+1
     
     ! Some output to the screen
     if (ntime.eq.1) write(*,'(a12,a2,6a12)') 'Step','  ','Time  ','CFLmax','Umax  ','Vmax  ','Divergence','Piterations'
     write(*,'(i12,a2,1ES12.5,1F12.4,3ES12.3,i12)') ntime,'  ',time,CFL,maxval(abs(U)),maxval(abs(V)),maxval(abs(div)),pit
     
     ! Levelset
     if (use_levelset) call levelset_step
     
     ! Velocity step
     call velocity_step
     
     ! Pressure step
      call pressure_step
     
     ! Dump data for visualization
     call visualize_dump
     
  end do timeloop
  
  ! Get final time
  call cpu_time(walltime)
  print*,'Time taken: ',walltime-walltime_ref
  
end program main


! ==========================================
! ========= PARAMETERS TO MODIFY ===========
! ==========================================
subroutine demoflow_setup
  use demoflow
  implicit none
  integer  :: i,j
  
  ! Initialize mesh
  do i=1,nx+1
     x(i)=real(i-1,WP)*Lx/real(nx,WP)-0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j)=real(j-1,WP)*Ly/real(ny,WP)-0.5_WP*Ly
  end do

  ! Create position of cell centers
  do i=1,nx
     xm(i)=0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)=0.5_WP*(y(j)+y(j+1))
  end do
  
  ! Mask out walls
  mask=0
  
  ! Initial conditions - solid body rotation
  do j=0,ny+1
     do i=0,nx+1
        U(i,j)=-2.0_WP*Pi*ym(j)
        V(i,j)=+2.0_WP*Pi*xm(i)
     end do
  end do
  U=0.0_WP
  V=0.0_WP
  P=0.0_WP
  
  return
end subroutine demoflow_setup
! ==========================================


! ========================= !
! Basic code initialization !
! - output directories      !
! - mesh                    !
! - FV coefficients         !
! - walls                   !
! - Laplacian operators     !
! ========================= !
subroutine demoflow_init
  use demoflow
  implicit none
  integer :: i,j
  integer, dimension(:), allocatable :: seed
  
  ! Seed the random number generator
  call RANDOM_SEED(size=i)
  allocate(seed(i))
  call system_clock(count=seed(1))
  seed(1:i)=seed(1)
  call RANDOM_SEED(put=seed)
  deallocate(seed)
  
  ! Apply Neumann condition on masks
  mask(   0,:)=mask( 1,:)
  mask(nx+1,:)=mask(nx,:)
  mask(:,   0)=mask(:, 1)
  mask(:,ny+1)=mask(:,ny)
  
  ! Store mesh size
  d=Lx/real(nx,WP)
  
  ! Check that cells are square
  if (abs(d-Ly/real(ny,WP)).gt.1.0e-15_WP) stop ('Improper mesh definition')

  ! Extend the mesh to provide ghost cells
  x(0)=x(1)-d; x(nx+2)=x(nx+1)+d; xm(0)=xm(1)-d; xm(nx+1)=xm(nx)+d
  y(0)=y(1)-d; y(ny+2)=y(ny+1)+d; ym(0)=ym(1)-d; ym(ny+1)=ym(ny)+d
  

  
  ! Create U Laplacian for viscous term
  ulap(:,:,1,-1)=+1.0_WP/d**2
  ulap(:,:,1, 0)=-2.0_WP/d**2
  ulap(:,:,1,+1)=+1.0_WP/d**2
  ulap(:,:,2,-1)=+1.0_WP/d**2
  ulap(:,:,2, 0)=-2.0_WP/d**2
  ulap(:,:,2,+1)=+1.0_WP/d**2
  
  ! Neumann conditions if non-periodic
  ulap(:,1,2, 0)=ulap(:,1,2,0)+ulap(:,1,2,-1)
  ulap(:,1,2,-1)=0.0_WP
  ulap(:,ny,2, 0)=ulap(:,ny,2,0)+ulap(:,ny,2,+1)
  ulap(:,ny,2,+1)=0.0_WP
  
  ! Handle walls in the domain
  do j=1,ny
     do i=2,nx
        ! If wall cell, zero out Laplacian
        if (mask(i,j).eq.1.or.mask(i-1,j).eq.1) then
           ulap(i,j,:,:)=0.0_WP
        else
           ! If neighbor to wall cell in y, change stencil
           if (mask(i,j+1).eq.1.or.mask(i-1,j+1).eq.1) then
              ulap(i,j,2, 0)=ulap(i,j,2,0)-ulap(i,j,2,+1)
              ulap(i,j,2,+1)=0.0_WP
           end if
           if (mask(i,j-1).eq.1.or.mask(i-1,j-1).eq.1) then
              ulap(i,j,2, 0)=ulap(i,j,2,0)-ulap(i,j,2,-1)
              ulap(i,j,2,-1)=0.0_WP
           end if
        end if
     end do
  end do
  
  ! Create V Laplacian for viscous term
  vlap(:,:,1,-1)=+1.0_WP/d**2
  vlap(:,:,1, 0)=-2.0_WP/d**2
  vlap(:,:,1,+1)=+1.0_WP/d**2
  vlap(:,:,2,-1)=+1.0_WP/d**2
  vlap(:,:,2, 0)=-2.0_WP/d**2
  vlap(:,:,2,+1)=+1.0_WP/d**2
  
  ! Handle walls in the domain
  do j=2,ny
     do i=1,nx
        ! If wall cell, zero out Laplacian
        if (mask(i,j).eq.1.or.mask(i,j-1).eq.1) then
           vlap(i,j,:,:)=0.0_WP
        else
           ! If neighbor to wall cell in x, change stencil
           if (mask(i+1,j).eq.1.or.mask(i+1,j-1).eq.1) then
              vlap(i,j,1, 0)=vlap(i,j,1,0)-vlap(i,j,1,+1)
              vlap(i,j,1,+1)=0.0_WP
           end if
           if (mask(i-1,j).eq.1.or.mask(i-1,j-1).eq.1) then
              vlap(i,j,1, 0)=vlap(i,j,1,0)-vlap(i,j,1,-1)
              vlap(i,j,1,-1)=0.0_WP
           end if
        end if
     end do
  end do
  
  ! Initialize residuals
  HU1=0.0_WP; HU2=0.0_WP
  HV1=0.0_WP; HV2=0.0_WP
  
  ! Clean up initial conditions
  do j=1,ny
     do i=1,nx
        if (maxval(mask(i-1:i,j)).eq.1) U(i,j)=0.0_WP
        if (maxval(mask(i,j-1:j)).eq.1) V(i,j)=0.0_WP
     end do
  end do
  
  return
end subroutine demoflow_init


! ============================ !
! Solve Navier-Stokes equation !
! ============================ !
subroutine velocity_step
  use demoflow
  implicit none
  integer  :: i,j
  real(WP) :: visc,conv
  
  ! Remember previous residual
  HU2=HU1; HU1=0.0_WP
  HV2=HV1; HV1=0.0_WP
  
  ! Velocity update
  conv=0.25_WP*dt/d
  visc=dt*knu
  do j=1,ny
     do i=2,nx
        ! Compute U velocity residual if no wall
        if (maxval(mask(i-1:i,j)).eq.0) then
           HU1(i,j)=-conv*((U(i,j  )+U(i+1,j))*(U(i,j  )+U(i+1,j  ))-(U(i-1,j)+U(i,j  ))*(U(i-1,j)+U(i  ,j))) &
                &   -conv*((U(i,j+1)+U(i  ,j))*(V(i,j+1)+V(i-1,j+1))-(U(i  ,j)+U(i,j-1))*(V(i  ,j)+V(i-1,j))) &
                &   +visc*(sum(ulap(i,j,1,-1:+1)*U(i-1:i+1,j))+sum(ulap(i,j,2,-1:+1)*U(i,j-1:j+1)))
        end if
     end do
  end do
  do j=2,ny
     do i=1,nx
        ! Compute V velocity residual if no wall
        if (maxval(mask(i,j-1:j)).eq.0) then
           HV1(i,j)=-conv*((V(i+1,j)+V(i,j  ))*(U(i+1,j)+U(i+1,j-1))-(V(i,j  )+V(i-1,j))*(U(i,j  )+U(i,j-1))) &
                &   -conv*((V(i  ,j)+V(i,j+1))*(V(i  ,j)+V(i  ,j+1))-(V(i,j-1)+V(i  ,j))*(V(i,j-1)+V(i,j  ))) &
                &   +visc*(sum(vlap(i,j,1,-1:+1)*V(i-1:i+1,j))+sum(vlap(i,j,2,-1:+1)*V(i,j-1:j+1)))
        end if
     end do
  end do
  
  ! Adams-Bashforth
  U=U+HU1+ABcoeff*(HU1-HU2)
  V=V+HV1+ABcoeff*(HV1-HV2)
  
  ! Apply Dirichlet in x- (nothing to do)
  ! Apply Neumann in x+ with negative velocity clipping
  do j=1,ny
     if (mask(nx,j).eq.0) then
        U(nx+1,j)=max(U(nx,j),0.0_WP)
     end if
     if (maxval(mask(nx,j-1:j)).eq.0) then
        V(nx+1,j)=V(nx,j)
     end if
  end do
  ! Apply Neumann in y+ and y- for U
  do i=1,nx+1
     if (maxval(mask(i-1:i,1)).eq.0) then
        U(i,0)=U(i,1)
     end if
     if (maxval(mask(i-1:i,ny)).eq.0) then
        U(i,ny+1)=U(i,ny)
     end if
  end do
  ! Apply no-penetration in y+ and y- for V
  do i=1,nx
     if (mask(i,1).eq.0) then
        V(i,1)=0.0_WP
     end if
     if (mask(i,ny).eq.0) then
        V(i,ny+1)=0.0_WP
     end if
  end do
  
  return
end subroutine velocity_step


! =============================== !
! Solve pressure Poisson equation !
! Uses a conjugate gradient with  !     
! a Gauss-Seidel preconditioner   !
! =============================== !
subroutine pressure_step
  use demoflow
  implicit none
  integer  :: i,j
  real(WP) :: left_mfr,right_mfr
  real(WP) :: left_area,right_area
  
  ! Enforce global mass conservation
  ! Left flow rate
  left_mfr=0.0_WP
  left_area=0.0_WP
  do j=1,ny
     if (mask(1,j).eq.0) then
        left_mfr=left_mfr+d*U(1,j)
        left_area=left_area+d
     end if
  end do
  if (left_area.gt.0.0_WP) then
     left_mfr=left_mfr/left_area
  else
     left_mfr=0.0_WP
  end if
  ! Right flow rate
  right_mfr=0.0_WP
  right_area=0.0_WP
  do j=1,ny
     if (mask(nx,j).eq.0) then
        right_mfr=right_mfr+d*U(nx+1,j)
        right_area=right_area+d
     end if
  end do
  if (right_area.gt.0.0_WP) then
     right_mfr=right_mfr/right_area
  else
     right_mfr=0.0_WP
  end if
  ! Correct right flow rate
  if (right_area.eq.0.0_WP.and.left_mfr.gt.0.0_WP) &
       stop ("Error: Cannot handle an inflow without outflow!")
  do j=1,ny
     if (mask(nx,j).eq.0) then
        U(nx+1,j)=U(nx+1,j)+(left_mfr*left_area/right_area-right_mfr)
     end if
  end do
  
  ! Compute divergence
  div=0.0_WP
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) div(i,j)=(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j))/d
     end do
  end do

  ! Pressure laplcian created in here
 
  call levelset_jump

  
  ! Neumann conditions if non-periodic
  plap(1,:,1, 0)=plap(1,:,1,0)+plap(1,:,1,-1)
  plap(1,:,1,-1)=0.0_WP
  plap(nx,:,1, 0)=plap(nx,:,1,0)+plap(nx,:,1,+1)
  plap(nx,:,1,+1)=0.0_WP
  plap(:,1,2, 0)=plap(:,1,2,0)+plap(:,1,2,-1)
  plap(:,1,2,-1)=0.0_WP
  plap(:,ny,2, 0)=plap(:,ny,2,0)+plap(:,ny,2,+1)
  plap(:,ny,2,+1)=0.0_WP
  
  ! Handle walls in the domain
  do j=1,ny
     do i=1,nx
        ! If wall cell, change Laplacian to identity
        if (mask(i,j).eq.1) then
           plap(i,j,:,:)=0.0_WP
           plap(i,j,1,0)=1.0_WP
        else
           ! If neighbor to wall cell, use Neumann
           if (mask(i+1,j).eq.1) then
              plap(i,j,1, 0)=plap(i,j,1,0)+plap(i,j,1,+1)
              plap(i,j,1,+1)=0.0_WP
           end if
           if (mask(i-1,j).eq.1) then
              plap(i,j,1, 0)=plap(i,j,1,0)+plap(i,j,1,-1)
              plap(i,j,1,-1)=0.0_WP
           end if
           if (mask(i,j+1).eq.1) then
              plap(i,j,2, 0)=plap(i,j,2,0)+plap(i,j,2,+1)
              plap(i,j,2,+1)=0.0_WP
           end if
           if (mask(i,j-1).eq.1) then
              plap(i,j,2, 0)=plap(i,j,2,0)+plap(i,j,2,-1)
              plap(i,j,2,-1)=0.0_WP
           end if
        end if
     end do
  end do

  ! Solve the pressure Poisson equation
  call pressure_solve
  
 !call sleep(1)
  
  ! Correct U velocity
  do j=1,ny
     do i=2,nx
        if (maxval(mask(i-1:i,j)).eq.0) then
           U(i,j)=U(i,j)-dt*(P(i,j)-P(i-1,j))/d
           if (G(i-1,j  )*G(i,j).lt.0) then
            icurv=curv(i-1,j)*(1+G(i-1,j)/(G(i,j)-G(i-1,j)))+curv(i,j)*(-G(i-1,j)/(G(i,j)-G(i-1,j)))
            !icurv=1.0_WP/0.15_WP
            add=sigma*icurv
            if (G(i-1,j).gt.G(i,j)) add=-add
            U(i,j)=U(i,j)+dt*add/d
           end if
        end if
     end do
  end do
  
  ! Correct V velocity
  do j=2,ny
     do i=1,nx
        if (maxval(mask(i,j-1:j)).eq.0) then
           V(i,j)=V(i,j)-dt*(P(i,j)-P(i,j-1))/d
           if (G(i,j-1)*G(i,j).lt.0) then
            icurv=curv(i,j-1)*(1+G(i,j-1)/(G(i,j)-G(i,j-1)))+curv(i,j)*(-G(i,j-1)/(G(i,j)-G(i,j-1)))
            !print *,icurv
            !icurv=1.0_WP/.15_WP
            add=sigma*icurv
            if (G(i,j-1).gt.G(i,j)) add=-add
            V(i,j)=V(i,j)+dt*add/d
           end if
        end if
     end do
  end do
  
  ! Recalculate divergence
  div=0.0_WP
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) div(i,j)=(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j))/d
     end do
  end do
  
end subroutine pressure_step


! ========================== !
! Time step size calculation !
! ========================== !
subroutine time_adjust
  use demoflow
  implicit none
  integer  :: i,j
  real(WP), parameter :: alpha=0.7_WP
  
  ! Get max CFL
  CFL=0.0_WP
  do j=1,ny
     do i=1,nx
        CFL=max(CFL,abs(U(i,j))/d)
        CFL=max(CFL,abs(V(i,j))/d)
        CFL=max(CFL,4.0_WP*knu/d**2)
     end do
  end do
  CFL=CFL*dt
  
  ! Adjust time step size
  dt_old=dt
  dt=min(maxCFL/(CFL+epsilon(1.0_WP))*dt_old,maxdt)
  if (dt.gt.dt_old) dt=alpha*dt+(1.0_WP-alpha)*dt_old
  
  ! Adams-Bashforth coefficient
  ABcoeff=0.5_WP*dt/dt_old
  if (ntime.eq.1) ABcoeff=0.0_WP
  
  return
end subroutine time_adjust


