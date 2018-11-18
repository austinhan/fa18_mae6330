module levelset
  use demoflow
  implicit none
  

   
  ! Normal vector
  real(WP), dimension(-1:nx+2,-1:ny+2) :: normx
  real(WP), dimension(-1:nx+2,-1:ny+2) :: normy

  ! Residual vector
  real(WP), dimension(-1:nx+2,-1:ny+2) :: HG1,HG2

  ! QUICK coefficients
  real(WP), dimension(-1:+1), parameter :: qp=(/-1.0_WP/6.0_WP,+5.0_WP/6.0_WP,+2.0_WP/6.0_WP/)
  real(WP), dimension( 0:+2), parameter :: qm=(/+2.0_WP/6.0_WP,+5.0_WP/6.0_WP,-1.0_WP/6.0_WP/)

  ! Reinitialization step
  integer , parameter :: maxreinit=10
  integer , parameter :: reinit_freq=10
  
end module levelset


! ==================================== !
! Initialization of the levelset field !
! ==================================== !
subroutine levelset_init
  use levelset
  implicit none
  integer :: i,j
  
  ! ==========================================
  ! ========= PARAMETERS TO MODIFY ===========
  ! ==========================================
  
  ! Loop over full field and create initial distance - Zalesak
  do j=1,ny
     do i=1,nx
        !G(i,j)=-sqrt(xm(i)**2+ym(j)**2)+0.15
        !if ((xm(i)**2+ym(j)**2).lt.0.15) G(i,j)=-G(i,j)
        !G(i,j)=init_zalesak((/xm(i),ym(j),0.0_WP/),(/0.0_WP,0.25_WP,0.0_WP/),0.15_WP,0.05_WP,0.25_WP)
        G(i,j)=init_deformed((/xm(i),ym(j),0.0_WP/),0.25_WP,0.05_WP)
     end do
  end do
  
  ! ==========================================
  
  ! Calculate normal and curvature
  call levelset_normal
  call levelset_curvature
  
contains
  
  ! Zalesak disk initialization
  function init_zalesak(xyz,center,radius,width,height)
    implicit none
    real(WP) :: init_zalesak
    real(WP), dimension(3), intent(in) :: xyz,center
    real(WP), intent(in) :: radius,width,height
    real(WP) :: c,b,b1,b2,h1,h2
    c = radius-sqrt(sum((xyz-center)**2))
    b1 = center(1)-0.5_WP*width
    b2 = center(1)+0.5_WP*width
    h1 = center(2)-radius*cos(asin(0.5_WP*width/radius))
    h2 = center(2)-radius+height
    if     (c>=0.0_WP.and.xyz(1)<=b1.and.xyz(2)<=h2) then
       b = b1-xyz(1)
       init_zalesak = min(c,b)
    elseif (c>=0.0_WP.and.xyz(1)>=b2.and.xyz(2)<=h2) then
       b = xyz(1)-b2
       init_zalesak = min(c,b)
    elseif (c>=0.0_WP.and.xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)>=h2) then
       b = xyz(2)-h2
       init_zalesak = min(c,b)
    elseif (c>=0.0_WP.and.xyz(1)<=b1.and.xyz(2)>=h2) then
       b = sqrt(sum((xyz-(/b1,h2,0.0_WP/))**2))
       init_zalesak = min(c,b)
    elseif (c>=0.0_WP.and.xyz(1)>=b2.and.xyz(2)>=h2) then
       b = sqrt(sum((xyz-(/b2,h2,0.0_WP/))**2))
       init_zalesak = min(c,b)
    elseif (xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)<=h2.and.xyz(2)>=h1) then
       init_zalesak = -min(abs(xyz(1)-b1),abs(xyz(1)-b2),abs(xyz(2)-h2)) 
    elseif (xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)<=h1) then
       init_zalesak = -min(sqrt(sum((xyz-(/b1,h1,0.0_WP/))**2)),sqrt(sum((xyz-(/b2,h1,0.0_WP/))**2)))
    else
       init_zalesak = c  
    endif
  end function init_zalesak

  function init_deformed(xyz,r0,alpha)
    implicit none
    real(WP) :: init_deformed
    real(WP), dimension(3), intent(in) :: xyz
    real(WP), intent(in) :: r0,alpha
    real(WP) :: r
    r = r0*(1.0_WP+alpha*cos(2.0_WP*atan2(xyz(2),xyz(1))))
    init_deformed = r-(xyz(1)**2+xyz(2)**2)**0.5_WP
  end function init_deformed
  
end subroutine levelset_init



! ==================================== !
! Initialization of the levelset field !
! ==================================== !
subroutine levelset_step
  use levelset
  implicit none
  
  integer  :: i,j
  real(WP) :: conv
  
  ! Remember previous residual
  HG2=HG1; HG1=0.0_WP
  
  ! Levelset update
  conv=dt/d
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) then
           HG1(i,j)=-conv*(max(U(i+1,j),0.0_WP)*sum(qp(-1:+1)*G(i-1:i+1,j))+min(U(i+1,j),0.0_WP)*sum(qm(0:+2)*G(i  :i+2,j))) &
                &   +conv*(max(U(i  ,j),0.0_WP)*sum(qp(-1:+1)*G(i-2:i  ,j))+min(U(i  ,j),0.0_WP)*sum(qm(0:+2)*G(i-1:i+1,j))) &
                &   -conv*(max(V(i,j+1),0.0_WP)*sum(qp(-1:+1)*G(i,j-1:j+1))+min(V(i,j+1),0.0_WP)*sum(qm(0:+2)*G(i,j  :j+2))) &
                &   +conv*(max(V(i,j  ),0.0_WP)*sum(qp(-1:+1)*G(i,j-2:j  ))+min(V(i,j  ),0.0_WP)*sum(qm(0:+2)*G(i,j-1:j+1)))
        end if
     end do
  end do
  
  ! Adams-Bashforth
  G=G+HG1+ABcoeff*(HG1-HG2)
  
  ! Apply Neumann in all directions
  do j=1,ny
     G(-1:0,j)=G(1,j)
     G(nx+1:nx+2,j)=G(nx,j)
  end do
  do i=-1,nx+2
     G(i,-1:0)=G(i,1)
     G(i,ny+1:ny+2)=G(i,ny)
  end do
  
  ! Calculate normal and curvature
  call levelset_normal
  call levelset_curvature

  ! Reinitialize the distance field
  if (mod(ntime,reinit_freq).eq.0) call levelset_reinit
  
  return
end subroutine levelset_step


! =================================== !
! Reinitialize levelset to a distance !
! This uses a 3rd order WENO scheme   !
! =================================== !
subroutine levelset_reinit
  use levelset
  implicit none
  integer  :: i,j,it
  real(WP), dimension(-1:nx+2,-1:ny+2) :: dGdt,Gold
  real(WP), dimension(-2:+1) :: weno3_xm,weno3_ym
  real(WP), dimension(-1:+2) :: weno3_xp,weno3_yp
  real(WP) :: dt_reinit,Gp,Gm,Gsign,dGdx_p,dGdx_m,dGdy_p,dGdy_m
  !real(WP) :: Gx,Gy,norm
  
  ! Remember G
  Gold=G
  
  ! Loop in pseudo-time
  do it=1,maxreinit

     ! Calculate dGdt
     do j=1,ny
        do i=1,nx
           if (mask(i,j).eq.0) then
              ! Naive centered attempt
              ! Magnitude of gradient
              !Gx=(G(i+1,j)-G(i-1,j))/(2.0_WP*d)
              !Gy=(G(i,j+1)-G(i,j-1))/(2.0_WP*d)
              !norm=sqrt(Gx**2+Gy**2)
              ! RHS of equation
              !dGdt(i,j)=sign(1.0_WP-norm,Gold(i,j))
              
              ! Prepare WENO3 coefficients
              call weno3_coeff
              
              ! Use WENO3 to calculate grad(G)
              dGdx_m=sum(weno3_xm(-2:+1)*G(i-2:i+1,j))
              dGdx_p=sum(weno3_xp(-1:+2)*G(i-1:i+2,j))
              dGdy_m=sum(weno3_ym(-2:+1)*G(i,j-2:j+1))
              dGdy_p=sum(weno3_yp(-1:+2)*G(i,j-1:j+2))
              
              ! Godunov gradient magnitude
              Gp=sqrt(max(max(dGdx_m,0.0_WP)**2,min(dGdx_p,0.0_WP)**2)+max(max(dGdy_m,0.0_WP)**2,min(dGdy_p,0.0_WP)**2))
              Gm=sqrt(max(min(dGdx_m,0.0_WP)**2,max(dGdx_p,0.0_WP)**2)+max(min(dGdy_m,0.0_WP)**2,max(dGdy_p,0.0_WP)**2))

              ! Smooth sign
              Gsign=Gold(i,j)/sqrt(Gold(i,j)**2+d**2+1.0e-9_WP)

              ! Compute dGdt
              dGdt(i,j)=-(max(Gsign,0.0_WP)*(Gp-1.0_WP)+min(Gsign,0.0_WP)*(Gm-1.0_WP))
              
           end if
        end do
     end do
     
     ! Update G
     dt_reinit=0.5_WP*d
     G=G+dt_reinit*dGdt
     
     ! Apply Neumann in all directions
     do j=1,ny
        G(-1:0,j)=G(1,j)
        G(nx+1:nx+2,j)=G(nx,j)
     end do
     do i=-1,nx+2
        G(i,-1:0)=G(i,1)
        G(i,ny+1:ny+2)=G(i,ny)
     end do
     
  end do

contains
  
  ! ================================================= !
  ! Set the coefficients for Hamilton-Jacobi operator !
  ! from Jiang and Peng, JCP, Vol. 21 No. 6, pp 2126  !
  ! ================================================= !
  subroutine weno3_coeff
    implicit none
    real(WP) :: r_minus,r_plus,w_minus,w_plus
    real(WP), parameter :: eps=1.0e-9_WP
    ! X DIRECTION ============================================
    ! Direction x - left biased stencil
    r_minus=(eps+(G(i  ,j)-2.0_WP*G(i-1,j)+G(i-2,j))**2.0_WP)&
         & /(eps+(G(i+1,j)-2.0_WP*G(i  ,j)+G(i-1,j))**2.0_WP)
    w_minus=1.0_WP/(1.0_WP+2.0_WP*r_minus**2.0_WP)
    weno3_xm(-2)=0.5_WP/d*(       +1.0_WP*w_minus)
    weno3_xm(-1)=0.5_WP/d*(-1.0_WP-3.0_WP*w_minus)
    weno3_xm( 0)=0.5_WP/d*(       +3.0_WP*w_minus)
    weno3_xm(+1)=0.5_WP/d*(+1.0_WP-1.0_WP*w_minus)
    ! Direction x - right biased stencil
    r_plus=(eps+(G(i+2,j)-2.0_WP*G(i+1,j)+G(i  ,j))**2.0_WP)&
         &/(eps+(G(i+1,j)-2.0_WP*G(i  ,j)+G(i-1,j))**2.0_WP)
    w_plus=1.0_WP/(1.0_WP+2.0_WP*r_plus**2.0_WP)
    weno3_xp(-1)=0.5_WP/d*(-1.0_WP+1.0_WP*w_plus)
    weno3_xp( 0)=0.5_WP/d*(       -3.0_WP*w_plus)
    weno3_xp(+1)=0.5_WP/d*(+1.0_WP+3.0_WP*w_plus)
    weno3_xp(+2)=0.5_WP/d*(       -1.0_WP*w_plus)
    
    ! Y DIRECTION ============================================
    ! Direction y - left biased stencil
    r_minus=(eps+(G(i,j  )-2.0_WP*G(i,j-1)+G(i,j-2))**2.0_WP)&
         & /(eps+(G(i,j+1)-2.0_WP*G(i,j  )+G(i,j-1))**2.0_WP)
    w_minus=1.0_WP/(1.0_WP+2.0_WP*r_minus**2.0_WP)
    weno3_ym(-2)=0.5_WP/d*(       +1.0_WP*w_minus)
    weno3_ym(-1)=0.5_WP/d*(-1.0_WP-3.0_WP*w_minus)
    weno3_ym( 0)=0.5_WP/d*(       +3.0_WP*w_minus)
    weno3_ym(+1)=0.5_WP/d*(+1.0_WP-1.0_WP*w_minus)
    ! Direction y - right biased stencil
    r_plus=(eps+(G(i,j+2)-2.0_WP*G(i,j+1)+G(i,j  ))**2.0_WP)&
         &/(eps+(G(i,j+1)-2.0_WP*G(i,j  )+G(i,j-1))**2.0_WP)
    w_plus=1.0_WP/(1.0_WP+2.0_WP*r_plus**2.0_WP)
    weno3_yp(-1)=0.5_WP/d*(-1.0_WP+1.0_WP*w_plus)
    weno3_yp( 0)=0.5_WP/d*(       -3.0_WP*w_plus)
    weno3_yp(+1)=0.5_WP/d*(+1.0_WP+3.0_WP*w_plus)
    weno3_yp(+2)=0.5_WP/d*(       -1.0_WP*w_plus)
    return
  end subroutine weno3_coeff
  
end subroutine levelset_reinit


! ===================================== !
! Calculation of levelset normal vector !
! ===================================== !
subroutine levelset_normal
  use levelset
  implicit none
  integer  :: i,j
  real(WP) :: norm
  
  ! Normal update
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) then
           ! Levelset gradient
           normx(i,j)=(G(i+1,j)-G(i-1,j))/(2.0_WP*d)
           normy(i,j)=(G(i,j+1)-G(i,j-1))/(2.0_WP*d)
           ! Normalize it
           norm=sqrt(normx(i,j)**2+normy(i,j)**2)
           normx(i,j)=normx(i,j)/(norm+epsilon(1.0_WP))
           normy(i,j)=normy(i,j)/(norm+epsilon(1.0_WP))
        end if
     end do
  end do
  
  ! Apply Neumann in all directions
  do j=1,ny
     normx(-1:0,j)=normx(1,j)
     normx(nx+1:nx+2,j)=normx(nx,j)
     normy(-1:0,j)=normy(1,j)
     normy(nx+1:nx+2,j)=normy(nx,j)
  end do
  do i=-1,nx+2
     normx(i,-1:0)=normx(i,1)
     normx(i,ny+1:ny+2)=normx(i,ny)
     normy(i,-1:0)=normy(i,1)
     normy(i,ny+1:ny+2)=normy(i,ny)
  end do
  
  return
end subroutine levelset_normal


! ================================= !
! Calculation of levelset curvature !
! ================================= !
subroutine levelset_curvature
  use levelset
  implicit none
  integer  :: i,j
  real(WP) :: Gx,Gy,Gxx,Gyy,Gxy
  
  ! Curvature update
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) then
           ! 2nd order FD compact scheme - 9 pt stencil
           Gx =(G(i+1,j)-G(i-1,j))/(2.0_WP*d)
           Gy =(G(i,j+1)-G(i,j-1))/(2.0_WP*d)
           Gxx=(G(i+1,j)-2.0_WP*G(i,j)+G(i-1,j))/d**2
           Gyy=(G(i,j+1)-2.0_WP*G(i,j)+G(i,j-1))/d**2
           Gxy=((G(i+1,j+1)-G(i+1,j-1))/(2.0_WP*d)-(G(i-1,j+1)-G(i-1,j-1))/(2.0_WP*d))/(2.0_WP*d)
           curv(i,j)=-(+Gx**2*Gyy-2.0_WP*Gx*Gy*Gxy+Gy**2*Gxx)/((Gx**2+Gy**2)**(1.5_WP)+epsilon(1.0_WP))
        end if
        curv(i,j)=1.0_WP/0.15_WP
     end do
  end do
  
  ! Apply Neumann in all directions
  do j=1,ny
     curv(-1:0,j)=curv(1,j)
     curv(nx+1:nx+2,j)=curv(nx,j)
  end do
  do i=-1,nx+2
     curv(i,-1:0)=curv(i,1)
     curv(i,ny+1:ny+2)=curv(i,ny)
  end do

  return
end subroutine levelset_curvature

subroutine levelset_jump
   use levelset

   implicit none
   integer i,j
 
 do i=1,nx
    do j=1,ny
       jcx(i,j)=0.0_WP
       jcy(i,j)=0.0_WP
 
       if ((mask(i+1,j  ).eq.1).or. &
       &   (mask(i-1,j  ).eq.1).or. &
       &   (mask(i  ,j-1).eq.1).or. &
       &   (mask(i  ,j+1).eq.1)) cycle
 
       if (G(i+1,j  )*G(i,j).lt.0) then 
          jcx(i,j)=sigma*curv(i,j)
          if (G(i+1,j).lt.G(i,j)) jcx(i,j)=-jcx(i,j)
       end if
 
       if (G(i-1,j  )*G(i,j).lt.0) then 
          jcx(i,j)=sigma*curv(i,j)
          if (G(i-1,j).lt.G(i,j)) jcx(i,j)=-jcx(i,j)
       end if
 
       if (G(i  ,j+1)*G(i,j).lt.0) then 
          jcy(i,j)=sigma*curv(i,j)
          if (G(i,j+1).lt.G(i,j)) jcy(i,j)=-jcy(i,j)
       end if
 
       if (G(i  ,j-1)*G(i,j).lt.0) then 
          jcy(i,j)=sigma*curv(i,j)
          if (G(i,j-1).lt.G(i,j)) jcy(i,j)=-jcy(i,j)
       end if
 
       div(i,j) = div(i,j)+(jcx(i,j)+jcy(i,j))*dt/d**2
    end do
 end do



   end subroutine levelset_jump
