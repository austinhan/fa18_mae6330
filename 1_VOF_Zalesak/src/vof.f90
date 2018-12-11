module vofmod
  use demoflow
  use vof_lookup
  implicit none
  

  
  ! Normal vector
  real(WP), dimension(0:nx+1,0:ny+1) :: normx
  real(WP), dimension(0:nx+1,0:ny+1) :: normy
    
  ! PLIC distance
  real(WP), dimension(0:nx+1,0:ny+1) :: dist
  
  ! ==========================================
  ! ========= PARAMETERS TO MODIFY ===========
  ! ==========================================
  
  ! Surface tension coefficient
  real(WP), parameter :: sigma=1.0e-2_WP
  
  ! Phase densities
  real(WP), parameter :: rho_l=10.0_WP
  real(WP), parameter :: rho_g=1.0_WP

  ! Limits to VOF values
  real(WP), parameter :: VOFlo=1.0e-6_WP
  real(WP), parameter :: VOFhi=1.0_WP-VOFlo
  
  ! ==========================================

contains
  
  ! ==================== !
  ! Get indices of point !
  ! ==================== !
  function get_indices(pt,ind_in) result(ind)
    implicit none
    real(WP), dimension(3), intent(in) :: pt
    integer,  dimension(3), intent(in) :: ind_in
    integer,  dimension(3) :: ind
    ! X direction
    ind(1)=ind_in(1)
    do while (pt(1).gt.x(ind(1)+1)); ind(1)=ind(1)+1; end do
    do while (pt(1).lt.x(ind(1)  )); ind(1)=ind(1)-1; end do
    ! Y direction
    ind(2)=ind_in(2)
    do while (pt(2).gt.y(ind(2)+1)); ind(2)=ind(2)+1; end do
    do while (pt(2).lt.y(ind(2)  )); ind(2)=ind(2)-1; end do
    ! Z direction
    ind(3)=ind_in(3)
    return
  end function get_indices
  
  ! ===================== !
  ! Calculate sign of tet !
  ! ===================== !
  function tet_sign(vert) result(s)
    implicit none
    real(WP) :: s
    real(WP), dimension(3,4), intent(in) :: vert
    real(WP), dimension(3) :: a,b,c
    a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
    s=sign(1.0_WP,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP)
    return
  end function tet_sign
  
  ! ========================================= !
  ! Move point p1 to p2 according to velocity !
  ! ========================================= !
  function project(p1,mydt,myi,myj,myk) result(p2)
    implicit none
    real(WP), dimension(3) :: p2
    real(WP), dimension(3), intent(in) :: p1
    real(WP),               intent(in) :: mydt
    integer,                intent(in) :: myi,myj,myk
    real(WP), dimension(3) :: v1,v2,v3,v4
    v1=get_velocity(p1               ,myi,myj,myk)
    v2=get_velocity(p1+0.5_WP*mydt*v1,myi,myj,myk)
    v3=get_velocity(p1+0.5_WP*mydt*v2,myi,myj,myk)
    v4=get_velocity(p1+       mydt*v3,myi,myj,myk)
    p2=p1+mydt/6.0_WP*(v1+2.0_WP*v2+2.0_WP*v3+v4)
  end function project
  
  ! ================================================ !
  ! Interpolate velocity to pos near cell (i0,j0,k0) !
  ! ================================================ !
  function get_velocity(pos,i0,j0,k0) result(vel)
    implicit none
    real(WP), dimension(3) :: vel
    real(WP), dimension(3), intent(in) :: pos
    integer, intent(in) :: i0,j0,k0
    integer :: i,j
    real(WP) :: wx1,wy1
    real(WP) :: wx2,wy2
    real(WP), dimension(2,2) :: ww
    ! U velocity
    ! Find right i index in [1,nx]
    i=max(min(nx,i0),1)
    do while (pos(1).ge.x (i+1).and.i+1.le.nx); i=i+1; end do
    do while (pos(1).lt.x (i  ).and.i-1.ge.1 ); i=i-1; end do
    ! Find right j index in [0,ny]
    j=max(min(ny,j0),0)
    do while (pos(2).ge.ym(j+1).and.j+1.le.ny); j=j+1; end do
    do while (pos(2).lt.ym(j  ).and.j-1.ge.0 ); j=j-1; end do
    ! Prepare bilinear interpolation coefficients
    wx1=(pos(1)-x (i))/(x (i+1)-x (i)); wx2=1.0_WP-wx1
    wy1=(pos(2)-ym(j))/(ym(j+1)-ym(j)); wy2=1.0_WP-wy1
    ! Combine in a matrix
    ww(1,1)=wx1*wy1
    ww(2,1)=wx2*wy1
    ww(1,2)=wx1*wy2
    ww(2,2)=wx2*wy2
    ! Bilinear interpolation of U
    vel(1)=sum(ww*U(i:i+1,j:j+1))
    ! V velocity
    ! Find right i index in [0,nx]
    i=max(min(nx,i0),0)
    do while (pos(1).ge.xm(i+1).and.i+1.le.nx); i=i+1; end do
    do while (pos(1).lt.xm(i  ).and.i-1.ge.0 ); i=i-1; end do
    ! Find right j index in [1,ny]
    j=max(min(ny,j0),1)
    do while (pos(2).ge.y (j+1).and.j+1.le.ny); j=j+1; end do
    do while (pos(2).lt.y (j  ).and.j-1.ge.1 ); j=j-1; end do
    ! Prepare bilinear interpolation coefficients
    wx1=(pos(1)-xm(i))/(xm(i+1)-xm(i)); wx2=1.0_WP-wx1
    wy1=(pos(2)-y (j))/(y (j+1)-y (j)); wy2=1.0_WP-wy1
    ! Combine in a matrix
    ww(1,1)=wx1*wy1
    ww(2,1)=wx2*wy1
    ww(1,2)=wx1*wy2
    ww(2,2)=wx2*wy2
    ! Bilinear interpolation of V
    vel(2)=sum(ww*V(i:i+1,j:j+1))
    ! W velocity
    vel(3)=0.0_WP
    return
  end function get_velocity
  
  ! =============================================== !
  ! Cut tet by computational mesh to compute fluxes !
  ! =============================================== !
  recursive function tet2flux(tet,ind) result(fluxes)
    implicit none
    real(WP), dimension(3,4), intent(in) :: tet
    integer,  dimension(3,4), intent(in) :: ind
    real(WP), dimension(2) :: fluxes
    integer :: dir,cut_ind
    integer :: n,nn,case
    integer :: v1,v2
    real(WP), dimension(4) :: myd
    real(WP), dimension(3,8) :: vert
    integer,  dimension(3,8,2) :: vert_ind
    real(WP) :: mu
    real(WP), dimension(3,4) :: newtet
    integer,  dimension(3,4) :: newind
    
    ! Cut by x planes
    if      (maxval(ind(1,:))-minval(ind(1,:)).gt.0) then
       dir=1
       cut_ind=maxval(ind(1,:))
       myd(:)=tet(1,:)-x(cut_ind)
    ! Cut by y planes
    else if (maxval(ind(2,:))-minval(ind(2,:)).gt.0) then
       dir=2
       cut_ind=maxval(ind(2,:))
       myd(:)=tet(2,:)-y(cut_ind)
    ! Cut by interface and compute fluxes
    else 
       fluxes=tet2flux_plic(tet,ind(1,1),ind(2,1),ind(3,1))    
       return
    end if
    
    ! Find case of cut
    case=1+int(0.5_WP+sign(0.5_WP,myd(1)))+&
         2*int(0.5_WP+sign(0.5_WP,myd(2)))+&
         4*int(0.5_WP+sign(0.5_WP,myd(3)))+&
         8*int(0.5_WP+sign(0.5_WP,myd(4)))   
    
    ! Get vertices and indices of tet
    do n=1,4
       vert    ( : ,n  )=tet(:,n)
       vert_ind( : ,n,1)=ind(:,n)
       vert_ind( : ,n,2)=ind(:,n)
       vert_ind(dir,n,1)=min(vert_ind(dir,n,1),cut_ind-1) ! Enforce boundedness
       vert_ind(dir,n,2)=max(vert_ind(dir,n,1),cut_ind  )
    end do
    
    ! Create interpolated vertices on cut plane
    do n=1,cut_nvert(case)
       v1=cut_v1(n,case); v2=cut_v2(n,case)
       mu=min(1.0_WP,max(0.0_WP,-myd(v1)/(sign(abs(myd(v2)-myd(v1))+tiny(1.0_WP),myd(v2)-myd(v1)))))
       vert(:,4+n)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
       ! Get index for interpolated vertex
       vert_ind(:,4+n,1)=get_indices(vert(:,4+n),vert_ind(:,v1,1))
       ! Enforce boundedness
       vert_ind(:,4+n,1)=max(vert_ind(:,4+n,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       vert_ind(:,4+n,1)=min(vert_ind(:,4+n,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
       ! Set +/- indices in cut direction
       vert_ind(:,4+n,2)=vert_ind(:,4+n,1)
       vert_ind(dir,4+n,1)=cut_ind-1
       vert_ind(dir,4+n,2)=cut_ind
    end do
    
    ! Create new tets
    fluxes=0.0_WP
    do n=1,cut_ntets(case)
       do nn=1,4
          newtet(:,nn)=vert    (:,cut_vtet(nn,n,case))
          newind(:,nn)=vert_ind(:,cut_vtet(nn,n,case),cut_side(n,case))
       end do
       ! Cut by next plane
       fluxes=fluxes+tet2flux(newtet,newind)
    end do
    
    return
  end function tet2flux
  
  ! ========================= !
  ! Cut tet by PLIC interface !
  ! ========================= !
  function tet2flux_plic(tet,i,j,k) result(fluxes)
    implicit none
    real(WP), dimension(3,4), intent(in) :: tet
    integer, intent(in) :: i,j,k
    real(WP), dimension(2) :: fluxes
    integer :: n,case,v1,v2
    real(WP) :: mu,my_vol
    real(WP), dimension(4) :: myd
    real(WP), dimension(3,8) :: vert
    real(WP), dimension(3) :: a,b,c
    
    ! Cut by PLIC
    myd(:)=normx(i,j)*tet(1,:)+normy(i,j)*tet(2,:)-dist(i,j)
    
    ! Find cut case
    case=1+int(0.5_WP+sign(0.5_WP,myd(1)))+&
         2*int(0.5_WP+sign(0.5_WP,myd(2)))+&
         4*int(0.5_WP+sign(0.5_WP,myd(3)))+&
         8*int(0.5_WP+sign(0.5_WP,myd(4)))
    
    ! Copy vertices
    do n=1,4
       vert(:,n)=tet(:,n)
    end do
    
    ! Create interpolated vertices on cut plane
    do n=1,cut_nvert(case)
       v1=cut_v1(n,case); v2=cut_v2(n,case)
       mu=min(1.0_WP,max(0.0_WP,-myd(v1)/(sign(abs(myd(v2)-myd(v1))+tiny(1.0_WP),myd(v2)-myd(v1)))))
       vert(:,4+n)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
    end do
    
    ! Zero fluxes
    fluxes=0.0_WP
        
    ! Analyze gas tets
    do n=1,cut_nntet(case)-1
       ! Compute volume
       a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
       my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
       ! Quantities in gas phase
       fluxes(1)=fluxes(1)+my_vol
    end do
    
    ! Analyze liquid tets
    do n=cut_ntets(case),cut_nntet(case),-1
       ! Compute volume
       a=vert(:,cut_vtet(1,n,case))-vert(:,cut_vtet(4,n,case)); b=vert(:,cut_vtet(2,n,case))-vert(:,cut_vtet(4,n,case)); c=vert(:,cut_vtet(3,n,case))-vert(:,cut_vtet(4,n,case))
       my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
       ! Quantities in liquid phase
       fluxes(1)=fluxes(1)+my_vol
       fluxes(2)=fluxes(2)+my_vol
    end do
    
    return
  end function tet2flux_plic
  
end module vofmod


! ==================================== !
! Initialization of the levelset field !
! ==================================== !
subroutine vof_init
  use vofmod
  implicit none
  integer  :: i,j,ii,jj
  integer, parameter :: nfine=50
  
  ! ==========================================
  ! ========= PARAMETERS TO MODIFY ===========
  ! ==========================================
  
  ! Loop over full field and create initial distance - Zalesak
  VOF=0.0_WP
  do j=1,ny
     do i=1,nx
        ! Get Zalesak distance on submesh
        do jj=1,nfine
           do ii=1,nfine
              if (init_zalesak((/x(i)+(real(ii-1,WP)+0.5_WP)*d/real(nfine,WP),&
                   &             y(j)+(real(jj-1,WP)+0.5_WP)*d/real(nfine,WP),0.0_WP/),&
                   &           (/0.0_WP,0.25_WP,0.0_WP/),0.15_WP,0.05_WP,0.25_WP).gt.0.0_WP) then
                 VOF(i,j)=VOF(i,j)+1.0_WP/(nfine**2)
              end if
           end do
        end do
     end do
  end do
  
  ! ==========================================
  
  ! Calculate normal and PLIC
  call vof_normal
  call vof_plic
  
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
  
end subroutine vof_init


! ============================ !
! Advancement of the vof field !
! ============================ !
subroutine vof_step
  use vofmod
  implicit none
  
  integer  :: i,j,n
  real(WP), dimension(3,8) :: pt
  real(WP), dimension(3,4) :: tet
  integer,  dimension(3,4) :: ind
  real(WP), dimension( 2 ) :: flux
  
  ! Perform pure semi-Lagrangian update
  do j=1,ny
     do i=1,nx
        ! Skip cells without normals as surrogate for band
        if (normx(i,j)**2+normy(i,j)**2.lt.1.0e-6_WP) cycle
        ! Form cell and project backward in time
        pt(:,1)=[x(i+1),y(j  ),0.0_WP]; pt(:,1)=project(pt(:,1),-dt,i,j,0)
        pt(:,2)=[x(i+1),y(j  ),1.0_WP]; pt(:,2)=project(pt(:,2),-dt,i,j,0)
        pt(:,3)=[x(i+1),y(j+1),1.0_WP]; pt(:,3)=project(pt(:,3),-dt,i,j,0)
        pt(:,4)=[x(i+1),y(j+1),0.0_WP]; pt(:,4)=project(pt(:,4),-dt,i,j,0)
        pt(:,5)=[x(i  ),y(j  ),0.0_WP]; pt(:,5)=project(pt(:,5),-dt,i,j,0)
        pt(:,6)=[x(i  ),y(j  ),1.0_WP]; pt(:,6)=project(pt(:,6),-dt,i,j,0)
        pt(:,7)=[x(i  ),y(j+1),1.0_WP]; pt(:,7)=project(pt(:,7),-dt,i,j,0)
        pt(:,8)=[x(i  ),y(j+1),0.0_WP]; pt(:,8)=project(pt(:,8),-dt,i,j,0)
        ! Split in 6 tets
        flux=0.0_WP
        do n=1,6
           ! Build tet and corresponding indices
           tet(:,1)=pt(:,tet_map2(1,n)); ind(:,1)=get_indices(tet(:,1),[i,j,0])
           tet(:,2)=pt(:,tet_map2(2,n)); ind(:,2)=get_indices(tet(:,2),[i,j,0])
           tet(:,3)=pt(:,tet_map2(3,n)); ind(:,3)=get_indices(tet(:,3),[i,j,0])
           tet(:,4)=pt(:,tet_map2(4,n)); ind(:,4)=get_indices(tet(:,4),[i,j,0])
           ! Add corresponding flux
           flux=flux+tet_sign(tet)*tet2flux(tet,ind)
        end do
        ! Update VOF field
        VOF(i,j)=flux(2)/flux(1)
        if      (VOF(i,j).lt.VOFlo) then
           VOF(i,j)=0.0_WP
        else if (VOF(i,j).gt.VOFhi) then
           VOF(i,j)=1.0_WP
        end if
     end do
  end do
  
  ! Apply Neumann in all directions
  do j=1,ny
     VOF(0,j)=VOF(1,j)
     VOF(nx+1,j)=VOF(nx,j)
  end do
  do i=0,nx+1
     VOF(i,0)=VOF(i,1)
     VOF(i,ny+1)=VOF(i,ny)
  end do
  
  ! Calculate normal and PLIC
  call vof_normal
  call vof_plic
  
  return
end subroutine vof_step


! ====================================== !
! Calculation of VOF PLIC reconstruction !
! 2D version of Scardovelli's method     !
! ====================================== !
subroutine vof_plic
  use vofmod
  implicit none
  integer  :: i,j
  real(WP) :: alpha,mm1,mm2,factor,norm,tmp,VOFo,V0,m1,m2
  
  ! Update distance to PLIC plane
  do j=0,ny+1
     do i=0,nx+1
        ! Skip wall cells
        if (mask(i,j).eq.1) cycle
        ! Form plane equation variables, dealing with non-cubic cells
        mm1=normx(i,j)*d
        mm2=normy(i,j)*d
        ! Calculate normalization for mm
        norm=abs(mm1)+abs(mm2)
        ! Branching treatment of normal vs. no normal
        if (norm.lt.epsilon(1.0_WP)) then
           ! We do not have normals, return a large distance
           dist(i,j)=sign(huge(1.0_WP),VOF(i,j)-0.5d0)
        else
           ! We do have normals, finish PLIC calculation
           mm1=mm1/norm; mm2=mm2/norm
           ! Deal with negative mm
           factor=0.0_WP
           if (mm1.lt.0.0_WP) factor=factor+mm1
           if (mm2.lt.0.0_WP) factor=factor+mm2
           ! Deal with VOF>0.5
           VOFo=VOF(i,j)
           if (VOFo.gt.0.5_WP) VOFo=1.0_WP-VOFo
           ! Compute alpha
           m1=min(abs(mm1),abs(mm2))
           m2=max(abs(mm1),abs(mm2))
           V0=0.5_WP*m1/m2
           if (VOFo.lt.V0) then
              alpha=sqrt(2.0_WP*m1*m2*VOFo)
           else
              alpha=0.5_WP*m1+m2*VOFo
           end if
           ! Deal with VOF>0.5
           if (VOF(i,j).gt.0.5_WP) alpha=1.0_WP-alpha
           ! Adjust alpha if mm is negative
           alpha=alpha+factor
           ! Write plane with barycenter as origin inside cube
           alpha=alpha-mm1*0.5_WP-mm2*0.5_WP
           ! Make distance consistant with original normal
           dist(i,j)=alpha*norm
           ! Write plane in a global reference frame
           dist(i,j)=dist(i,j)+normx(i,j)*xm(i)+normy(i,j)*ym(j)
        end if
     end do
  end do
  
  return
end subroutine vof_plic


! ================================ !
! Calculation of VOF normal vector !
! Based on simple gradient...      !
! ================================ !
subroutine vof_normal
  use vofmod
  implicit none
  integer  :: i,j
  real(WP) :: norm
  real(WP), dimension(0:nx+1,0:ny+1) :: tmp

  ! First smooth VOF field
  call smoother
  
  ! Normal update using Young's method
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) then
           ! VOF gradient
           normx(i,j)=(tmp(i+1,j)-tmp(i-1,j))/(2.0_WP*d)
           normy(i,j)=(tmp(i,j+1)-tmp(i,j-1))/(2.0_WP*d)
           ! Normalize it
           norm=sqrt(normx(i,j)**2+normy(i,j)**2)
           if (norm.gt.0.0_WP) then
              normx(i,j)=-normx(i,j)/norm
              normy(i,j)=-normy(i,j)/norm
           else
              normx(i,j)=0.0_WP
              normy(i,j)=0.0_WP
           end if
        end if
     end do
  end do
  
  ! Apply Neumann in all directions
  do j=1,ny
     normx(0,j)=normx(1,j)
     normx(nx+1,j)=normx(nx,j)
     normy(0,j)=normy(1,j)
     normy(nx+1,j)=normy(nx,j)
  end do
  do i=0,nx+1
     normx(i,0)=normx(i,1)
     normx(i,ny+1)=normx(i,ny)
     normy(i,0)=normy(i,1)
     normy(i,ny+1)=normy(i,ny)
  end do
  
contains
  
  ! =========================================== !
  ! Smoothing of VOF to improve Young's normals !
  ! =========================================== !
  subroutine smoother
    implicit none
    integer :: i,j,ii,jj,count
    
    ! Perform local average
    tmp=0.0_WP
    do j=1,ny
       do i=1,nx
          ! Skip wall cells
          if (mask(i,j).eq.1) cycle
          ! Smooth local value
          count=0
          do jj=j-1,j+1
             do ii=i-1,i+1
                ! Skip wall cells
                if (mask(ii,jj).eq.1) cycle
                ! Average
                tmp(i,j)=tmp(i,j)+VOF(ii,jj)
                count=count+1
             end do
          end do
          tmp(i,j)=tmp(i,j)/real(count,WP)
       end do
    end do

    ! Handle BC
    do j=1,ny
       tmp(0,j)=tmp(1,j)
       tmp(nx+1,j)=tmp(nx,j)
    end do
    do i=0,nx+1
       tmp(i,0)=tmp(i,1)
       tmp(i,ny+1)=tmp(i,ny)
    end do
    
    return
  end subroutine smoother
  
end subroutine vof_normal

