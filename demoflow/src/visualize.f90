module visualize
  use demoflow
  implicit none
  
  ! Ensight output information
  real(WP) :: nviz
  real(SP), dimension(:), allocatable :: tviz
  
  ! Single-precision buffer for output
  real(SP), dimension(1:nx,1:ny) :: rarray
  
end module visualize


! ======================== !
! Initialize visualization !
! ======================== !
subroutine visualize_init
  use visualize
  implicit none
  integer  :: unit,i,j,ibuffer,ierr
  real(SP) :: rbuffer
  character(len=80) :: filename,cbuffer
  integer,  dimension(1:nx+1,1:ny+1,2) :: iblank
  
  ! Initialize visualization flag
  nviz=-1.0_WP
  
  ! Create directories
  call system('mkdir -p viz')
  call system('mkdir -p viz/V')
  call system('mkdir -p viz/P')
  call system('mkdir -p viz/phi')
  call system('mkdir -p viz/S')
  call system('mkdir -p viz/div')
  if (lpttrack.eq.1) call system('mkdir -p viz/particles')
  
  ! Output problem geometry
  filename='viz/geometry'
  open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
  cbuffer='C Binary';                          write(unit) cbuffer
  cbuffer='Ensight Gold Geometry File';        write(unit) cbuffer
  cbuffer='Structured Geometry from DemoFlow'; write(unit) cbuffer
  cbuffer='node id off';                       write(unit) cbuffer
  cbuffer='element id off';                    write(unit) cbuffer
  cbuffer='part';                              write(unit) cbuffer
  ibuffer=1;                                   write(unit) ibuffer
  cbuffer='Complete geometry';                 write(unit) cbuffer
  cbuffer='block uniform iblanked';            write(unit) cbuffer
  ibuffer=nx+1;                                write(unit) ibuffer
  ibuffer=ny+1;                                write(unit) ibuffer
  ibuffer=2;                                   write(unit) ibuffer ! This should be set to 1
  rbuffer=x(1);                                write(unit) rbuffer
  rbuffer=y(1);                                write(unit) rbuffer
  rbuffer=-0.5_SP*d;                           write(unit) rbuffer
  rbuffer=d;                                   write(unit) rbuffer
  rbuffer=d;                                   write(unit) rbuffer
  rbuffer=d;                                   write(unit) rbuffer ! This should really be 0
  ! Add blanking: all nodes assumed exterior at first
  iblank=0
  do j=1,ny
     do i=1,nx
        if (mask(i,j).eq.0) then
           iblank(i:i+1,j:j+1,:)=1
        end if
     end do
  end do
  write(unit) iblank
  close(unit)
  
  ! Ouput initial case file
  filename='viz/demoflow.case'
  open(newunit=unit,file=filename,form='formatted',iostat=ierr,status='replace')
  cbuffer='FORMAT';             write(unit,'(a80)') cbuffer
  cbuffer='type: ensight gold'; write(unit,'(a80)') cbuffer
  cbuffer='GEOMETRY';           write(unit,'(a80)') cbuffer
  cbuffer='model: geometry';    write(unit,'(a80)') cbuffer
  close(unit)
  
  ! Allocate tviz
  allocate(tviz(ceiling(maxtime/viztime)+1)); tviz=0.0_SP
  
  ! Perform first visualization dump
  call visualize_dump
  
  return
end subroutine visualize_init


! ============================================= !
! Output simulation data in Ensight Gold format !
! ============================================= !
subroutine visualize_dump
  use visualize
  implicit none
  integer  :: unit,i,j,ibuffer,ierr
  real(SP) :: rbuffer
  character(len=80) :: filename,cbuffer,buff
  
  ! Output visualization data
  if (floor(time/viztime).ne.nviz) then
     
     ! New file index
     nviz=floor(time/viztime)
     write(buff,'(i6.6)') int(nviz)
     tviz(int(nviz)+1)=time
     
     ! Create velocity data file
     filename='viz/V/V.'//trim(adjustl(buff))
     open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
     cbuffer='V';     write(unit) cbuffer
     cbuffer='part';  write(unit) cbuffer
     ibuffer=1;       write(unit) ibuffer
     cbuffer='block'; write(unit) cbuffer
     do j=1,ny
        do i=1,nx
           rarray(i,j)=0.5_WP*(U(i,j)+U(i+1,j))
        end do
     end do
     write(unit) rarray
     do j=1,ny
        do i=1,nx
           rarray(i,j)=0.5_WP*(V(i,j)+V(i,j+1))
        end do
     end do
     write(unit) rarray
     rarray=0.0_SP
     write(unit) rarray
     close(unit)
     
     ! Create pressure data file
     filename='viz/P/P.'//trim(adjustl(buff))
     open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
     cbuffer='P';     write(unit) cbuffer
     cbuffer='part';  write(unit) cbuffer
     ibuffer=1;       write(unit) ibuffer
     cbuffer='block'; write(unit) cbuffer
     rarray=P(1:nx,1:ny)
     write(unit) rarray
     close(unit)

     ! Create levelset data file
     filename='viz/phi/phi.'//trim(adjustl(buff))
     open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
     cbuffer='phi';     write(unit) cbuffer
     cbuffer='part';  write(unit) cbuffer
     ibuffer=1;       write(unit) ibuffer
     cbuffer='block'; write(unit) cbuffer
     rarray=phi(1:nx,1:ny)
     write(unit) rarray
     close(unit)
     
     ! Create divergence data file
     filename='viz/div/div.'//trim(adjustl(buff))
     open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
     cbuffer='div';   write(unit) cbuffer
     cbuffer='part';  write(unit) cbuffer
     ibuffer=1;       write(unit) ibuffer
     cbuffer='block'; write(unit) cbuffer
     rarray=div(1:nx,1:ny)
     write(unit) rarray
     close(unit)
     
     ! Output case file
     filename='viz/demoflow.case'
     open(newunit=unit,file=filename,form='formatted',iostat=ierr,status='replace')
     cbuffer='FORMAT';                                   write(unit,'(a80)') cbuffer
     cbuffer='type: ensight gold';                       write(unit,'(a80)') cbuffer
     cbuffer='GEOMETRY';                                 write(unit,'(a80)') cbuffer
     cbuffer='model: geometry';                          write(unit,'(a80)') cbuffer
     if (lpttrack.eq.1) cbuffer='measured: 1 particles/particles.******'; write(unit,'(a80)') cbuffer
     cbuffer='VARIABLE';                                 write(unit,'(a80)') cbuffer
     cbuffer='vector per element: 1 V V/V.******';       write(unit,'(a80)') cbuffer
     cbuffer='scalar per element: 1 P P/P.******';       write(unit,'(a80)') cbuffer
     cbuffer='scalar per element: 1 phi phi/phi.******';       write(unit,'(a80)') cbuffer
     cbuffer='scalar per element: 1 div div/div.******'; write(unit,'(a80)') cbuffer
     if (lpttrack.eq.1) cbuffer='vector per measured node: 1 PartVel particles/velocity.******'; write(unit,'(a80)') cbuffer
     cbuffer='TIME';                                     write(unit,'(a80)') cbuffer
     cbuffer='time set: 1';                              write(unit,'(a80)') cbuffer
     cbuffer='number of steps:';                   write(unit,'(a16,x,i12)') cbuffer,int(nviz)+1
     cbuffer='filename start number: 0';                 write(unit,'(a80)') cbuffer
     cbuffer='filename increment: 1';                    write(unit,'(a80)') cbuffer
     cbuffer='time values:';   write(unit,'(a12,x,10000000(3(ES12.5,x),/))') cbuffer,tviz(1:int(nviz)+1)
     close(unit)

     ! Create particles datafiles
     if (lpttrack.eq.1) then

        ! Particle positions
        filename='viz/particles/particles.'//trim(adjustl(buff))
        open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
        cbuffer='C Binary';                          write(unit) cbuffer
        cbuffer='Particle positions from demoflow';  write(unit) cbuffer
        cbuffer='particle coordinates';              write(unit) cbuffer
        ibuffer=max(np,1);                           write(unit) ibuffer
        do i=1,np ! Does not conform to Ensight's specifications (Ensight is wrong here)
           ibuffer=i;                                write(unit) ibuffer
        end do
        do i=1,np
           rbuffer=xp(i);                            write(unit) rbuffer
           rbuffer=yp(i);                            write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer ! Need a z position
        end do
        if (np.eq.0) then ! Zero particle requires special handling (Ensight is wrong again)
           ibuffer=1;                                write(unit) ibuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer ! Need a z position
        end if
        close(unit)

        ! Particle velocity
        filename='viz/particles/velocity.'//trim(adjustl(buff))
        open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
        cbuffer='Particle velocity';                 write(unit) cbuffer
        do i=1,np
           rbuffer=up(i);                            write(unit) rbuffer
           rbuffer=vp(i);                            write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer ! Need a w velocity
        end do
        if (np.eq.0) then ! Zero particle requires special handling (Ensight is wrong again)
           rbuffer=0.0_SP;                           write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer
           rbuffer=0.0_SP;                           write(unit) rbuffer ! Need a w velocity
        end if
        close(unit)
        
     end if
     
  end if
  
  return
end subroutine visualize_dump
