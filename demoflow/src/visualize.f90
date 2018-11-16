module visualize
  use demoflow
  use levelset
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
  call system('mkdir -p viz/S')
  call system('mkdir -p viz/div')
  if (use_levelset) then
     call system('mkdir -p viz/G')
     call system('mkdir -p viz/norm')
     call system('mkdir -p viz/curv')
  end if
  
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
     cbuffer='VARIABLE';                                 write(unit,'(a80)') cbuffer
     cbuffer='vector per element: 1 V V/V.******';       write(unit,'(a80)') cbuffer
     cbuffer='scalar per element: 1 P P/P.******';       write(unit,'(a80)') cbuffer
     cbuffer='scalar per element: 1 div div/div.******'; write(unit,'(a80)') cbuffer
     if (use_levelset) then
        cbuffer='scalar per element: 1 G G/G.******'; write(unit,'(a80)') cbuffer
        cbuffer='vector per element: 1 norm norm/norm.******'; write(unit,'(a80)') cbuffer
        cbuffer='scalar per element: 1 curv curv/curv.******'; write(unit,'(a80)') cbuffer
     end if
     cbuffer='TIME';                                     write(unit,'(a80)') cbuffer
     cbuffer='time set: 1';                              write(unit,'(a80)') cbuffer
     cbuffer='number of steps:';                   write(unit,'(a16,x,i12)') cbuffer,int(nviz)+1
     cbuffer='filename start number: 0';                 write(unit,'(a80)') cbuffer
     cbuffer='filename increment: 1';                    write(unit,'(a80)') cbuffer
     cbuffer='time values:';   write(unit,'(a12,x,10000000(3(ES12.5,x),/))') cbuffer,tviz(1:int(nviz)+1)
     close(unit)

     ! Create levelset datafiles
     if (use_levelset) then

        ! Level set
        filename='viz/G/G.'//trim(adjustl(buff))
        open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
        cbuffer='G';     write(unit) cbuffer
        cbuffer='part';  write(unit) cbuffer
        ibuffer=1;       write(unit) ibuffer
        cbuffer='block'; write(unit) cbuffer
        rarray=G(1:nx,1:ny)
        write(unit) rarray
        close(unit)

        ! Level set
        filename='viz/curv/curv.'//trim(adjustl(buff))
        open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
        cbuffer='curv';  write(unit) cbuffer
        cbuffer='part';  write(unit) cbuffer
        ibuffer=1;       write(unit) ibuffer
        cbuffer='block'; write(unit) cbuffer
        rarray=curv(1:nx,1:ny)
        write(unit) rarray
        close(unit)
        
        ! Normal vector
        filename='viz/norm/norm.'//trim(adjustl(buff))
        open(newunit=unit,file=filename,form='unformatted',access='stream',iostat=ierr,status='replace')
        cbuffer='norm';  write(unit) cbuffer
        cbuffer='part';  write(unit) cbuffer
        ibuffer=1;       write(unit) ibuffer
        cbuffer='block'; write(unit) cbuffer
        rarray=normx(1:nx,1:ny)
        write(unit) rarray
        rarray=normy(1:nx,1:ny)
        write(unit) rarray
        rarray=0.0_SP
        write(unit) rarray
        close(unit)
        
     end if
     
  end if
  
  return
end subroutine visualize_dump
