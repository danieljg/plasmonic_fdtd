module engine
real :: cz,cx,dx,dz,dt,c,pi
parameter (c=29979245800,pi=4.0*atan(1.0))
integer :: nstep,nwrite,nx,nz,nfield
real :: center_frequency
real, allocatable :: Ex(:,:),Ez(:,:),H(:,:),&
                     epsx(:,:),epsz(:,:)

contains
 integer function read_integer(k,variable)
 implicit none
 integer :: k, variable
 character(12) :: value
  call get_command_argument(k,value)
  read(value,*) variable
  read_integer=k+1
 end function read_integer

 integer function read_real(k,variable)
 implicit none
 integer :: k
 real    :: variable
 character(12) :: value
  call get_command_argument(k,value)
  read(value,*) variable
  read_real=k+1
 end function read_real

 subroutine impose_initial_conditions
 implicit none
 real :: wavelength=138.0e-7, width=1.0e-4, length=0.5e-4
 integer :: i,k,imax,kmax
 imax=width/dx
 kmax=length/dz
   Ex(500,500)=1.0
 end subroutine impose_initial_conditions

 subroutine randomize_randomness
 implicit none
 integer :: M=2,date_time(8),time_seed(2)
 character (len=12) :: real_clock(3)
  call date_and_time(real_clock(1), real_clock(2), &
                     real_clock(3), date_time)
  time_seed(1) = date_time(6)+date_time(7)+date_time(8)
  time_seed(2) = date_time(1)+date_time(2)+date_time(8)
  call random_seed(size=M)
  call random_seed(put=time_seed(1:M))
 end subroutine randomize_randomness

 subroutine initialize_and_read_parameters
 implicit none
  nstep=4
  nwrite=1
  nfield=2400
  if(command_argument_count().eq.0)then
   call read_default_values
  else
   call read_command_line_arguments
  endif
  call allocate_fields
  dt=min(dx,dz)
  dt=dx/(10*c)
  cx=c*dt/dx
  cz=c*dt/dz
 end subroutine initialize_and_read_parameters

 subroutine allocate_fields
 implicit none
  allocate(Ex(nx,nz+1))
  allocate(Ez(nx+1,nz))
  allocate(epsx(nx,nz+1))
  allocate(epsz(nx+1,nz))
  allocate(H(nx,nz))
  Ex(:,:)=0.0;Ez(:,:)=0.0;H(:,:)=0.0;
 end subroutine allocate_fields

 subroutine read_default_values
 implicit none
  nx=1000
  nz=1000
  dx=5e-7
  dz=dx
 end subroutine read_default_values

 subroutine read_command_line_arguments
 implicit none
 integer :: k=1
  k=read_integer(k,nx)
  k=read_integer(k,nz)
  k=read_real(k,dx)
  k=read_real(k,dz)
 end subroutine read_command_line_arguments

 subroutine build_physical_space
 implicit none
 epsx=1.0
 epsz=1.0
 end subroutine build_physical_space

 subroutine propagate_H
 implicit none
 integer :: i,k
 do i=1,nx
  do k=1,nz
   H(i,k)=H(i,k)+( cx*(Ez(i+1,k)-Ez(i,k))&
                  -cz*(Ex(i,k+1)-Ex(i,k)))
  end do
 end do
 end subroutine propagate_H

 subroutine propagate_E
 implicit none
 integer :: i,k
 ! these are needed...
 do k=2,nz
  Ex(1,k)=Ex(1,k)+cz*(H(1,k)-H(1,k-1))/epsx(1,k)
 end do
 do i=2,nx
  Ez(i,1)=Ez(i,1)+cx*(H(i,1)-H(i,1))/epsz(i,1)
 end do
 ! ...so that these share the loop
 do i=2,nx
  do k=2,nz
   Ex(i,k)=Ex(i,k)+cz*(H(i,k)-H(i,k-1))/epsx(i,k)
   Ez(i,k)=Ez(i,k)+cx*(H(i,k)-H(i-1,k))/epsz(i,k)
  end do
 end do
 end subroutine propagate_E

 subroutine welcome_message
 implicit none
 write(*,*)'Welcome to the Transverse-Magnetic solver'
 write(*,*)'Private release, Daniel Jimenez (2013)'
 end subroutine welcome_message
end module engine

program fdtd
use engine
implicit none

 call welcome_message
 call initialize_and_read_parameters
 call build_physical_space
 call impose_initial_conditions
 call propagation_cycle

contains
subroutine propagation_cycle
implicit none
integer :: i,cont
 open(unit=12,file='E.dat')
 call randomize_randomness
 call dump_out
 cont=1
 do i=1,nstep
  call propagate_H
  call propagate_E
  if(cont==nwrite)then
   call dump_out
   cont=1
  else
   cont=cont+1
  endif
 end do
 close(12)
end subroutine propagation_cycle

subroutine dump_out
implicit none
real :: xran,yran,tol,eex,eez,edx,edz
integer :: iran,kran,counter
integer :: i,k
 counter=0
! do
!  call random_number(xran)
!  call random_number(yran)
!  iran=ceiling(xran*nx)
!  kran=ceiling(yran*nz)
!  tol=4*epsilon(xran)
!  if( abs(H(iran,kran)).ge.tol )then
!  if( (abs(Ex(iran,kran)).ge.tol).or.&
!      (abs(Ez(iran,kran)).ge.tol) )then
!write(*,*)'db',counter,iran,kran,Ex(iran,kran)
!   eex=iran*dx
!   eez=kran*dz
!   edx=(Ex(iran,kran+1)+Ex(iran,kran))/2
!   edz=(Ez(iran+1,kran)+Ez(iran,kran))/2
!   write(12,*)eex,eez,4*edx*dx,4*edz*dz
!   counter=counter+1
!  endif
!  if(counter.ge.nfield)exit
! end do
do i=1,nx,1
 do k=1,nz,1
   edx=(Ex(i,k+1)+Ex(i,k))/2
   edz=(Ez(i+1,k)+Ez(i,k))/2
   write(12,*)i,k,edx,edz
 end do
end do
 write(12,*)
 write(12,*)
end subroutine dump_out
end program fdtd
