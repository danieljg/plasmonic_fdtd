module engine
real :: econst,hconst,cz,cx
integer :: nstep,nwrite,nx,nz
real, allocatable :: Ex(:,:),Ez(:,:),H(:,:)

contains
 integer function read_integer(k,variable)
 implicit none
 integer :: k, variable
 character(12) :: value
  call get_command_argument(k,value)
  read(value,*) variable
  read_integer=k+1
 end function read_integer

 subroutine initialize_and_read_parameters
 implicit none
  nstep=128
  nwrite=32
  if(command_argument_count().eq.0)then
   call read_default_values
  else
   call read_command_line_arguments
  endif
  allocate(Ex(nx+1,nz))
  allocate(Ez(nx,nz+1))
  allocate(H(nx,nz))
  Ex=0.0;Ez=0.0;H=0.0;
 end subroutine initialize_and_read_parameters

 subroutine read_default_values
 implicit none
  nx=1000
  nz=1000
 end subroutine read_default_values

 subroutine read_command_line_arguments
 implicit none
 integer :: k=1
  k=read_integer(k,nx)
  k=read_integer(k,nz)
 end subroutine read_command_line_arguments

 subroutine build_physical_space
 implicit none
 end subroutine build_physical_space

 subroutine impose_initial_conditions
 implicit none
 end subroutine impose_initial_conditions

 subroutine propagate_H
 implicit none
 integer :: i,k
 do i=1,nx
  do k=1,nz
   H(i,k)=H(i,k)+hconst*( cx*(Ez(i+1,k)-Ez(i,k))&
                         -cz*(Ex(i,k+1)-Ex(i,k)))
  end do
 end do
 end subroutine propagate_H

 subroutine propagate_E
 implicit none
 integer :: i,k
 ! these are needed...
 do k=2,nz
  Ex(1,k)=Ex(1,k)+econst*cz*(H(1,k)-H(1,k-1))
 end do
 do i=2,nx
  Ez(i,1)=Ez(i,1)+econst*cx*(H(i,1)-H(i,1))
 end do
 ! ...so that these share the loop
 do i=2,nx
  do k=2,nz
   Ex(i,k)=Ex(i,k)+econst*cz*(H(i,k)-H(i,k-1))
   Ez(i,k)=Ez(i,k)+econst*cx*(H(i,k)-H(i-1,k))
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
 call dump_out
 cont=1
 do i=1,nstep
  if(cont==nwrite)then
   call dump_out
   cont=1
  else
   cont=cont+1
  endif
  call propagate_H
  call propagate_E
 end do
 close(12)
end subroutine propagation_cycle
subroutine dump_out
implicit none

end subroutine dump_out
end program fdtd
