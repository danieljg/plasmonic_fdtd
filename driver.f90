module engine
real :: 
integer :: nstep,nwrite
real, allocatable :: Ex(:,:)Ez(:,:),H(:,:)
contains
 subroutine intialize_and_read_paramters
 implicit none
 nstep=128
 nwrite=32
 end subroutine intialize_and_read_parameters
 subroutine build_physical_space
 implicit none
 end subroutine build_physical_space
 subroutine impose_initial_conditions
 implicit none
 end subroutine impose_intial_conditions
 subroutine propagate_H
 implicit none
 end subroutine propagate_H
 subroutine propagate_Ex
 implicit none
 end subroutine propagate_Ex
 subroutine propagate_Ez
 implicit none
 end subroutine propagate_Ez

end module engine

program fdtd
use engine
implicit none

 call initialize_and_read_parameters
 call build_physical_space
 call impose_initial_conditions
 call propagation_cycle

contains
subroutine propagation_cycle
implicit none
integer :: i,cont
 call write_data
 cont=1
 do i=1,nstep
  if(cont==nwrite)then
   call write_data
   cont=1
  else
   cont=cont+1
  endif
  call propagate_H
  call propagate_Ex
  call propagate_Ez
 end do
end subroutine propagation_cycle
end program fdtd
