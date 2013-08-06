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
 call randomize_randomness
 open(unit=12,file='tmp/E.dat')
 open(unit=14,file='tmp/H.dat')
 open(unit=16,file='tmp/Psi.dat')
! open(unit=20,file='tmp/Evec.dat')
 call dump_and_plot(0)
 cont=1
 do i=1,nstep
  call propagate_H
  call propagate_E
  call plot_or_not(i,cont)
 end do
 close(12)
 close(14)
 close(16)
! close(20)
end subroutine propagation_cycle
subroutine plot_or_not(i,cont)
implicit none
integer :: i,cont
  if(cont==nwrite)then
   call dump_and_plot(i)
   cont=1
  else
   cont=cont+1
  endif
end subroutine
subroutine dump_and_plot(i)
implicit none
integer :: i
character(len=4):: flenumber,charsnapshots
 call int_to_char(flenumber,i/nwrite)
 call int_to_char(charsnapshots,snapshots)
 call dump_out(i)
 call system('gnuplot -e "idx='//flenumber//';snapshots='//charsnapshots//';" plot.p')
 rewind(12)
 rewind(14)
 rewind(16)
! rewind(20)
! call system('rm tmp/*.dat')
end subroutine dump_and_plot
subroutine dump_out(j)
implicit none
real :: xran,yran,tol,eex,eez,edx,edz,check,norm
!integer :: iran,kran,counter
integer :: i,k,j
! counter=0
! do
!  call random_number(xran)
!  call random_number(yran)
!  iran=ceiling(xran*nx)
!  kran=ceiling(yran*nz)
!  tol=2*epsilon(xran)
!  if( abs(H(iran,kran)).ge.tol )then
!  if( (abs(Ex(iran,kran)).ge.tol).or.&
!      (abs(Ez(iran,kran)).ge.tol) )then
!write(*,*)'db',counter,iran,kran,Ex(iran,kran)
!   eex=iran*dx
!   eez=kran*dz
!   edx=(Ex(iran,kran+1)+Ex(iran,kran))/2
!   edz=(Ez(iran+1,kran)+Ez(iran,kran))/2
!   norm=sqrt(edx**2.+edz**2.)
!   write(20,*)eex,eez,edx*dx,edz*dz
!   counter=counter+1
!  endif
!  if(counter.ge.nfield)exit
! end do
 write(*,*)'writing step',j,'of',nstep
! do k=0,nz+1,4
!  do i=0,nx+1,4
!   if( (i.eq.0.or.i.eq.(nx+1))&
!   .or.(k.eq.0.or.k.eq.(nz+1)) )then
!    edx=0
!    edz=0
!   else
!    edx=(Ex(i-1,k)+Ex(i-1,k-1))/2.
!    edz=(Ez(i,k-1)+Ez(i-1,k-1))/2.
!   endif
!   norm=sqrt(edx**2.+edz**2.)
!   write(20,*)i*dx,k*dz,edx*dx/norm,edz*dz/norm
!  end do
!  write(20,*)
! end do
! write(20,*)
! write(20,*)
do k=0,nz+1,4
 do i=0,nx+1,4
  if( (i.eq.0.or.i.eq.(nx+1))&
  .or.(k.eq.0.or.k.eq.(nz+1)) )then
   edx=0
   edz=0
  else
   edx=(Ex(i-1,k)+Ex(i-1,k-1))/2.
   edx=(Ez(i,k-1)+Ez(i-1,k-1))/2.
   write(12,*)i*dx,k*dz,(edx**2+edz**2)**0.5
  endif
 end do
 write(12,*)
end do
do k=1,nz,4
 do i=1,nx,4
   write(14,*)i*dx,k*dz,H(i,k)
 end do
 write(14,*)
end do
do k=nz/2-nz/4,nz/2+nz/4,1
 do i=material_limits(1)-20,nx,1
  if(inside_metal(i,k).eq..true.)then
   write(16,*)i*dx,k*dz,sqrt(psix(i,k)**2+psiz(i,k)**2)
  else
   write(16,*)i*dx,k*dz,0.0
  endif
 end do
 write(16,*)
end do
end subroutine dump_out
subroutine int_to_char(tmp,n)
implicit none
 integer n
 character(len=4) :: tmp
 write(tmp,'(I4)')n
 tmp = adjustl(tmp)
end subroutine int_to_char
end program fdtd
