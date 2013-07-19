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
 open(unit=12,file='tmp/E.dat')
 open(unit=14,file='tmp/H.dat')
 open(unit=20,file='tmp/Evec.dat')
 call randomize_randomness
 call dump_out(0)
 cont=1
 do i=1,nstep
  call propagate_H
  call propagate_E
  if(cont==nwrite)then
   call dump_out(i)
   cont=1
  else
   cont=cont+1
  endif
 end do
 close(12)
 close(13)
 close(14)
 close(20)
end subroutine propagation_cycle

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
 do k=0,nz+1,4
  do i=0,nx+1,4
   if( (i.eq.0.or.i.eq.(nx+1))&
   .or.(k.eq.0.or.k.eq.(nz+1)) )then
    edx=0
    edz=0
   else
    edx=(Ex(i-1,k)+Ex(i-1,k-1))/2.
    edz=(Ez(i,k-1)+Ez(i-1,k-1))/2.
   endif
   norm=sqrt(edx**2.+edz**2.)
   write(20,*)i*dx,k*dx,edx*dx,edz*dx
  end do
  write(20,*)
 end do
 write(20,*)
 write(20,*)
do k=0,nz+1,4
 do i=0,nx+1,4
  if( (i.eq.0.or.i.eq.(nx+1))&
  .or.(k.eq.0.or.k.eq.(nz+1)) )then
   edx=0
   edz=0
  else
   edx=(Ex(i-1,k)+Ex(i-1,k-1))/2.
   edx=(Ez(i,k-1)+Ez(i-1,k-1))/2.
  endif
  write(12,*)i*dx,k*dx,(edx**2+edz**2)**0.5
 end do
 write(12,*)
end do
 write(12,*)
 write(12,*)
do k=1,nz,2
 do i=1,nx,2
   write(14,*)i*dx,k*dx,H(i,k)
 end do
 write(14,*)
end do
 write(14,*)
 write(14,*)
end subroutine dump_out
end program fdtd
