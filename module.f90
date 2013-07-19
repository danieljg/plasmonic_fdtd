module engine
real :: cz,cx,dx,dz,dt,c,pi,czabc,cxabc,s
parameter (c=29979245800,pi=4.0*atan(1.0))
!parameter (c=1,pi=4.0*atan(1.0))
integer :: nstep,nwrite,snapshots,nx,nz,nfield,material_limits(4)
!THIS MUST MATCH WITH THE NUMBER OF SNAPSHOTS ON PLOT
parameter (snapshots=240,nwrite=13)
real :: wavelength,packetwidth,packetlength,plasma_freq,drude_gamma,&
        drude_sigma,xi_zero,delta_xi,egamma,inverse_k
real, allocatable :: Ex(:,:),Ez(:,:),H(:,:),&
                     epsx(:,:),epsz(:,:),&
                     psix(:,:),psiz(:,:),&
                     Exsav1(:),Ezsav1(:),&
                     Exsav2(:),Ezsav2(:)

contains
 subroutine read_default_values
 implicit none
  nx=360+3   !puntos en x
  nz=300+3   !puntos en z
  dx=10e-7    !espaciado en x
  dz=dx      !espaciado en z
  s=0.2      !parametro de estabilidad
  wavelength=400e-7*0.707  !longitud de onda central del paquete incidente
  packetwidth=3.0    !ancho del paquete (en longitudes de onda)
  packetlength=2.0   !largo del paquete 
  plasma_freq=1.36e16!frecuencia angular de plasma
  drude_gamma=3.6e13 !frecuencia de colisiones
  material_limits(1)=220
  material_limits(2)=280
  material_limits(3)=100
  material_limits(4)=200
 end subroutine read_default_values

 subroutine gaussian_packet(center_x,center_z,angle,wln,inwidth,inlength)
 implicit none
 integer, intent(in) :: center_x,center_z
 real, intent (in) :: angle, wln, inwidth, inlength
 real :: width, length, denomx, denomz
 real :: Np, wavevec(2),unit_vec(2),rex(2),rez(2),rh(2)
 integer :: i, k
  width=wln*inwidth
  length=wln*inlength
  Np=wln/dx
  wavevec(1)=2.0*pi*cos(angle)/wln
  wavevec(2)=2.0*pi*sin(angle)/wln
  denomz  = (width*cos(angle))**2+(length*sin(angle))**2
  denomx  = (width*sin(angle))**2+(length*cos(angle))**2
  do i=1,nx
   do k=1,nz
    rex(1)=i*dx
    rex(2)=(k+0.5)*dz
    rez(1)=(i+0.5)*dx
    rez(2)=k*dz
    rh(1)=(i+0.5+0.5*(1.0-s))*dx
    rh(2)=(k+0.5+0.5*(1.0-s))*dz
    Ex(i,k)= sin(angle)*cos(dot_product(wavevec,rex))*&
       exp(-(dx*center_x-rex(1))**2/denomx)*&
       exp(-(dz*center_z-rex(2))**2/denomz)
    Ez(i,k)=-cos(angle)*cos(dot_product(wavevec,rez))*&
       exp(-(dx*center_x-rez(1))**2/denomx)*&
       exp(-(dz*center_z-rez(2))**2/denomz)
    H(i,k) = cos(dot_product(wavevec,rh))*&
       exp(-(dx*center_x-rh(1))**2/denomx)*&
       exp(-(dz*center_z-rh(2))**2/denomz)
   end do
  end do
 end subroutine gaussian_packet

 subroutine impose_initial_conditions
 implicit none
  call gaussian_packet(nint(nx/3.),nint(nz/2.),0.,&
                         wavelength,packetwidth,packetlength)
 end subroutine impose_initial_conditions

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
  nstep=snapshots*nwrite
  nfield=5120
  if(command_argument_count().eq.0)then
   call read_default_values
  else
   call read_command_line_arguments
  endif
  call allocate_fields
  dt=min(dx,dz)
  dt=s*dt/c
  cx=c*dt/dx
  cz=c*dt/dz
  czabc=(c*dt-dz)/(c*dt+dz)
  cxabc=(c*dt-dx)/(c*dt+dx)
  egamma=exp(-drude_gamma*dt)
  xi_zero=-(plasma_freq/drude_gamma)**2*(1-egamma)
  delta_xi=-(plasma_freq*(1-egamma)/drude_gamma)**2
  drude_sigma=plasma_freq**2*drude_gamma/(4.0*pi)
  inverse_k=1+4*pi*xi_zero+4*pi*drude_sigma*dt
  inverse_k=1.0/inverse_k
  write(*,*)'delta t: ', dt
  write(*,*)'total simulation time: ', dt*nwrite*snapshots
 end subroutine initialize_and_read_parameters

 subroutine allocate_fields
 implicit none
  allocate(Ex(nx,nz+1))
  allocate(Exsav1(nx))
  allocate(Exsav2(nx))
  allocate(Ez(nx+1,nz))
  allocate(Ezsav1(nz))
  allocate(Ezsav2(nz))
  allocate(epsx(nx,nz+1))
  allocate(epsz(nx+1,nz))
  allocate(H(nx,nz))
  Ex(:,:)=0.0;Ez(:,:)=0.0;H(:,:)=0.0;
 end subroutine allocate_fields

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
 integer :: imin,imax,kmin,kmax
 epsx=1.0
 epsz=1.0
 imin=material_limits(1)
 imax=material_limits(2)
 kmin=material_limits(3)
 kmax=material_limits(4)
 !call slab_of_glass(imin,imax,kmin,kmax)
 call slab_of_metal(imin,imax,kmin,kmax)
 end subroutine build_physical_space

 subroutine slab_of_glass(imin,imax,kmin,kmax)
 implicit none
 integer, intent(in) :: imin,imax,kmin,kmax
 integer :: i,k
 do i=imin,imax
  do k=kmin,kmax
   epsx(i,k)=(1.55)**2.
   epsz(i,k)=(1.55)**2.
  end do
 end do 
 end subroutine slab_of_glass

 subroutine slab_of_metal(imin,imax,kmin,kmax)
 implicit none
 integer, intent(in) :: imin,imax,kmin,kmax
 allocate(psix(imin:imax,kmin:kmax))
 allocate(psiz(imin:imax,kmin:kmax))
 psix(:,:)=0.
 psiz(:,:)=0.
 end subroutine slab_of_metal

 function inside_metal(i,k)
 implicit none
 logical :: inside_metal
 integer, intent(in) :: i,k
  if(  ((i.gt.material_limits(1)).and.(i.lt.material_limits(2)))&
  .and.((k.gt.material_limits(3)).and.(k.lt.material_limits(4))) )then
  inside_metal=.true.
  else
  inside_metal=.false.
  endif
  ! hexagono!!!!!!!
  
 end function inside_metal

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
 Ezsav1(:)=Ez(2,:)
 Ezsav2(:)=EZ(nx,:)
 Exsav1(:)=Ex(:,2)
 Exsav2(:)=Ex(:,nz)
 ! these are needed...
 do k=2,nz
  Ex(1,k)=Ex(1,k)-cz*(H(1,k)-H(1,k-1))/epsx(1,k)
 end do
 do i=2,nx
  Ez(i,1)=Ez(i,1)+cx*(H(i,1)-H(i-1,1))/epsz(i,1)
 end do
 ! ...so that these share the loop...
 do i=2,nx
  do k=2,nz
   if(inside_metal(i,k).eq..true.)then
    Ex(i,k)=inverse_k*(Ex(i,k)-cz*(H(i,k)-H(i,k-1))+4*pi*psix(i,k))
    Ez(i,k)=inverse_k*(Ez(i,k)+cx*(H(i,k)-H(i-1,k))+4*pi*psiz(i,k))
    psix(i,k)=Ex(i,k)*delta_xi+egamma*psix(i,k)
    psiz(i,k)=Ez(i,k)*delta_xi+egamma*psix(i,k)
   else
    Ex(i,k)=Ex(i,k)-cz*(H(i,k)-H(i,k-1))/epsx(i,k)
    Ez(i,k)=Ez(i,k)+cx*(H(i,k)-H(i-1,k))/epsz(i,k)
   endif
  end do
 end do
 ! Absorbing boundary conditions Mur, first order
 do k=1,nz
  Ez(1,k)=Ezsav1(k)+czabc*(Ez(2,k)-Ez(1,k))
  Ez(nx+1,k)=Ezsav2(k)+czabc*(Ez(nx,k)-Ez(nx+1,k))
 end do
 do i=1,nx
  Ex(i,1)=Exsav1(i)+cxabc*(Ex(i,2)-Ex(i,1))
  Ex(i,nz+1)=Exsav2(i)+cxabc*(Ex(i,nz)-Ex(i,nz+1))
 end do
 end subroutine propagate_E

 subroutine welcome_message
 implicit none
 write(*,*)'Welcome to the Transverse-Magnetic solver'
 write(*,*)'Private release, Daniel Jimenez (2013)'
 write(*,*)'Total number of steps ',nwrite*snapshots
 end subroutine welcome_message
end module engine
