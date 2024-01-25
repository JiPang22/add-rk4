program a
implicit none
integer :: i,j
integer,parameter :: imax=50000 
real :: h,t,x,v,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,sumi,sumr,om,z1,z2,u1,u2
real, dimension(imax) :: xt,noise
!!!!!!!!!!!!!!!Box-Muller Method로 정규분포난수생성!!!!!!!!!!!!!
call random_seed()  ! 시드 설정
do i = 1, imax / 2
!!!!!!!!!!!!!!!!!! 두 개의 독립적인 균등분포 난수 생성!!!!!!!!!!!!!
call random_number(u1)
call random_number(u2)

! !!!!!!!!!!!!!!!!!Box-Muller 변환 수행!!!!!!!!!!!!!!!
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)

!!!!!!!!!!!!!! 생성된 정규 분포 난수 저장!!!!!!!!!!!!!
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do
!!!!!!!!!!!!!!!난수 생성 끝!!!!!!!!!!!!!!!!!!

open(1,file='xt')
!!!!!!!!초기조건!!!!!!!!!!!!
t=0.;x=0.;v=1.;h=1.e-2
write(1,*) v,x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! rk4 method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,imax
k1x=fx(t,x,v);k1v=fv(t,x,v)+0.1*noise(i)

K2x=fx(t+h/2.,x+h*k1x/2.,v+h*k1v/2.);k2v=fv(t+h/2.,x+h*k1x/2.,v+h*k1v/2.)+0.1*noise(i)

k3x=fx(t+h/2.,x+h*k2x/2.,v+h*k2v/2.);k3v=fv(t+h/2.,x+h*k2x/2.,v+h*k2v/2.)+0.1*noise(i)

k4x=fx(t+h,x+h*k3x,v+h*k3v);k4v=fv(t+h,x+h*k3x,v+h*k3v)+0.1*noise(i)

t=t+h;x=x+h*(k1x+2.*k2x+2.*k3x+k4x)/6.;v=v+h*(k1v+2.*k2v+2.*k3v+k4v)/6.
write(1,*) v,x
xt(i)=x
end do
open(2,file='wx')

!DFT
do j=1,100
om=6.28*j/(imax*h)

!reset
t=0.
sumr=0.
sumi=0.

!sum
do i=1,imax
t=i*h
sumr=sumr+xt(i)*cos(om*t)*h
sumi=sumr-xt(i)*sin(om*t)*h
enddo
write(2,*) om, sqrt(sumi**2+sumr**2)
enddo

!!!!!!!!!!!!!!!!함수 정의!!!!!!!!!!!!!!!!!!!!
contains
real function fx(t,x,v)
real t,x,v
fx=v
end function fx

!!!!!!!!!!!!!!!!!운동방정식!!!!!!!!!!!!!!!!!
real function fv(t,x,v)
real t,x,v,om

om=0.9 ! rad*kHz
fv=-0.1*v-x+cos(om*t)+sin(om*t)
end function fv
!!!!!!!!!!!!!!!함수 정의 끝!!!!!!!!!!!!!!!!!!!!
end program 
