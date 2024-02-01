program ab
implicit none
integer :: i,j,k
integer, parameter :: imax=500000
real :: t,x,v,dx,dv,sumi,sumr,om,z1,z2,u1,u2,f,xbar,dxbar
real, dimension(imax) :: xt,noise,signal,adap,pulse
real, parameter :: dt=1.e-2,tau=300.


!!!!!!!!!!!!!!!!!!!펄스 신호 생성!!!!!!!!!!!!!!!!!!!!!!!
t=0.
do i=1,imax
t=i*dt
!펄스 높이 늘리기
if ((t .gt. 1.) .and. (t .lt. 10.)) then
pulse(i) = 0.001

elseif ((t .gt. 520.) .and. (t .lt. 530.)) then
pulse(i) = 0.01

elseif ((t .gt. 1040.) .and. (t .lt. 1050.)) then
pulse(i) = 0.1

elseif ((t .gt. 1560.) .and. (t .lt. 1570.)) then
pulse(i) = 1.0

elseif ((t .gt. 2080.) .and. (t .lt. 2090.)) then
pulse(i) = 10.

!펄스 폭 늘리기
elseif ((t .gt. 3120.) .and. (t .lt. 3120+dt)) then
pulse(i) = 1.

elseif ((t .gt. 3620.) .and. (t .lt. 3620+dt*10.)) then
pulse(i) = 1.

elseif ((t .gt. 4120.) .and. (t .lt. 4120+dt*100.)) then
pulse(i) = 1.

elseif ((t .gt. 4620.) .and. (t .lt. 4620+dt*1000.)) then
pulse(i) = 1.

else
pulse(i) = 0.
end if
enddo
!!!!!!!!!!!!!!!!!!!펄스신호 생성 끝!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!노이즈 생성 !!!!!!!!!!!!!!!!!!
call random_seed()  
do i = 1, imax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do
!!!!!!!!!!!!!!!!!!!!!노이즈  생성 끝!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!운동방정식 오일러 계산!!!!!!!!!!!!!!!
open(1,file='xt')
open(3,file='sig2')

!!!!!!!!초기조건!!!!!!!!!!!!
t=0.;x=1.;v=1.
xbar=1. 
!!!!!!!!!!!!!초기조건 끝!!!!!!!!!!!!

do i=1,imax
!!!!!!!!!!!!!!!!!!신호 고르기!!!!!!!!!!!!!!!!!!!!
signal(i)=pulse(i)
!signal(i)=(f/m)*sin(0.5*t+3.14)
!!!!!!!!!!!!!!!!!신호 고르기 끝!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!adaptation force 생성!!!!!!!!!!!!!!!
dxbar=(x-xbar)/tau
xbar=xbar+dxbar*dt

adap(i)=sign(1.,x-xbar)
!!!!!!!!!!!!!!!!adaptation force 생성 끝 !!!!!!!!!!!!!!!!!!

! 자발 진동 상황
!dv=-0.1*v-x+adap(i)+noise(i)

! 자발 진동 상황 노이즈 제거
!dv=-0.1*v-x+adap(i)+0.*noise(i)


! 일반적 상황
dv=-0.1*v-x+adap(i)+signal(i)+noise(i)

dx=v
v=v+dv*dt
x=x+dx*dt
t=i*dt
write(1,*) t,x
write(3,*) t,noise(i)+pulse(i)
xt(i)=x
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!운동방정식 계산 끝!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!노이즈 vs 펄스 신호 비교!!!!!!!!!!!!!!!!!
!t=0.
!write(3,*) ''
!do i=1,imax
!t=i*dt
!write(3,*) t,pulse(i)
!enddo
!!!!!!!!!!!!!!!!!!노이즈 vs 펄스 신호 비교 끝 !!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!xt vs 펄스 신호 비교!!!!!!!!!!!!!!!!!
t=0.
write(1,*) ''
do i=1,imax
t=i*dt
write(1,*) t,pulse(i)
enddo
!!!!!!!!!!!!!!!!!!xt vs 펄스 신호 비교 끝 !!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DFT!!!!!!!!!!!!!!!!!!!!!!!!!!
open(2,file='wx')
do j=1,50
om=6.28*j/(imax*dt)
!reset
t=0.
sumr=0.
sumi=0.
!sum
do i=1,imax
t=i*dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo
write(2,*) om, sqrt(sumi**2+sumr**2)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!DFT 끝!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program
