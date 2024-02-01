program a
implicit none
integer i
real h,t,x,v,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,gam,om0

open(1,file='xt')
!test
open(2,file='aa')

t=0.;x=0.;v=1.;h=1.e-2
write(1,*) t,x

!test
write(2,*) t,x

do i=0,10000
k1x=fx(t,x,v);k1v=fv(t,x,v)
k2x=fx(t+h/2.,x+h*k1x/2.,v+h*k1v/2.);k2v=fv(t+h/2.,x+h*k1x/2.,v+h*k1v/2.)
k3x=fx(t+h/2.,x+h*k2x/2.,v+h*k2v/2.);k3v=fv(t+h/2.,x+h*k2x/2.,v+h*k2v/2.)
k4x=fx(t+h,x+h*k3x,v+h*k3v);k4v=fv(t+h,x+h*k3x,v+h*k3v)
t=t+h;x=x+h*(k1x+2.*k2x+2.*k3x+k4x)/6.;v=v+h*(k1v+2.*k2v+2.*k3v+k4v)/6.
write(1,*) t,x

!test
gam=0.1;om0=2.
write(2,*) t,-(1./2.)*exp(-gam*t)*cos(om0*t+acos(0.))
end do


contains

real function fx(t,x,v)
real t,x,v
fx=v
end function fx


real function fv(t,x,v)
real t,x,v,gam,om0
gam=0.1;om0=2.
fv=-2.*gam*v-om0**2*x
end function fv
end program 
