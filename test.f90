program aaaa
implicit none
integer i
real t,dt
real, dimension(imax) :: pulse
open(1,file='aa')
t=0.
dt=1.e-2
do i=1,1000
t=i*dt

    if ((t .gt. 1.) .and. (t .lt. 2.)) then
      pulse(i) = 0.01
write(1,*) t, x(i)
    else
      pulse(i) = 0.
write(1,*) t, x(i)
    end if
enddo
end
