program aa
integer i
real, dimension(1000) :: a
open(1,file='ab')
do i=1,1000
a(i)=rand()
write(1,*) i,a(i)
enddo
end
