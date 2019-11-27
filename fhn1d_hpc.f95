program fhn1d2
implicit none

integer :: NN,M
real*8 :: dt,dx,L,T
integer :: i,j
integer :: p
character*4 :: arg

call get_command_argument(1, arg)
read(arg, '(i4)') p

L=1800;
T=2000;
NN=2000000;
M=600;
dt=T/NN;
dx=L/M;

call pdesolver(NN,M,dt,dx,L,T,p)

end program fhn1d2
!-----------------------------------------------
subroutine pdesolver(NN,M,dt,dx,L,T,p)

integer :: NN,M
real*8 :: dt,dx,L,T,timescale
integer :: p,nr
integer :: i,j,k,n
real*8 :: D,h,diff,const,db
real*8, dimension(1:M+1,1:401) :: u,v
real*8 :: r(M+1),b(M+1),s(M+1),c(M+1),a(M+1),eps(M+1)
integer, allocatable :: seed(:)
integer*4 :: q(12)
character*48 :: filename
integer*2 :: z(441,2)

do i=1,21
	z(1+21*(i-1):21*i,1)=i
	do j=1,21
		z(21*(i-1)+j,2)=j
	end do
end do

do j=1,12
	q(j)=j**8+5788994+10000000*p
end do
call random_seed(size = n)
call random_seed(put=q(1:n))
call random_number(r)
call random_number(s)

eps=0.01
D=1
h=0.5

diff=0
const=0.0005*(p-1)

timescale=3

do i=1,M+1
	!a(i)=-0.85
	b(i)=0.05 !4.3+0.2*z(p,1) !-0.1+0.1*p !+0.2*r(i)
	c(i)=1.2 !4.3+0.2*z(p,1) !-0.1+0.1*p !+0.2*r(i)
end do

db=-0.851+0.001*p
do i=1,61
	!a(i)=db !-0.862+0.036/20*(i-1) !-0.851+0.001*p 
end do
do i=62,282
	!a(i)=db-((db+0.85)/222)*(i-61) !-0.826-0.024/80*(i-21)
end do
do i=283,319
	!a(i)=-0.75 !0.795+0.005*p
end do
do i=320,540
	!a(i)=-0.85+((db+0.85)/222)*(i-319) !-0.85+0.024/80*(i-501)
end do
do i=541,601
	!a(i)=db !-0.826-0.036/20*(i-581) !-0.851+0.001*p
end do

do i=1,290
	!a(i)=-0.85+const*(1-(i-1)/290.)
end do
do i=291,311
	!a(i)=-0.8 !0.795+0.005*p
end do
do i=312,601
	!a(i)=-0.85+const*((i-311)/290.)
end do

do i=1,601
	a(i)=-0.882!-0.772
	!eps(i)=0.03
end do
do i=126,146
	a(i)=-0.59
end do
do i=26,46
	!a(i)=-0.882+0.292*r(1)
	!a(i+50)=-0.882+0.292*r(2)
	!a(i+150)=-0.882+0.292*r(3)
	!a(i+200)=-0.882+0.292*r(4)
	!a(i+250)=-0.882+0.292*r(5)
	!a(i+300)=-0.882+0.292*r(6)
	!a(i+350)=-0.882+0.292*r(7)
	!a(i+400)=-0.882+0.292*r(8)
	!a(i+450)=-0.882+0.292*r(9)
	!a(i+500)=-0.882+0.292*r(10)
end do

do i=1,21
	!a(i)=-0.832+(0.041/20)*(i-1)
end do
do i=22,100
	!a(i)=-0.791-(0.028/80)*(i-21)
end do
do i=101,290
	!a(i)=-0.819
end do
do i=291,311
	!a(i)=-0.78+0.0005*(p-1)
end do
do i=312,501
	!a(i)=-0.819
end do
do i=502,580
	!a(i)=-0.819+(0.028/80)*(i-501)
end do
do i=581,601
	!a(i)=-0.791-(0.041/20)*(i-581)
end do

do i=1,M+1
	u(i,1)=0 !-0.05+0.1*s(i)
	v(i,1)=0
end do

u(1,1)=u(2,1)
v(1,1)=v(2,1)
u(M+1,1)=u(M,1)
v(M+1,1)=v(M,1)

write(filename,'("output/output1d_2_sim053",I4,".txt")') p
open(unit=p, file=filename)
write(p,*) a(1:M+1)
write(p,*) u(1:M+1,1)

do k=1,40*NN/40000
	do n=1,400
		do i=2,M
			u(i,n+1)=u(i,n)+dt*timescale*(diff/timescale*((u(i+1,n)-2*u(i,n)+u(i-1,n))/dx**2)-(u(i,n))**3+c(i)*(u(i,n))**2+h*u(i,n)-v(i,n))
			v(i,n+1)=v(i,n)+dt*timescale*(diff/timescale*D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/dx**2)+eps(i)*(u(i,n)-b(i)*v(i,n)+a(i)))
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,401)
		v(i,1)=v(i,401)
	end do
	!write(p,*) u(1:M+1,1)
end do
do k=1,40*NN/40000
	do n=1,400
		do i=2,M
			u(i,n+1)=u(i,n)+dt*timescale*(diff/timescale*((u(i+1,n)-2*u(i,n)+u(i-1,n))/dx**2)-(u(i,n))**3+c(i)*(u(i,n))**2+h*u(i,n)-v(i,n))
			v(i,n+1)=v(i,n)+dt*timescale*(diff/timescale*D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/dx**2)+eps(i)*(u(i,n)-b(i)*v(i,n)+a(i)))
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,401)
		v(i,1)=v(i,401)
	end do
	!write(p,*) u(1:M+1,1)
end do
do k=1,20*NN/40000
	do n=1,400
		do i=2,M
			u(i,n+1)=u(i,n)+dt*timescale*(diff/timescale*((u(i+1,n)-2*u(i,n)+u(i-1,n))/dx**2)-(u(i,n))**3+c(i)*(u(i,n))**2+h*u(i,n)-v(i,n))
			v(i,n+1)=v(i,n)+dt*timescale*(diff/timescale*D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/dx**2)+eps(i)*(u(i,n)-b(i)*v(i,n)+a(i)))
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,401)
		v(i,1)=v(i,401)
	end do
	write(p,*) u(1:M+1,1)
end do

close(p)

end subroutine pdesolver