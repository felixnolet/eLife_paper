program cell1d_hpc
implicit none

integer :: NN,M
real*8 :: dt,dx,L,T
integer :: i,j
integer*4 :: p
integer*2 :: z(441,2)
character*4 :: arg

call get_command_argument(1, arg)
read(arg, '(i4)') p

do i=1,21
	z(1+21*(i-1):21*i,1)=i
	do j=1,21
		z(21*(i-1)+j,2)=j
	end do
end do

L=1800;
T=2000;
NN=2000000;
M=600;
dt=T/NN;
dx=L/M;

call pdesolver(NN,M,dt,dx,L,T,p)

end program cell1d_hpc
!-----------------------------------------------
subroutine pdesolver(NN,M,dt,dx,L,T,p)

integer :: NN,M
real*8 :: dt,dx,L,T 
integer*4 :: p,nr
integer :: i,j,k,n
real*8 :: a2,a3,b2,b3,e1,e2,e3,kk,D,sigma,pi,n1,n2,n3,dba,dbb
real*8, dimension(1:M+1,1:201) :: u,v
real*8 :: a1(M+1),b1(M+1),row(M+1),r(M+1),s(M+1)
real*8 :: nuclei(15)
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

pi=3.1415927;
sigma=50

a2=0.4
a3=0.01
b2=2
b3=0.06
e1=0.35
e2=0.30
e3=0.32
D=0
kk=0.015
n1=11
n2=3.5
n3=17

do i=1,M+1
	a1(i)=0.8
	b1(i)=4
	do k=1,15
	nuclei(k)=1/(sigma*sqrt(2*pi))*exp(-(dx*(i-1)-150*k)**2/(2*sigma**2))
	end do
	row(i)=sum(nuclei)+100/(sigma*sqrt(2*pi))*exp(-(dx*(i-1)-1200)**2/(2*sigma**2))
end do 

!dba=(1+0.01*p)*0.8
!dbb=(1+0.01*p)*4
do i=1,61
	!a1(i)=dba !-0.862+0.036/20*(i-1) !-0.851+0.001*p 
	!b1(i)=dbb
end do
do i=62,282
	!a1(i)=dba-((dba-0.8)/222)*(i-61) !-0.826-0.024/80*(i-21)
	!b1(i)=dbb-((dbb-4)/222)*(i-61)
end do
do i=283,319
	!a1(i)=1.05*0.8 !0.795+0.005*p
	!b1(i)=1.05*4
end do
do i=320,540
	!a1(i)=0.8+((dba-0.8)/222)*(i-319) !-0.85+0.024/80*(i-501)
	!b1(i)=4+((dbb-4)/222)*(i-319)
end do
do i=541,601
	!a1(i)=dba !-0.826-0.036/20*(i-581) !-0.851+0.001*p
	!b1(i)=dbb
end do

do i=1,M+1
	a1(i)=0.4904
	b1(i)=2.452 !a1(i)+0.1*row(i)
end do
do i=126,146
	a1(i)=0.6
	b1(i)=3.0
end do
do i=26,46
	!a1(i)=0.4904+0.1096*r(1)
	!a1(i+50)=0.4904+0.1096*r(2)
	!a1(i+150)=0.4904+0.1096*r(3)
	!a1(i+200)=0.4904+0.1096*r(4)
	!a1(i+250)=0.4904+0.1096*r(5)
	!a1(i+300)=0.4904+0.1096*r(6)
	!a1(i+350)=0.4904+0.1096*r(7)
	!a1(i+400)=0.4904+0.1096*r(8)
	!a1(i+450)=0.4904+0.1096*r(9)
	!a1(i+500)=0.4904+0.1096*r(10)
	!b1(i)=2.452+0.548*r(1)
	!b1(i+50)=2.452+0.548*r(2)
	!b1(i+150)=2.452+0.548*r(3)
	!!b1(i+200)=2.452+0.548*r(4)
	!b1(i+250)=2.452+0.548*r(5)
	!b1(i+300)=2.452+0.548*r(6)
	!b1(i+350)=2.452+0.548*r(7)
	!b1(i+400)=2.452+0.548*r(8)
	!b1(i+450)=2.452+0.548*r(9)
	!b1(i+500)=2.452+0.548*r(10)
end do

dba=0.8
dbb=4
do i=1,21
	!a1(i)=(0.986+(0.042/20)*(i-1))*dba
	!b1(i)=(0.986+(0.042/20)*(i-1))*dbb
end do
do i=22,100
	!a1(i)=(1.028-(0.028/80)*(i-21))*dba
	!b1(i)=(1.028-(0.028/80)*(i-21))*dbb
end do
do i=101,290
	!a1(i)=dba
	!b1(i)=dbb
end do
do i=291,311
	!a1(i)=(1.025+0.0003*(p-1))*dba
	!b1(i)=(1.025+0.0003*(p-1))*dbb
end do
do i=312,501
	!a1(i)=dba
	!b1(i)=dbb
end do
do i=502,580
	!a1(i)=(1+(0.028/80)*(i-501))*dba
	!b1(i)=(1+(0.028/80)*(i-501))*dbb
end do
do i=581,601
	!a1(i)=(1.028-(0.042/20)*(i-581))*dba
	!b1(i)=(1.028-(0.042/20)*(i-581))*dbb
end do

do i=1,M+1
	u(i,1)=0 !-0.9 !-0.05+0.1*s(i)
	v(i,1)=0 !-0.2
end do

write(filename,'("output/output1d_cell_2_sim030",I4,".txt")') p
open(unit=p, file=filename)
write(p,*) a1(1:M+1)
write(p,*) u(1:M+1,1)

do k=1,2*NN/2000
	do n=1,200
		do i=2,M
			u(i,n+1)=u(i,n)+dt*(D*((u(i+1,n)-2*u(i,n)+u(i-1,n))/(dx**2)) &
				+(a1(i)+b1(i)*(u(i,n)**n1)/(e1**n1+u(i,n)**n1))*(v(i,n)-u(i,n)) &
				-(a2+b2*(e2**n2)/(e2**n2+u(i,n)**n2))*u(i,n) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*u(i,n)+kk)
			v(i,n+1)=v(i,n)+dt*(D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/(dx**2)) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*v(i,n)+kk)
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,201)
		v(i,1)=v(i,201)
	end do
	!write(p,*) u(1:M+1,1)
end do
do k=1,6*NN/2000
	do n=1,200
		do i=2,M
			u(i,n+1)=u(i,n)+dt*(D*((u(i+1,n)-2*u(i,n)+u(i-1,n))/(dx**2)) &
				+(a1(i)+b1(i)*(u(i,n)**n1)/(e1**n1+u(i,n)**n1))*(v(i,n)-u(i,n)) &
				-(a2+b2*(e2**n2)/(e2**n2+u(i,n)**n2))*u(i,n) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*u(i,n)+kk)
			v(i,n+1)=v(i,n)+dt*(D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/(dx**2)) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*v(i,n)+kk)
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,201)
		v(i,1)=v(i,201)
	end do
	!write(p,*) u(1:M+1,1)
end do
do k=1,2*NN/2000
	do n=1,200
		do i=2,M
			u(i,n+1)=u(i,n)+dt*(D*((u(i+1,n)-2*u(i,n)+u(i-1,n))/(dx**2)) &
				+(a1(i)+b1(i)*(u(i,n)**n1)/(e1**n1+u(i,n)**n1))*(v(i,n)-u(i,n)) &
				-(a2+b2*(e2**n2)/(e2**n2+u(i,n)**n2))*u(i,n) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*u(i,n)+kk)
			v(i,n+1)=v(i,n)+dt*(D*((v(i+1,n)-2*v(i,n)+v(i-1,n))/(dx**2)) &
				-(a3+b3*(u(i,n)**n3)/(e3**n3+u(i,n)**n3))*v(i,n)+kk)
		end do
		u(1,n+1)=u(2,n+1)
		v(1,n+1)=v(2,n+1)
		u(M+1,n+1)=u(M,n+1)
		v(M+1,n+1)=v(M,n+1)
	end do
	do i=1,M+1
		u(i,1)=u(i,201)
		v(i,1)=v(i,201)
	end do
	write(p,*) u(1:M+1,1)
	!write(p,*) v(1:M+1,1)
end do

close(p)

end subroutine pdesolver