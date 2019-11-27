program toy1d
implicit none

integer :: NN,M,dN,NC
real*8 :: dt,dx,L,T,dist
integer :: i,j
integer*4 :: p,nr 

nr=7;
do p=1,nr 

L=2400;
T=200;
NN=1000000;
dN=1000;
M=600;
NC=1;
dist=L/(NC+1)
dt=T/NN;
dx=L/M;

call pdesolver(NN,dN,M,dt,dx,L,T,NC,dist,p)

end do
end program toy1d
!-----------------------------------------------
subroutine pdesolver(NN,dN,M,dt,dx,L,T,NC,dist,p)

integer :: NN,M,dN,NC
real*8 :: dt,dx,L,T,dNC,dist
integer*4 :: p,nr
integer :: i,j,k,n
real*8 :: D,pi,sigma,time,F,per,beta,eps,alpha
real*8, dimension(1:M+1,1:dN+1) :: c
real*8 :: dV(M+1),ddV(M+1) 
real*8 :: r(NC),noise(NC),row1(NC),row2(NC)
integer, allocatable :: seed(:)
integer*4 :: q(12)
character*48 :: filename
integer*2 :: z(121,2)

D=600
pi=3.1415927;
sigma=60
time=0
eps=40000
per=40
alpha=0.7
beta=1

!dNC=(L-2*dist)/(NC-1);

do j=1,12
	q(j)=j**8+5788994+10000000*2
end do
call random_seed(size = n)
call random_seed(put=q(1:n))
call random_number(r)

do i=1,NC
	noise(i)=0 !-(-5+5*p)+2*(-5+5*p)*r(i)
end do

do i=1,11
	z(1+11*(i-1):11*i,1)=i
	do j=1,11
		z(11*(i-1)+j,2)=j
	end do
end do

do i=1,M+1
	c(i,1)=1
end do

do j=1,M+1
	do k=1,NC
		row1(k)=1/(sigma**3*sqrt(2*pi))*((j-1)*dx-(dist+dNC*(k-1))+noise(k))*exp(-((j-1)*dx &
			-(dist+dNC*(k-1))+noise(k))**2/(2*sigma**2))
		row2(k)=1/(sigma**3*sqrt(2*pi))*exp(-((j-1)*dx-(dist+dNC*(k-1))+noise(k))**2/(2*sigma**2)) &
			-1/(sigma**5*sqrt(2*pi))*((j-1)*dx-(dist+dNC*(k-1))+noise(k))**2*exp(-((j-1)*dx &
			-(dist+dNC*(k-1))+noise(k))**2/(2*sigma**2))
	end do
	dV(j)=sum(row1)
	ddV(j)=sum(row2)
end do

write(filename,'("output/output1d_toymodel_1_sim034",I4,".txt")') p
open(unit=p, file=filename)
write(p,*) c(1:M+1,1)

do k=1,9*NN/(10*dN)
	do n=1,dN
		if (mod(time,per)<alpha*per) then
			F=eps
		else
			F=0
		endif
		do i=2,M
			c(i,n+1)=c(i,n)+dt*(D*((c(i+1,n)-2*c(i,n)+c(i-1,n))/dx**2) &
				+F*(dV(i)*(c(i+1,n)-c(i-1,n))/(2*dx) + ddV(i)*c(i,n)))
		end do
		c(1,n+1)=(D/(D-dx*F*dV(1)))*c(2,n+1)
		c(M+1,n+1)=(D/(D+dx*F*dV(M+1)))*c(M,n+1)
		time=time+dt
	end do
	do i=1,M+1
		c(i,1)=c(i,dN+1)
	end do
	write(p,*) c(1:M+1,1)
end do
do k=1,NN/(10*dN)
	do n=1,dN
		if (mod(time,per)<alpha*per) then
			F=eps
		else
			F=0
		endif
		do i=2,M
			c(i,n+1)=c(i,n)+dt*(D*((c(i+1,n)-2*c(i,n)+c(i-1,n))/dx**2) &
				+F*(dV(i)*(c(i+1,n)-c(i-1,n))/(2*dx) + ddV(i)*c(i,n)))
		end do
		c(1,n+1)=(D/(D-dx*F*dV(1)))*c(2,n+1)
		c(M+1,n+1)=(D/(D+dx*F*dV(M+1)))*c(M,n+1)
		time=time+dt
	end do
	do i=1,M+1
		c(i,1)=c(i,dN+1)
	end do
	write(p,*) c(1:M+1,1)
end do

close(p)

end subroutine pdesolver