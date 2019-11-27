program toy2d
implicit none

integer :: NN,Mx,My,dN,columns,rows
real*8 :: dt,dx,dy,Lx,Ly,T
integer :: i,j
character*4 :: arg1
integer*4 :: p

call get_command_argument(1, arg1)
read(arg1, '(i4)') p

columns=15
rows=p

Lx=150*(columns+1);
Ly=150*(rows+1);
T=200;
NN=1000000;
dN=5000;
Mx=Lx/10;
My=Ly/10;
dt=T/NN;
dx=Lx/Mx;
dy=Ly/My;

call pdesolver(NN,dN,Mx,My,dt,dx,dy,Lx,Ly,T,columns,rows,p)

end program toy2d
!-----------------------------------------------
subroutine pdesolver(NN,dN,Mx,My,dt,dx,dy,Lx,Ly,T,columns,rows,p)

integer :: NN,Mx,My,dN,rows,columns
real*8 :: dt,dx,dy,Lx,Ly,T 
integer*4 :: p
integer :: i,j,k,n,m
real*8 :: D,pi,sigma,time,F,per,eps,distx,disty
real*8, dimension(1:Mx+1,1:My+1,1:dN+1) :: c
real*8 :: dxV(Mx+1,My+1),dyV(Mx+1,My+1),ddxV(Mx+1,My+1),ddyV(Mx+1,My+1) 
real*8 :: rowx(columns,rows),rowy(columns,rows),rowxx(columns,rows),rowyy(columns,rows),r(columns,rows,2),noise(columns,rows,2)
integer, allocatable :: seed(:)
integer*4 :: q(12)
character*48 :: filename
integer*2 :: z(110,2)

do j=1,12
	q(j)=j**8+5788994+10000000*2
end do
call random_seed(size = n)
call random_seed(put=q(1:n))
call random_number(r)

do i=1,columns
	do j=1,rows
		noise(i,j,1)=-15+30*r(i,j,1)
		noise(i,j,2)=-15+30*r(i,j,2)
	end do
end do

D=600
pi=3.1415927;
sigma=60
time=0
eps=40000
per=40
distx=150
disty=150

do i=1,Mx+1
	do j=1,My+1
		c(i,j,1)=1
	end do
end do

do i=1,Mx+1
	do j=1,My+1
		do k=1,columns
			do m=1,rows
				rowx(k,m)=1/(sigma**3*sqrt(2*pi))*((i-1)*dx-(distx+150*(k-1))+noise(k,m,1)) &
					*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2))
				rowy(k,m)=1/(sigma**3*sqrt(2*pi))*((j-1)*dy-(disty+150*(m-1))) &
					*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2))
				rowxx(k,m)=1/(sigma**3*sqrt(2*pi))*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2 &
					+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2)) &
					-1/(sigma**5*sqrt(2*pi))*((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2 &
					*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2))
				rowyy(k,m)=1/(sigma**3*sqrt(2*pi))*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2 &
					+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2)) & 
					-1/(sigma**5*sqrt(2*pi))*((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2 &
					*exp(-(((i-1)*dx-(distx+150*(k-1))+noise(k,m,1))**2+((j-1)*dy-(disty+150*(m-1))+noise(k,m,2))**2)/(2*sigma**2))
			end do
		end do
		dxV(i,j)=sum(rowx)
		dyV(i,j)=sum(rowy)
		ddxV(i,j)=sum(rowxx)
		ddyV(i,j)=sum(rowyy)
	end do
end do

write(filename,'("output/output2d_toymodel_1_sim013",I4,".txt")') p
open(unit=p, file=filename)
write(p,*) c(1:Mx+1,1:My+1,1)

do k=1,9*NN/(10*dN)
	do n=1,dN
		if (mod(time,per)<0.7*per) then
			F=eps
		else
			F=0
		endif
		do i=2,Mx
			do j=2,My
				c(i,j,n+1)=c(i,j,n)+dt*(D*((c(i+1,j,n)-2*c(i,j,n)+c(i-1,j,n))/dx**2 &
					+(c(i,j+1,n)-2*c(i,j,n)+c(i,j-1,n))/dy**2) &
					+F*(dxV(i,j)*(c(i+1,j,n)-c(i-1,j,n))/(2*dx) + dyV(i,j)*(c(i,j+1,n)-c(i,j-1,n))/(2*dy) &
					+c(i,j,n)*(ddxV(i,j)+ddyV(i,j))))
			end do
		end do
		do j=2,My
			c(1,j,n+1)=(D/(D-dx*F*dxV(1,j)))*c(2,j,n+1)
			c(Mx+1,j,n+1)=(D/(D+dx*F*dxV(Mx+1,j)))*c(Mx,j,n+1)
		end do
		do i=1,Mx+1
			c(i,1,n+1)=(D/(D-dy*F*dyV(i,1)))*c(i,2,n+1)
			c(i,My+1,n+1)=(D/(D+dy*F*dyV(i,My+1)))*c(i,My,n+1)
		end do
		time=time+dt
	end do
	do i=1,Mx+1
		do j=1,My+1
			c(i,j,1)=c(i,j,dN+1)
		end do
	end do
	write(p,*) c(1:Mx+1,1:My+1,1)
end do
do k=1,NN/(10*dN)
	do n=1,dN
		if (mod(time,per)<0.7*per) then
			F=eps
		else
			F=0
		endif
		do i=2,Mx
			do j=2,My
				c(i,j,n+1)=c(i,j,n)+dt*(D*((c(i+1,j,n)-2*c(i,j,n)+c(i-1,j,n))/dx**2 &
					+(c(i,j+1,n)-2*c(i,j,n)+c(i,j-1,n))/dy**2) &
					+F*(dxV(i,j)*(c(i+1,j,n)-c(i-1,j,n))/(2*dx) + dyV(i,j)*(c(i,j+1,n)-c(i,j-1,n))/(2*dy) &
					+c(i,j,n)*(ddxV(i,j)+ddyV(i,j))))
			end do
		end do
		do j=2,My
			c(1,j,n+1)=(D/(D-dx*F*dxV(1,j)))*c(2,j,n+1)
			c(Mx+1,j,n+1)=(D/(D+dx*F*dxV(Mx+1,j)))*c(Mx,j,n+1)
		end do
		do i=1,Mx+1
			c(i,1,n+1)=(D/(D-dy*F*dyV(i,1)))*c(i,2,n+1)
			c(i,My+1,n+1)=(D/(D+dy*F*dyV(i,My+1)))*c(i,My,n+1)
		end do
		time=time+dt
	end do
	do i=1,Mx+1
		do j=1,My+1
			c(i,j,1)=c(i,j,dN+1)
		end do
	end do
	write(p,*) c(1:Mx+1,1:My+1,1)
end do


close(p)

end subroutine pdesolver
