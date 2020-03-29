program toy1d
implicit none

integer :: NN,M,dN
real*8 :: dt,dx,L,T,dist,dNC,eps,sigma,D,beta,gamma
integer :: i,j
character*8 :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9
integer*4 :: id,NC

call get_command_argument(1, arg1)
read(arg1, '(i4)') id

call get_command_argument(2, arg2)
read(arg2, '(i4)') NC

call get_command_argument(3, arg3)
read(arg3,*) dist

call get_command_argument(4, arg4)
read(arg4,*) dNC

call get_command_argument(5, arg5)
read(arg5,*) eps

call get_command_argument(6, arg6)
read(arg6,*) sigma

call get_command_argument(7, arg7)
read(arg7,*) D

call get_command_argument(8, arg8)
read(arg8,*) beta

call get_command_argument(9, arg9)
read(arg9,*) gamma

!definition of time discretization
T=4000;
NN=400000;
dN=40;
dt=T/NN;

!definition of space discretization, L defined indirectly
L=2*dist+(NC-1)*dNC;
M=nint(L/5);
dx=L/M;

call pdesolver(NN,dN,M,dt,dx,L,T,id,NC,dist,dNC,eps,sigma,D,beta,gamma)

end program toy1d
!-----------------------------------------------
subroutine pdesolver(NN,dN,M,dt,dx,L,T,id,NC,dist,dNC,eps,sigma,D,beta,gamma)

integer :: NN,M,dN,mid
real*8 :: dt,dx,L,T,dNC,dist,beta,gamma
integer*4 :: id,NC
integer :: i,j,k,n,s
real*8 :: D,pi,sigma,per,eps,alpha,cn,Tn,time,v1,v2
real*8, dimension(1:M+1,1:dN+1) :: c
real*8 :: dV(M+1),ddV(M+1) 
real*8 :: r(NC),noise(NC),row1(NC),row2(NC),pos(NC),F(NC),bf(NC)
real*8 :: timeSold(NC),timeMold(NC),timeS(NC),timeM(NC)
real*8 :: xt(2)
integer, allocatable :: seed(:)
integer*4 :: q(12)
character*64 :: filename,filename1,filename2
integer :: unit2

pi=3.1415927;
alpha=0.7
per=40

do j=1,12
	q(j)=j**8+5788994+10000000*2
end do
call random_seed(size = n)
call random_seed(put=q(1:n))
call random_number(r)

do i=1,M+1
	c(i,1)=1
end do

do k=1,NC
	timeS(k)=0
	timeM(k)=0
	noise(k)=0 !-40+80*r(k)
	bf(k)=1
end do

mid=(NC+1)/2
bf(mid)=beta
bf(mid-1)=beta
bf(mid+1)=beta

do k=1,NC
	pos(k)=dist+dNC*(k-1) !+noise(k)
end do

v1=0.01
v2=0.01

write(filename1,'("output/output1d_toymodel_1_T_sim006_id",i4,"_params.txt")') id
open(unit=0, file=filename1)
write(0,*) NC
write(0,*) pos(1:NC)
write(0,*) dist
write(0,*) dNC
write(0,*) eps
write(0,*) sigma
write(0,*) D
write(0,*) alpha
write(0,*) beta
write(0,*) gamma
write(0,*) v1
write(0,*) v2
write(0,*) L
write(0,*) M
write(0,*) T
write(0,*) NN
close(0)

unit2=10000+id;

write(filename,'("output/output1d_toymodel_1_T_sim006_id",i4,".txt")') id
write(filename2,'("output/output1d_toymodel_1_T_sim006_id",i4,"_xt.txt")') id
open(unit=id, file=filename)
open(unit=unit2, file=filename2)
!write(id,*) c(1:M+1,1)

time=0
do s=1,NN/(100*dN)*95
	do n=1,dN
		do k=1,NC
			timeSold(k)=timeS(k)
			timeMold(k)=timeM(k)
			cn=c(nint(pos(k)/dx)+1,n) !concentration at nucleus k
			Tn=per/(1+gamma*(cn-1.1)) !local period of nucleus k
			if (timeS(k)<alpha*Tn .and. timeM(k)==0) then
				if (timeS(k)<v1*per) then
					F(k)=timeS(k)/(v1*per)
				else
					F(k)=1
				endif
				timeS(k)=timeS(k)+dt
				timeM(k)=0
			else if (timeS(k)>=alpha*Tn .and. timeM(k)==0) then
				F(k)=1
				timeS(k)=0
				timeM(k)=timeM(k)+dt
			else if (timeS(k)==0 .and. timeM(k)<(1-alpha)*Tn) then
				if (timeM(k)<v2*per) then
					F(k)=1-timeM(k)/(v2*per)
				else
					F(k)=0
				endif
				timeS(k)=0
				timeM(k)=timeM(k)+dt
			else if (timeS(k)==0 .and. timeM(k)>=(1-alpha)*Tn) then
				F(k)=0
				timeS(k)=timeS(k)+dt
				timeM(k)=0
			endif
		end do
		do k=1,NC
			if (timeS(k)==0 .and. timeSold(k)>0) then
				xt(1)=pos(k)
				xt(2)=time
				write(unit2,*) xt
			endif
		end do
		do j=1,M+1
			do k=1,NC
				row1(k)=bf(k)*F(k)/(sigma**3*sqrt(2*pi))*((j-1)*dx-pos(k))*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2))
				row2(k)=bf(k)*F(k)/(sigma**3*sqrt(2*pi))*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2)) &
					-bf(k)*F(k)/(sigma**5*sqrt(2*pi))*((j-1)*dx-pos(k))**2*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2))
			end do
			dV(j)=eps*sum(row1)
			ddV(j)=eps*sum(row2)
		end do
		do i=2,M
			c(i,n+1)=c(i,n)+dt*(D*((c(i+1,n)-2*c(i,n)+c(i-1,n))/dx**2) &
				+ (dV(i)*(c(i+1,n)-c(i-1,n))/(2*dx) + ddV(i)*c(i,n)))
		end do
		c(1,n+1)=(D/(D-dx*dV(1)))*c(2,n+1)
		c(M+1,n+1)=(D/(D+dx*dV(M+1)))*c(M,n+1)
		time=time+dt
	end do
	do i=1,M+1
		c(i,1)=c(i,dN+1)
	end do
	!write(id,*) c(1:M+1,1)
end do
do s=1,NN/(100*dN)*5
	do n=1,dN
		do k=1,NC
			timeSold(k)=timeS(k)
			timeMold(k)=timeM(k)
			cn=c(nint(pos(k)/dx)+1,n) !concentration at nucleus k
			Tn=per/(1+gamma*(cn-1.1)) !local period of nucleus k
			if (timeS(k)<alpha*Tn .and. timeM(k)==0) then
				if (timeS(k)<v1*per) then
					F(k)=timeS(k)/(v1*per)
				else
					F(k)=1
				endif
				timeS(k)=timeS(k)+dt
				timeM(k)=0
			else if (timeS(k)>=alpha*Tn .and. timeM(k)==0) then
				F(k)=1
				timeS(k)=0
				timeM(k)=timeM(k)+dt
			else if (timeS(k)==0 .and. timeM(k)<(1-alpha)*Tn) then
				if (timeM(k)<v2*per) then
					F(k)=1-timeM(k)/(v2*per)
				else
					F(k)=0
				endif
				timeS(k)=0
				timeM(k)=timeM(k)+dt
			else if (timeS(k)==0 .and. timeM(k)>=(1-alpha)*Tn) then
				F(k)=0
				timeS(k)=timeS(k)+dt
				timeM(k)=0
			endif
		end do
		do k=1,NC
			if (timeS(k)==0 .and. timeSold(k)>0) then
				xt(1)=pos(k)
				xt(2)=time
				write(unit2,*) xt
			endif
		end do
		do j=1,M+1
			do k=1,NC
				row1(k)=bf(k)*F(k)/(sigma**3*sqrt(2*pi))*((j-1)*dx-pos(k))*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2))
				row2(k)=bf(k)*F(k)/(sigma**3*sqrt(2*pi))*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2)) &
					-bf(k)*F(k)/(sigma**5*sqrt(2*pi))*((j-1)*dx-pos(k))**2*exp(-((j-1)*dx-pos(k))**2/(2*sigma**2))
			end do
			dV(j)=eps*sum(row1)
			ddV(j)=eps*sum(row2)
		end do
		do i=2,M
			c(i,n+1)=c(i,n)+dt*(D*((c(i+1,n)-2*c(i,n)+c(i-1,n))/dx**2) &
				+ (dV(i)*(c(i+1,n)-c(i-1,n))/(2*dx) + ddV(i)*c(i,n)))
		end do
		c(1,n+1)=(D/(D-dx*dV(1)))*c(2,n+1)
		c(M+1,n+1)=(D/(D+dx*dV(M+1)))*c(M,n+1)
		time=time+dt
	end do
	do i=1,M+1
		c(i,1)=c(i,dN+1)
	end do
	write(id,*) c(1:M+1,1)
end do

close(id)
close(unit2)

end subroutine pdesolver