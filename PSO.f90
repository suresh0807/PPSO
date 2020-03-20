program PSO
implicit none
include 'mpif.h'

integer rank, size, ierror, ppc
integer:: swarmsize, maxiter, ndim, i, j, k, fitnum, fdim
real :: randn, rand1, rand2, c1, c2, dt
character*32::initialpos, formatt
double precision, dimension(:,:), allocatable::limits, pos, pos_n, posmpi, vel, fitval, fitvalloc, lbest, gbest, flbest
double precision, dimension(:), allocatable:: fgbest
integer, dimension(:), allocatable:: loc

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

open(10,file="input-pso",status="old")
read(10,*) swarmsize, maxiter, ndim
allocate(limits(ndim,3))
allocate(pos(ndim,swarmsize))
allocate(pos_n(ndim,swarmsize))
allocate(vel(ndim,swarmsize))
read(10,*) limits(:,1)
read(10,*) limits(:,2)
limits(:,3) = limits(:,2)-limits(:,1)
read(10,*) initialpos
if(initialpos == "yes") then
        do i=1,swarmsize,1
                read(10,*) pos(:,i)
        enddo
else
        pos(:,:)=0.0
endif
read(10,*) fitnum
allocate(fitval(fitnum,swarmsize))
allocate(lbest(ndim,swarmsize))
allocate(gbest(ndim,fitnum))
allocate(flbest(fitnum,swarmsize))
allocate(fgbest(fitnum))
allocate(loc(fitnum))
loc(:)=0
pos_n(:,:)=0.0
vel(:,:)=0.0
fitval(:,:)=1000
flbest(:,:)=1000
fgbest(:)=1000
lbest(:,:)=0.0
gbest(:,:)=0.0
read(10,*) c1, c2, fdim, dt
close(10)

ppc=swarmsize/size
if(size*ppc .NE. swarmsize) then
        write(*,*) "Number of particles must be divisible by the number of processors"
        call MPI_ABORT(MPI_COMM_WORLD,ierror)
endif

allocate(posmpi(ndim,ppc))
allocate(fitvalloc(fitnum,ppc))
posmpi(:,:)=0.0
fitvalloc(:,:)=1000
!write(*,*) swarmsize, maxiter, ndim
!write(*,'(22 F8.4)') limits(:,1)
!write(*,'(22 F8.4)') limits(:,2)
!write(*,'(22 F8.4)') limits(:,3)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

if(rank==0) then
        !Initialize Particles
        if(initialpos == "no") then
                do i=1,swarmsize,1
                        do j=1,ndim,1
                                call RANDOM_NUMBER(randn)
                                !write(*,*) randn
                                pos(j,i)=limits(j,1)+(randn * limits(j,3))
                        end do
                        !write(*,'(22 F8.4)') pos(:,i)
                enddo
        endif
        write(*,'(I3,A)') swarmsize," particles are initialized"
        write(*,'(I3,A)') size," processors are requested"
        lbest(:,:)=pos(:,:)
        do j=1,fitnum,1
                flbest(j,:)=fitval(j,:)
                fgbest(j)=fitval(j,1)
                gbest(:,j)=pos(:,1)
        enddo
        
        write(*,'(A,I3,A)') " Each processor will evaluate ",ppc," particles" 
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!scatter the position and fitval to all the available processors


do i=1,maxiter,1
call MPI_SCATTER(pos, ppc*ndim, MPI_DOUBLE_PRECISION, posmpi, ppc*ndim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
call MPI_SCATTER(fitval, ppc*fitnum, MPI_DOUBLE_PRECISION, fitvalloc, ppc*fitnum, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        
        !Evaluating the fitness function for the scattered data in parallel
        do j=1,ppc,1
                call RANDOM_NUMBER(randn)
                call fitness(posmpi(:,j),ndim,fitnum,randn,fitvalloc(:,j))
                write(*,*) "particle ", j+(rank*ppc), " is processed by rank ", rank, " with fitness values:", fitvalloc(fdim,j)
                !write(*,'(A,I3,22 F8.4,A,2 F8.4)') "particle ", j+(rank*ppc), posmpi(:,j),"  Fitness: ", fitvalloc(:,j)
                call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        
        !gather the scattered output into large fitval
call MPI_GATHER(fitvalloc, ppc*fitnum, MPI_DOUBLE_PRECISION, fitval, ppc*fitnum, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
        
        !check fitness in root
        if(rank==0) then
                !do j=1,swarmsize,1
                !        write(*,*) pos(:,j), fitval(:,j)
                !enddo
                if (i == 1) then
                        do j=1,swarmsize,1
                                do k=1,fitnum,1
                                        flbest(k,j)= fitval(k,j)
                                enddo
                                lbest(:,j)=pos(:,j)
                        enddo
                else
                        do j=1,swarmsize,1
                                if(fitval(fdim,j) < flbest(fdim,j)) then !choose the local best position of particle based on fdim
                                        flbest(fdim,j) = fitval(fdim,j)
                                        lbest(:,j)=pos(:,j)
                                endif
                        enddo
                endif

                !write(*,*) "EPOCH: ",i, " Completed!"
                
                loc=MINLOC(fitval,DIM=2)
                fgbest=MINVAL(fitval, DIM=2)
                do j=1,fitnum,1 
                        gbest(:,j)=pos(:,loc(j))
                        !write(*,'(A,I1,A,22 F8.4,A,F8.4)') "Global best for ",j," : ", gbest(:,j), "  Fitness: ",fgbest(j)
                enddo
                !write(*,*) "EPOCH ",i,": Global best for ",fdim," : ", gbest(:,fdim), "  Fitness: ",fgbest(fdim)
                write(formatt,*) "(A,I3,A,I3,A,",ndim,"F10.3,A,",fitnum,"F8.4)"
               !write(*,*) ADJUSTL(formatt)
                write(*,formatt) "EPOCH ",i,": Global best for ",fdim," : ",gbest(:,fdim), "  Fitness: ",fgbest
                
                !update velocity
               
                do j=1,swarmsize,1
                        do k=1,ndim,1
                                CALL RANDOM_NUMBER(rand1)
                                CALL RANDOM_NUMBER(rand2)
                                vel(k,j)= c1*rand1*((lbest(k,j)-pos(k,j))/dt) + c2*rand2*((gbest(k,fdim)-pos(k,j))/dt) !can be
                                !updated to include any number of response variables
                                pos_n(k,j)= pos(k,j) + vel(k,j)
                                if(pos_n(k,j) < limits(k,2) .AND. pos_n(k,j) > limits(k,1)) then
                                        pos(k,j) = pos_n(k,j)
                                elseif(pos_n(k,j) > limits(k,2)) then
                                        CALL RANDOM_NUMBER(rand1)
                                        !pos(k,j) = limits(k,2)
                                        pos(k,j) = limits(k,2)-(rand1 * limits(k,3))
                                elseif(pos_n(k,j) < limits(k,1)) then        
                                        !pos(k,j) = limits(k,1)
                                        CALL RANDOM_NUMBER(rand1)
                                        pos(k,j) = limits(k,1)+(rand1 * limits(k,3))
                                endif 
                        enddo
                enddo
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

call MPI_FINALIZE(ierror)
deallocate(limits,pos, posmpi, vel, fitval, fitvalloc, lbest, gbest, flbest, fgbest, loc)
end program      

subroutine fitness(pos,ndim,fitnum,randn,fitval)
integer:: ndim, ndim1, fitnum
real :: randn
double precision,dimension(ndim),intent(in)::pos
double precision,dimension(fitnum),intent(out)::fitval
character(len=1000) :: command
character(len=100) :: filename
character(len=20) :: FMT 
character(len=11) ::  FMT1
!write(*,*)
!write(*,'(24 F8.4)') pos, fitval
!fitval(1)=pos(1) + pos(2)
!fitval(2)=pos(1) + pos(2)
!fitval(1)=( -47.0669 + (0.1633*pos(1)) - (0.978167*pos(2)) + (54.0917*pos(3)) - (0.000084*pos(1)*pos(1)) - (13.5625*pos(3)*pos(3)) &
!& + (0.00065*pos(1)*pos(2)) - (0.023*pos(1)*pos(3)))
!fitval(1)=abs((pos(1)-1526.0) * (pos(2)-163.5) * (pos(3)-0.86))
ndim1=ndim+1
write(FMT,*) ndim1
write(FMT1,'(f11.9)') randn
!write(*,*) FMT
!write(*,*) FMT1,"chk"
write(command,"(A, " // ADJUSTL(FMT) // " F16.9)") "./Activation.sh ",randn,pos
write(filename,*)"runs/run_",ADJUSTL(FMT1),"/response"
!write(*,*) ADJUSTL(command)
!write(*,*) filename
call EXECUTE_COMMAND_LINE(command)
open(10,file=ADJUSTL(filename),status="old")
read(10,*) fitval(1), fitval(2)
read(10,*) fitval(3), fitval(4)
close(10)
end subroutine fitness
