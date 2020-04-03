program MOPSO
implicit none
include 'mpif.h'

integer rank, size, ierror, ppc
integer:: swarmsize, maxiter, ndim, i, ii, j, k, fitnum, paretoupdate
real :: randn, rand1, rand2, c1, c2, dt
character(len=100)::initialpos, formatt
double precision, dimension(:,:), allocatable::limits, gbest, fgbest, posmpi,fitvalloc
double precision, dimension(:), allocatable:: paretobest, fparetobest
double precision, dimension(:,:,:), allocatable:: pos, pos_n, vel, lbest, fitval, flbest
integer, dimension(:,:), allocatable:: loc

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

open(10,file="input-mopso",status="old")
read(10,*) swarmsize, maxiter, ndim             !Number of particles, epochs, number of parameters
allocate(limits(ndim,3))
read(10,*) limits(:,1)
read(10,*) limits(:,2)
limits(:,3) = limits(:,2)-limits(:,1)
read(10,*) initialpos
read(10,*) fitnum                               !Number of objective values

allocate(paretobest(ndim))
paretobest(:)=0.0
paretoupdate=0
allocate(fparetobest(fitnum))
fparetobest(:)=1000

allocate(pos(ndim,swarmsize,fitnum))            !position of particles per objective
allocate(pos_n(ndim,swarmsize,fitnum))           !updated position of particles per objective
allocate(vel(ndim,swarmsize,fitnum))            !velocity of particles per objective

if(initialpos == "yes") then
        do i=1,swarmsize,1
                read(10,*) pos(:,i,1)
                do j=2,fitnum,1
                        pos(:,i,j) = pos(:,i,1)
                enddo
        enddo
else
        pos(:,:,:)=0.0
endif

allocate(fitval(fitnum,swarmsize,fitnum))       !fitness values per particle per objective
allocate(lbest(ndim,swarmsize,fitnum))          !Best position of the particle per objective
allocate(gbest(ndim,fitnum))                    !Best position of the swarm per objective
allocate(flbest(fitnum,swarmsize,fitnum))       !Best fitness values of the particle per objective
allocate(fgbest(fitnum,fitnum))                 !Best fitness values of the swarm per objective
allocate(loc(fitnum,fitnum))
loc(:,:)=0
pos_n(:,:,:)=0.0
vel(:,:,:)=0.0
fitval(:,:,:)=1000
flbest(:,:,:)=1000
fgbest(:,:)=1000
lbest(:,:,:)=0.0
gbest(:,:)=0.0
read(10,*) c1, c2, dt
close(10)

ppc=swarmsize*fitnum/size
if(size*ppc .NE. swarmsize*fitnum) then
        write(*,*) "Number of particles * number of population must be divisible by the number of processors"
        call MPI_ABORT(MPI_COMM_WORLD,ierror)
endif

allocate(posmpi(ndim,ppc))               !position of particles per core per objective
allocate(fitvalloc(fitnum,ppc))          !fitness values of particles per core per objective
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
                                do k=1,fitnum,1
                                        call RANDOM_NUMBER(randn)
                                        !write(*,*) randn
                                        pos(j,i,k)=limits(j,1)+(randn * limits(j,3))
                                enddo
                        end do
                        !write(*,'(22 F8.4)') pos(:,i)
                enddo
        endif
        write(*,'(I3,A)') swarmsize," particles are initialized"
        write(*,'(I3,A)') fitnum," objectives will be tracked"
        write(*,'(I3,A)') size," processors are requested"
        lbest(:,:,:)=pos(:,:,:)
        flbest(:,:,:)=fitval(:,:,:)
        gbest(:,:)=pos(:,1,:)
        fgbest(:,:)=fitval(:,1,:)
        
        write(*,'(A,I3,A)') " Each processor will evaluate ",ppc," particles" 
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!scatter the position and fitval to all the available processors


do i=1,maxiter,1
!call MPI_SCATTER(pos,ppc*ndim*fitnum,MPI_DOUBLE_PRECISION,posmpi,ppc*ndim*fitnum,MPI_DOUBLE_PRECISION,&
call MPI_SCATTER(pos,ppc*ndim,MPI_DOUBLE_PRECISION,posmpi,ppc*ndim,MPI_DOUBLE_PRECISION,&
                & 0,MPI_COMM_WORLD,ierror)
call MPI_SCATTER(fitval,ppc*fitnum,MPI_DOUBLE_PRECISION,fitvalloc,ppc*fitnum,&
                & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        
        !Evaluating the fitness function for the scattered data in parallel
        do j=1,ppc,1
                call RANDOM_NUMBER(randn)
                call fitness(posmpi(:,j),ndim,fitnum,randn,fitvalloc(:,j))
                write(*,*) "Epoch ",i,": Particle ",j+(rank*ppc), posmpi(:,j), " is processed by rank ", rank, &
                                & " with fitness values:",fitvalloc(:,j)
                !write(*,'(A,I3,22 F8.4,A,2 F8.4)') "particle ", j+(rank*ppc), posmpi(:,j),"  Fitness: ", fitvalloc(:,j)
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        
        !gather the scattered output into large fitval
call MPI_GATHER(fitvalloc,ppc*fitnum,MPI_DOUBLE_PRECISION,fitval,ppc*fitnum,&
                &MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
        
        !check fitness in root and update master swarm
        if(rank==0) then
                !do j=1,swarmsize,1
                !        write(*,*) pos(:,j), fitval(:,j)
                !enddo
                if (i == 1) then
                        flbest=fitval
                        lbest=pos
                else
                        do j=1,swarmsize,1
                                do k=1,fitnum,1
                                        if(fitval(k,j,k) < flbest(k,j,k)) then !choose the local best position of particle based on fdim
                                                flbest(k,j,k) = fitval(k,j,k)
                                                lbest(:,j,k)=pos(:,j,k)
                                        endif
                                enddo
                        enddo
                endif

                !write(*,*) "EPOCH: ",i, " Completed!"
                write(formatt,*) "(A,I3,A,I3,A,",ndim,"F10.3,A,",fitnum,"F10.3)"
                do j=1,fitnum,1 ! per slave population
                        loc(:,j)=MINLOC(fitval(:,:,j),DIM=2)
                        fgbest(:,j)=MINVAL(fitval(:,:,j), DIM=2)
                        !write(*,*) loc(:,j)
                        !write(*,*) fgbest(:,j)
                !do j=1,fitnum,1 ! per fitness value
                        !do k=1,fitnum,1 ! per objective
                        gbest(:,j)=pos(:,loc(j,j),j)
                        !write(*,'(A,I1,A,22 F8.4,A,F8.4)') "Global best for ",j," : ", gbest(:,j), "  Fitness: ",fgbest(j)
                        !enddo
                        write(*,formatt) "EPOCH ",i,": Global best for ",j," : ",gbest(:,j), "  Fitness: ",fitval(:,loc(j,j),j)
                enddo
                !write(*,*) "EPOCH ",i,": Global best for ",fdim," : ", gbest(:,fdim), "  Fitness: ",fgbest(fdim)
                !write(formatt,*) "(A,I3,A,I3,A,",ndim,"F10.3,A,",fitnum,"F10.3)"
                !write(*,*) ADJUSTL(formatt)
                !do j=1,fitnum,1
                !        write(*,formatt) "EPOCH ",i,": Global best for ",j," : ",gbest(:,j), "  Fitness: ",fitval(:,loc(j,j),j)
                !enddo

                !update velocity of slave population
                do ii=1,fitnum,1                        !each slave population
                        do j=1,swarmsize,1
                                paretoupdate=0
                                do k=1,fitnum,1         !each objective value within the swarm of any population
                                        if(fparetobest(ii) >= fitval(k,j,ii)) then
                                                paretoupdate=paretoupdate+1
                                        endif
                                enddo
                                if(paretoupdate==fitnum) then
                                        fparetobest=fitval(:,j,ii)
                                        paretobest=pos(:,j,ii)
                                endif
                                do k=1,ndim,1
                                        CALL RANDOM_NUMBER(rand1)
                                        CALL RANDOM_NUMBER(rand2)
                                        vel(k,j,ii)= c1*rand1*((lbest(k,j,ii)-pos(k,j,ii))/dt) + &
                                                        &c2*rand2*((gbest(k,ii)-pos(k,j,ii))/dt) !can be
                                        !updated to include any number of response variables
                                        pos_n(k,j,ii)= pos(k,j,ii) + vel(k,j,ii)
                                        if(pos_n(k,j,ii) < limits(k,2) .AND. pos_n(k,j,ii) > limits(k,1)) then
                                                pos(k,j,ii) = pos_n(k,j,ii)
                                        elseif(pos_n(k,j,ii) > limits(k,2)) then
                                                CALL RANDOM_NUMBER(rand1)
                                                !pos(k,j) = limits(k,2)
                                                pos(k,j,ii) = limits(k,2)-(rand1 * limits(k,3))
                                        elseif(pos_n(k,j,ii) < limits(k,1)) then        
                                                !pos(k,j) = limits(k,1)
                                                CALL RANDOM_NUMBER(rand1)
                                                pos(k,j,ii) = limits(k,1)+(rand1 * limits(k,3))
                                        endif 
                                enddo
                        enddo
                enddo
                write(*,formatt) "EPOCH ",i,": Pareto best for ",fitnum," obj: ",paretobest, "  Fitness: ",fparetobest

        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
enddo
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

call MPI_FINALIZE(ierror)
!deallocate(limits,pos, posmpi, vel, fitval, fitvalloc, lbest, gbest, flbest, fgbest, loc)
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
fitval(1)=( -47.0669 + (0.1633*pos(1)) - (0.978167*pos(2)) + (54.0917*pos(3)) &
& -(0.000084*pos(1)*pos(1)) - (13.5625*pos(3)*pos(3)) &
& + (0.00065*pos(1)*pos(2)) - (0.023*pos(1)*pos(3)))
fitval(2)=( -37.0669 + (0.1633*pos(1)) - (0.978167*pos(2)) + (54.0917*pos(3)) &
& -(0.000084*pos(1)*pos(1)) - (13.5625*pos(2)*pos(3)) &
& + (0.00065*pos(3)*pos(2)) - (0.023*pos(2)*pos(3)))
!ndim1=ndim+1
!write(FMT,*) ndim1
!write(FMT1,'(f11.9)') randn
!!write(*,*) FMT
!!write(*,*) FMT1,"chk"
!write(command,"(A, " // ADJUSTL(FMT) // " F16.9)") "./Activation.sh ",randn,pos
!write(filename,*)"runs/run_",ADJUSTL(FMT1),"/response"
!!write(*,*) ADJUSTL(command)
!!write(*,*) filename
!call EXECUTE_COMMAND_LINE(command)
!open(10,file=ADJUSTL(filename),status="old")
!read(10,*) fitval(1), fitval(2)
!read(10,*) fitval(3), fitval(4)
!close(10)
end subroutine fitness
