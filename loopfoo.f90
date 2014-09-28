

  program main
  implicit none
!
! gfortran loopfoo.f90 TMC.f90 -fopenmp -o foo.out 
! export OMP_NUM_THREADS=8  #allocate 8 threads
! ./foo.out  nmax  r
!
! compute all the nonzero TMC(N,M,n,m; n1,m1,n2,m2,r=b2/b1), 
! for given nmax and r
!
  integer :: n1, m1, n2, m2, nn, mm, n, m
  integer :: nmax, jmax, ncore;
  integer :: i, j, k, ierr
  double precision :: r, val, norm, sm, mxdv
  double precision :: t0, t1, t2, t3, t4
  character(len=32) :: chnmax, chr, tag, arg

  integer, external :: omp_get_num_threads, omp_get_max_threads
  double precision, external :: TMC, omp_get_wtime

  call get_command_argument(1, chnmax, ierr)
  read(chnmax, '(i3)'), nmax
  call get_command_argument(2, chr, ierr)
  read(chr, '(f6.2)'), r 

  print *, 'Nmax = ', nmax, ', r = ', r

  call cpu_time(t0)
  t1 = omp_get_wtime()

  ncore = omp_get_max_threads()

  i = 0
  j = 0
  sm = 0d0
  mxdv = 0d0

!$OMP  parallel do default(shared)              &
!$OMP& private(m1,m2,jmax,n1,n2,m,mm,n,nn,val,norm)  &
!$OMP& schedule(static)                         &
!$OMP& reduction(+:i,j,sm,mxdv)
  do m1 = -nmax, nmax
  do m2 = -nmax+abs(m1), nmax-abs(m1)
    jmax = m1+m2
  do n1 = 0, (nmax-abs(m1)-abs(m2))/2
     n2 = (nmax-abs(m1)-abs(m2))/2-n1
	 if(2*(n1+n2)+abs(m1)+abs(m2) == nmax) then
       i = i + 1
       norm = 0d0
       do m = -ceiling((nmax+abs(jmax))/2.0), ceiling((nmax+abs(jmax))/2.0)
	      mm = jmax - m
          do n = 0, (nmax-abs(m)-abs(mm))/2
	         nn = (nmax-abs(m)-abs(mm))/2-n
	         if(2*(nn+n)+abs(m)+abs(mm) == nmax) then
			   j = j + 1
	           val = TMC(nn, mm, n, m, n1, m1, n2, m2, r);
		       norm = norm + val*val;
		     endif
	      enddo
        enddo
		sm = sm + (norm-1d0)**2
		mxdv = max(mxdv, abs(norm-1d0))
	 endif
  enddo
  enddo
  enddo
!$OMP  end parallel do

  call cpu_time(t3)
  t2 = omp_get_wtime()
  
  print *, 'count: lhs/total = ', i, j
  print *, 'ncore = ', ncore
  print *, 'deviation from unity: SD/max = ' , sqrt(sm/max(1,(j-1))), mxdv
  print *, 'total time, all/wall=', (t3-t0), t2-t1;

  open(unit=73, file='tmc_stat_omp.dat', access='append')
  write(73, *), nmax, j, ncore, r, t3-t0, t2-t1, sqrt(sm/max(1,(j-1))), mxdv
  close(73)
  end program

