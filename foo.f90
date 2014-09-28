

  program main
  implicit none
!
! gfortran foo.f90 transverseintegral.f90 -o foo.out
! ./foo.out  n1  m1  n2  m2  r
!
! compute all the nonzero TMC(N,M,n,m; n1,m1,n2,m2,r=b2/b1), 
! for given n1, m1, n2, m2, r
!
  integer, parameter :: imax = 5000;
  integer :: n1, m1, n2, m2, nn, mm, n, m, ierr, nmax, jmax, e1, e2, ee, e;
  integer :: i, j, k;
  double precision :: r, val, norm, t1, t2, t3, t4;
  double precision, external :: TMC;
  character(len=32) :: arg

  integer :: q(8,imax);
  double precision :: tm(imax);

  call get_command_argument(1, arg, ierr)
  read(arg, '(i3)'), n1
  call get_command_argument(2, arg, ierr)
  read(arg, '(i3)'), m1
  call get_command_argument(3, arg, ierr)
  read(arg, '(i3)'), n2
  call get_command_argument(4, arg, ierr)
  read(arg, '(i3)'), m2
  call get_command_argument(5, arg, ierr)
  read(arg, '(f6.2)'), r 

  norm = 0d0;
  nmax = 2*n1+2*n2+abs(m1)+abs(m2)
  jmax = m1+m2

  print 28, nmax, jmax;
28 format(1x, 'Nmax = ', i4, 4x, 'Jmax = ', i4);

  i = 0
  call cpu_time(t1)
  do m = -ceiling((nmax+abs(jmax))/2.0), ceiling((nmax+abs(jmax))/2.0)
	mm = jmax - m
    do n = 0, (nmax-abs(m)-abs(mm))/2
	  !do nn = 0, (nmax-abs(m)-abs(mm))/2-n
	    nn = (nmax-abs(m)-abs(mm))/2-n
	    if(nn>=0) then
	      i = i + 1
	      val = TMC(nn, mm, n, m, n1, m1, n2, m2, r);
		  norm = norm + val*val;
		  if(i.le.imax) then
			tm(i) = val
			q(:,i) = (/nn, mm, n, m, n1, m1, n2, m2/)
		  endif
		endif
	  !enddo
	enddo
  enddo
  call cpu_time(t2)
  
  print *, 'count = ', i
  print *, 'unitarity: ' , norm, ', error = ', abs(norm-1d0)
  print *, 'mean, total time', (t2-t1)/i, (t2-t1);
  open(unit=55, file="tmc.dat", action='write')
  write(55,*), i, r
  do j = 1, i
    write(55,*), q(:,j), tm(j)
  enddo
  close(55)

  end program

