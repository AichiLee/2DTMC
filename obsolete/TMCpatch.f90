!
! patch of the Talmi-Moshinsky coefficient subroutine
! improve the numerical stability as well as the numerical efficiency 
! 
  double precision function TMC(nn, mm, n, m, n1, m1, n2, m2, x1, x2)
  implicit none
  integer, intent(in) :: nn, mm, n, m, n1, m1, n2, m2;
  double precision, intent(in) :: x1, x2;
  double precision, external :: TMBracket;

  integer :: CapN, CapM, litn, litm, na, ma, nb, mb;

  if(2*n1+m1+abs(m1) .le. 2*n2+m2+abs(m2)) then
	  if( 2*nn+mm+abs(mm) .le. 2*n+m+abs(m) ) then
	    if( 2*n1+m1+abs(m1) .le. 2*nn+mm+abs(mm) ) then
		  TMC = TMBracket(nn, mm, n, m, n1, m1, n2, m2, x1, x2);
	    else
		  TMC = TMBracket(n1, m1, n2, m2, nn, mm, n, m, x1, x2);
		endif
	  else
	    if( 2*n1+m1+abs(m1) .le. 2*n+m+abs(m) ) then
		  TMC = TMBracket(nn, mm, n, m, n1, m1, n2, m2, x1, x2);
	    else
		  TMC = TMBracket(n1, m1, n2, m2, n, m, nn, mm, x2, x1)*(1-2*mod(abs(m2),2));
		endif
	  endif
  else
	  if( 2*nn+mm+abs(mm) .le. 2*n+m+abs(m) ) then
	    if( 2*n2+m2+abs(m2) .le. 2*nn+mm+abs(mm) ) then
		  TMC =TMBracket(nn, mm, n, m, n2, m2, n1, m1, x2, x1)*(1-2*mod(abs(m),2));
	    else
		  TMC =TMBracket(n2, m2, n1, m1, nn, mm, n, m, x2, x1)*(1-2*mod(abs(m),2));
		endif
	  else
	    if( 2*n2+m2+abs(m2) .le. 2*n+m+abs(m) ) then
		  TMC =TMBracket(nn, mm, n, m, n2, m2, n1, m1, x2, x1)*(1-2*mod(abs(m),2));
	    else
		  TMC =TMBracket(n2, m2, n1, m1, n, m, nn, mm, x1, x2)*(1-2*mod(abs(m+m1),2));
		endif
	  endif
  endif

  end function

  
  double precision function TMBracket2(nn, mm, n, m, n1, m1, n2, m2, x1, x2)
  implicit none
  integer, intent(in) :: nn, mm, n, m, n1, m1, n2, m2;
  double precision, intent(in) :: x1, x2;
  integer :: a, b, r, s, g3, g, rho3, rho
  double precision :: sindelta, cosdelta, tandelta;
  double precision, external :: LogNM, LogMultinomial3, LogBinomial;

  if(2*nn+abs(mm)+2*n+abs(m).eq.2*n1+abs(m1)+2*n2+abs(m2).and.mm+m.eq.m1+m2) then

     sindelta = sqrt(x2/(x1+x2));
     cosdelta = sqrt(x1/(x1+x2));
     tandelta = sqrt(x2/x1);
	 TMBracket2 = 0d0;
	 do a = -abs(m), abs(m), 2
	 do b = -abs(mm), abs(mm), 2
	    r = abs(sign(1,m)*a+sign(1,mm)*b+m2-m1)/2
	 do g3 = max(r-n, 0), nn
	 do g = g3-nn, nn-g3, 2
	    rho = n1-n2-g+(abs(m1)-abs(m2)-a-b)/2;
	 do rho3 = max(0,r-g3)+mod(abs(max(0,r-g3)+n+rho),2), min(n,rho+n,n-rho), 2
	 
	    s = 2*(rho3+g3)+sign(1,m)*a+sign(1,mm)*b+m2-m1
		if( mod(s, 4) .eq. 0 ) then 

		  TMBracket2 = TMBracket2 +                                            &
		    (1d0-2*mod(abs(rho3+(abs(m)-a)/2),2)) * tandelta**(a+2*rho)        &
		    * exp(  LogMultinomial3((nn-g3+g)/2, (nn-g3-g)/2, g3)              &
			      + LogMultinomial3((n-rho3+rho)/2, (n-rho3-rho)/2, rho3)      &
			      + LogBinomial(abs(mm), (abs(mm)+b)/2)                        &
			      + LogBinomial(abs(m), (abs(m)+a)/2)                          &
		          + LogBinomial(rho3+g3, s/4) )

	    endif
	
	 enddo
	 enddo
	 enddo
	 enddo
     enddo

	 TMBracket2 = TMBracket2 * (1d0-2*mod(nn+n+n1+n2,2))                       &
	            * sindelta**(2*n2+abs(m2)) * cosdelta**(2*n1+abs(m1))          &
				* exp( LogNM(n1,m1)+LogNM(n2,m2)-LogNM(nn,mm)-LogNM(n,m) );

  else
	  TMBracket2 = 0d0;
  endif
 
  end function
 
    Double Precision Function LogMultinomial3(a, b, c)
    implicit none
! 4-multinomial coefficient: (a+b+c)!/(a!*b!*c!)
! notice that if any argument is negative, the total value should be zero.
    
    integer, intent(in) :: a, b, c;
	integer :: i
    double precision :: loga, logb, logc;

    loga = 0D0;
    do i = 2, a
        loga = loga - log(dble(i));
    end do

    logb = 0D0;
    do i = 2, b
        logb = logb - log(dble(i));
    end do

    logc = 0D0;
    do i = 2, c
        logc = logc - log(dble(i));
    end do
    
    LogMultinomial3 = loga + logb + logc;
    do i = 2, a+b+c
        LogMultinomial3 = LogMultinomial3 + log(dble(i));
    end do
    
    End Function


