 !----------------------------------------------------------------------------
 !             Tamli-Moshinsky Transformation Coefficients                   ! 
 !                            for                                            !
 !                  2D Harmanic Oscillator Basis                             !
 !                                                                           !
 !    Yang, Li <leeyoung@iastate.edu>, Nov. 18, 2013                         !
 !    Contributors: Xingbo Zhao, Paul Wiecki                                 !
 !--                                                                         !
 !  definition of the Talmi-Moshinsky (TM) transformation:                   !
 !  (see tmc.pdf for more details)                                           !
 !- Let psi(n,m, p) be the normalized 2D HO functions. Then TM transformation!
 !  relates the HO function in the single-particle coordinate with the       !
 !  relative coordinate:                                                     ! 
 !                                                                           !
 !  psi(n1,m1,p1;b1)*psi(n2,m2,p2;b2) =                                      !
 !                      sum_{N,M,n,m} TMC(N,M,n,m;n1,m1,n2,m2;b1,b2)         !
 !                                          * psi(N,M,P;B)*psi(n,m,p,b)      !
 !  where B = sqrt(b1^2+b2^2), b = b1*b2/B; P = p1 + p2 [vector],            !
 !  p = (b2^2*p1 - b1^2*p2)/(b1^2+b2^2) [vector].                            !
 !                                                                           !
 !  (By doing the TM transformation reccursively, one  can relate the single-!
 !  particle coordinate HO function to the Jacobi coordinate HO function.)   !
 !
 !- The HO function, expressed in the momentum space, reads,                 !
 !  b*psi(n,m,p) = sqrt(4*pi*n!/(n+|m|)!) * exp(i*m*phi-rho^2/2) * rho^|m| * !
 !                LaguerreL(n, |m|, rho^2)                                   !
 !   where LaguerreL(n, a, x) is the generalized Laguerre polynomial.        !
 !   rho = |p|/b, phi = arg p, b is the basis energy scale.                  !
 !   [recall for H = p^2/(2*m) + 1/2*m*w^2*r^2, b = sqrt(m*w), so            !
 !    in the transformation, B = sqrt((m1+m2)*w), b = sqrt(m1*m2/(m1+m2)*w)] !
 !                                                                           !
 !  psi(n,m,p) psisatisfy the orthonormality relations:                      !
 !                                                                           !
 !  int d^2 p/(2*pi)^2 psi(n,m,p) psi^*(n',m',p) = delta(n,n') delta(m,m')   !
 !  sum_{n,m} psi(n,m,p) psi(n,m,p') -> (2*pi)^2*delta^2(p-p')               !
 !                                                                           !
 !  where "^*" means the complex conjugate.                                  ! 
 !                                                                           !
 !- Properties of the TM coefficients (TMC):                                 !
 !                                                                           !
 !  * TMC only depends on delta = arctan(b2/b1);                             !
 !  * TMC is proportional to delta(2*n1+|m1|+2*n2+|m2|, 2*N+|M|+2*n+|m|);    !
 !  * TMC is proportional to delta(m1+m2, M+m);                              !
 !                                                                           !
 !---------------------------------------------------------------------------!

  double precision function TMC(nn, mm, n, m, n1, m1, n2, m2, tandelta)
  implicit none
  !
  ! TMC(nn, mm, n, m; n1, m1, n2, m2; b2/b1)
  !
  ! we actuall only need Binomial coefficient for n>0,0<=m<=n here
  integer, intent(in) :: nn, mm, n, m, n1, m1, n2, m2
  double precision, intent(in) :: tandelta
  integer :: a, b, r, s, g3, g, rho3, rho
  double precision :: sindelta, cosdelta
  double precision, external :: LogNM, LogMultinomial3, LogBinomial

  if(2*nn+abs(mm)+2*n+abs(m)==2*n1+abs(m1)+2*n2+abs(m2).and.mm+m==m1+m2) then

	 sindelta = sin(atan(tandelta))
	 cosdelta = cos(atan(tandelta))

	 TMC = 0d0;
	 do a = -abs(m), abs(m), 2
	 do b = -abs(mm), abs(mm), 2
	    r = abs(sign(1,m)*a+sign(1,mm)*b+m2-m1)/2
	 do g3 = max(r-n, 0), nn
	 do g = g3-nn, nn-g3, 2
	    rho = n1-n2-g+(abs(m1)-abs(m2)-a-b)/2;
	 do rho3 = max(0,r-g3)+mod(abs(max(0,r-g3)+n+rho),2), min(n,rho+n,n-rho), 2
	 
	    s = 2*(rho3+g3)+sign(1,m)*a+sign(1,mm)*b+m2-m1
		if( mod(s, 4) .eq. 0 ) then 

		  TMC = TMC +                                                          &
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

	 TMC = TMC * (1d0-2*mod(nn+n+n1+n2,2))                                     &
	            * sindelta**(2*n2+abs(m2)) * cosdelta**(2*n1+abs(m1))          &
				* exp( LogNM(n1,m1)+LogNM(n2,m2)-LogNM(nn,mm)-LogNM(n,m) );

  else
	  TMC = 0d0;
  endif
 
  end function
 
!==============================================================================
!   Auxiliary functions: factorials, binomials, multinomials etc              !
!==============================================================================
 
    Function LogNM(n, m)
    implicit none
! 	log(sqrt((n+abs(m))!n!))
    
    double precision :: logNM, logn
    integer :: n, m, i

    logn = 0D0;
    do i = 2, n
        logn = logn + log(dble(i));
    end do

    logNM = logn;
    do i = n+1, n+abs(m)
        logNM = logNM + log(dble(i));
    end do 
    
    logNM = 0.5D0 * (logNM + logn);
    
    End Function logNM

    Function LogBinomial(n,m)
    implicit none
    ! generalized binomial coefficients
    ! We follow the definition of Kronenburg ( arXiv:1105.3689v1 [math.CO], 2011)
    !
    ! for n >=0, 
    !	if 0 <= m <= n, Binomial(n, m) = n!/(m! * (n-m)!); 
    !	otherwise Binomial(n, m) = 0;
    ! for n < 0, 
    !	if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    !	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
    !	otherwise Binomial(n, m) = 0; 
    !
    ! the symmetry Binomial(n, m) = Binomial(n, n-m) remains valid.
    !
    ! the binomial theorem remains valid:
    ! ( x + y )^n = sum_{k=0}^{infty} Binomial(n, k) x^k y^{n-k};
    !
    ! for n >=0, the summation terminates at k = n;
    ! for n < 0, the summation keeps going to infinity.
    ! of course, if n < 0, the above series only converges when |x| < |y|,
    ! when |x| > |y|, we should use:
    ! ( x + y )^n = sum_{k=0}^{infty} Binomial(n, n-k) x^{n-k} y^k 
    ! 		  = sum_{k'=-infty}^{n} Binomial(n,k') x^k' y^{n-k'}
    ! 
    ! This function evaluate log(|Binomial(n, m)|); namely, the sign of 
    ! binomial(n,m) is dropped;
    ! sign is given by a separate function BinomialSign(n, m);
    ! 
    ! tested via mathematica Binomial[] function, up to n, m = -100, 500.
    !
    ! Young Li, Feb. 2, 2012, Groundhog Day. Phil saw his shadow, six more weeks
    ! winter.
          	    
        double precision :: LogBinomial
        integer :: n, m, k, ii, jj
!        double precision, parameter :: iota = -1D12;
        
        LogBinomial = 0D0;
        if ( n >= 0 ) then
			if( m >= 0 .and. m <= n) then
				
        		do ii =  1, m
            		LogBinomial = LogBinomial + 					&
            	&		log(dble(n - m + ii)) - log(dble(ii));
        		end do
				
			else
				return;
			end if
        else 
        
    !	n < 0, if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    		if ( m >= 0 ) then
    		
        		do ii =  1, m
            		LogBinomial = LogBinomial + 					&
            	&		log(dble( ii - 1 - n )) - log(dble(ii));
        		end do
   
    !	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
        	else if ( m <= n ) then
        		
        		do ii = 1, -1 - n
        			LogBinomial = LogBinomial + 					&
        		&		log(dble( n - m + ii )) - log(dble(ii));
        		end do 

    !	otherwise Binomial(n, m) = 0;          		        		
        	else
        		return;
        	end if
        
        end if
        
    End Function

	Function BinomialSign(n, m)
	implicit none
    ! sign of generalized binomial coefficients
    ! We follow the definition of Kronenburg ( arXiv:1105.3689v1 [math.CO], 2011)
    !
    ! for n >=0, 
    !	if 0 <= m <= n, Binomial(n, m) = n!/(m! * (n-m)!); 
    !	otherwise Binomial(n, m) = 0;
    ! for n < 0, 
    !	if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    !	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
    !	otherwise Binomial(n, m) = 0; 	
    ! see comments in LogBinomial() for more information.
    integer :: n, m, BinomialSign;
    
    if( n >= 0 ) then 
    	
    	if( m >= 0 .and. m <= n ) then
    		
    		BinomialSign = 1;
    		
    	else
    		BinomialSign = 0;
    	end if
    else ! n < 0
    	
    	if ( m >= 0 ) then 
    		
    		BinomialSign = (-1)**(mod(m,2));
    	
    	else if ( m <= n) then
    		
    		BinomialSign = (-1)**(mod(n-m,2));
    		
    	else
    		
    		BinomialSign = 0;
    	end if
    	
    end if
    
    End Function            
    

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


