
!-------------------------------------------------------------------------------!
!																				!
!           Maris basis wavefunction integrals									!
! 			Yang Li, Jul 24, 2012, leeyoung@iastate.edu							!
!																				!
!-------------------------------------------------------------------------------!

!																				!
!===============================================================================!
!																				!
!			PART I, ELEMENTARY INTEGRALS (MarisBasis.pdf Eq. 25)				!
!																				!
!===============================================================================!
!																				!
    Function Inmk(n, m, k) 
    implicit none
! 	Maris Basis, Eq.(25a)
    
    integer :: n, m, k
    integer :: ii, jj
    double precision :: Inmk, logk, LogN
    double precision,parameter :: OneOverSqrtPi = 0.56418958354775629D0;
!	LogN(n,k) = log( sqrt(n!/(n+k)!) )

	if(m+k == 0) then
        
            logk = 0D0;
            do ii = 1, k
                logk = logk + Dlog(dble(ii));
            end do
            
            Inmk = 	(-1)**mod(n,2) * 2D0**dble(k) * OneOverSqrtPi 	&
            &		* exp(-logk - LogN(n,k)); 
	else
			
			Inmk = 0D0;
    endif

    End Function Inmk
    
    
    Function Inmnm1(n, m, np, mp)
    implicit none
! 	Maris Basis, Eq.(25c)
	
	integer ::	n, m, np, mp;
	double precision ::	Inmnm1;
	
	Inmnm1 = 0D0;
	
	if(mp == m+1) then
		
		if(m >= 0) then
			
			if(n == np) then
				
				Inmnm1 = sqrt(dble(n+abs(m+1)));
			
			else if(np == n - 1) then
				
				Inmnm1 = -1D0 * sqrt(dble(n));
			
			end if
		
		else	! m < 0
			
			if(n == np) then
				
				Inmnm1 = sqrt(dble(n+abs(m)));
			
			else if(np == n + 1) then
				
				Inmnm1 = -1D0 * sqrt(dble(n+1));
			
			end if
		
		endif
	
	endif    
    
    End Function Inmnm1

!																				!
!===============================================================================!
!																				!
!			PART II. 2D Talmi-Moshinsky Transform								!
!																				!
!===============================================================================!
!																				!

    Function TMBracket(nn1, mm1, nn2, mm2, n1, m1, n2, m2, x1, x2) 
    implicit none
    
    ! talmi-moshinsky transformation, Eq. 2.20 - Eq.2.21;
    ! < N1, M1, N2, M2, sqrt(2), sqrt(1/2) | n1, m1, n2, m2, 1, 1 >
    !
    ! here I use the version coded in mathematica, slitly different 
    ! from Eq. 2.21 
    ! document this algorithm (See MarisBasis.pdf)
    ! tested with 1. Mathematica code using the same algorithm by Young, 
    !             2. Fortran code using different algorithm by Xingbo Zhao
    ! up to Nmax = 20.

    integer :: nn1, mm1, nn2, mm2, n1, m1, n2, m2, BinomialSign, MultinomialSign;
    double precision :: xx1, xx2, x1, x2, TMBracket;
    integer :: 	gamma1, gamma2, beta1, beta2, beta3, beta, v1, v2, vv1, vv2, e1,&
    &		e2, ee1, ee2, sgn;
    double precision :: sindelta, cosdelta, tandelta, summ;
    double precision :: LogMultinomial, LogBinomial, LogNM, detector;
    double precision, parameter :: Pi = 3.141592653589793D0, 					&
    &				TwoPi = 6.283185307179586D0, 								&
    &				TwoPi2 = 39.47841760435743D0,								&
    &               SqrtTwoPi = 2.506628274631001D0, 							&
    &				Sqrt2 = 1.414213562373095D0, 								&
    &				iota = 1.0D-9;

!    x1 = 1D0;
!    x2 = 1D0;
    
    xx1 = x1+x2;
    !xx2 = sqrt(x1*x2/xx1);
    
    if( 2*nn1 + abs(mm1) + 2*nn2 + abs(mm2) /= 2*n1 + abs(m1) + 2*n2 + abs(m2) 	&
    &	.or. mm1 + mm2 /= m1 + m2) then
        TMBracket = 0D0;
        return;
    end if
    
!    if(abs(ll1 - sqrt(l1**2 + l2**2)) > iota .or. abs(ll2 - l1*l2/sqrt(l1**2 + l2**2)) > iota ) then
!        TMBracket = 0D0;
!        return;
!    end if
    
	v1 = n1 + (abs(m1) - m1)/2;
	v2 = n2 + (abs(m2) - m2)/2;
	vv1=nn1 + (abs(mm1)-mm1)/2;
	vv2=nn2 + (abs(mm2)-mm2)/2;
	e1 = 2 * n1 + abs(m1);
	e2 = 2 * n2 + abs(m2);
	ee1 =2 * nn1+abs(mm1);
	ee2 =2 * nn2+abs(mm2);

	sindelta = sqrt(x2/xx1);
	cosdelta = sqrt(x1/xx1);
	tandelta = sqrt(x2/x1); 
   
	!print *, ' ------------------------ ';
	!print *, v1, v2, vv1, vv2, sindelta, cosdelta, tandelta
	!print *, ;
   if(mm2 >= 0) then
   if(mm1 >= 0) then
	!print *, 'sec 1';
   ! M1 >= 0 and M2 >= 0
   	summ = 0D0;
        do gamma1 = 0, v1
            do gamma2 = max(0, v2 - v1 + gamma1 - m1), v2
                do beta1 = 0, min(vv2, gamma1)
                    do beta2 = 0, min(vv2, gamma2)
                       do beta3 = max(0, vv2 - v2 - beta1 - beta2 + gamma2), 	&
                       &	min(min(vv2, v1 - gamma1), vv2 - beta1 - beta2)
                            do beta = max(0, v1-v2+m1+gamma2-gamma1-mm1),  		&
                       &	min(v1 - v2 + m1 + gamma2 - gamma1, mm2)
                               
                               summ = summ + (1-2*mod(beta1+beta2+beta,2)) * exp( &
                               LogMultinomial(gamma1 - beta1,                     &
                               gamma2 - beta2,                                    &
                               v1 - gamma1 - beta3,                               &
                               v2 - gamma2 -vv2 + beta1 + beta2 + beta3)          &
                            +  LogMultinomial(beta1, beta2, beta3, vv2 - beta1 - beta2 - beta3) &
                            +  LogBinomial(mm1, v1 - v2 + m1 + gamma2 - gamma1 - beta)          &
                            +  LogBinomial(mm2, beta)   )                         &
                            *  tandelta ** (2*(beta1 - beta2 + beta));
                            !	 signs of Binomial coefficients are +1, guarrenteed 
                            !   by the bounds of beta.
                            !    min(min( ... )) in beta bound is used to guarrentee 
                            !   the positivity of arguments in logmultinomial().
                            !   summ = summ + detector;
                            !   print *, detector, gamma1, gamma2, beta1, beta2, beta3, beta;

                            end do
                        end do
                    end do
                end do
            end do
        end do
    else !mm2 >= 0, mm1 < 0
	!print *, 'sec 2';
    	summ = 0D0;
        do gamma1 = 0, v1
           do gamma2 = max(0, v2 - v1 + gamma1 - m1),  v2
                do beta1 = 0, min(vv2, gamma1)
                    do beta2 = 0, min(vv2, gamma2)
                       do beta3 = max(0, vv2 - v2 - beta1 - beta2 + gamma2), 	&
                       &	min(min(vv2, v1 - gamma1), vv2 - beta1 - beta2)
					        sgn = (1-2*mod(v1+v2+gamma1+gamma2+abs(m1)+abs(mm1),2))
                            do beta = 0,  min(v1 - v2 + m1 + gamma2 - gamma1, mm2 ) 

                            !   detector = 
                               !detector = (1-2*mod(beta1+beta2+beta,2)) * exp( &
                               detector = (1-2*mod(beta1+beta2,2)) * exp( &
                               LogMultinomial(gamma1 - beta1,                     &
                               gamma2 - beta2,                                    &
                               v1 - gamma1 - beta3,                               &
                               v2 - gamma2 -vv2 + beta1 + beta2 + beta3)          &
                            +  LogMultinomial(beta1, beta2, beta3, vv2 - beta1 - beta2 - beta3) &
                            +  LogBinomial(mm1, v1 - v2 + m1 + gamma2 - gamma1 - beta)          & 
                            +  LogBinomial(mm2, beta)   )                         &
                            !*  BinomialSign(mm1,v1 - v2 + m1 + gamma2 - gamma1 - beta )			&
							*  sgn * tandelta ** (2*(beta1 - beta2 + beta));
                            !  mm1 < 0, thus may have a negative sign. 
                            ! min(min( ... )) guarrentees the positivity of arguments
                            ! in multinomial/binomial coefficients.
                                summ = summ + detector;
                                !print *, detector, summ, gamma1, gamma2, beta1, beta2, beta3, beta;
                                !print 218, detector, summ, 0, v1, v2-v1-m1+gamma1, &
							         !v2, 0, min(vv2,gamma1), 0, min(vv2,gamma2), &
									 !vv2-v2-beta1-beta2+gamma2,                &
									 !min(v1-gamma1,vv2-beta1-beta2),0,         &
									 !min(mm2,v1-v2+m1+gamma2-gamma1);          &

218 format(1x, 2g12.6, 6('[', i4, ',', i4, ']') )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    else ! mm2 < 0
    if(mm1 >= 0) then ! mm2 < 0; mm1 >= 0
	!print *, 'sec 3', mm1, mm2;
    	summ = 0D0;
        do gamma1 = 0, v1
            do gamma2 = max(0, v2 - v1 + gamma1 - m1), v2
                do beta1 = 0, min(vv2, gamma1)
                    do beta2 = 0, min(vv2, gamma2)
                       do beta3 = max(0, vv2 - v2 - beta1 - beta2 + gamma2), 	&
                       &	min(min(vv2, v1 - gamma1), vv2 - beta1 - beta2)
                            do beta = max(0, v1 - v2 + m1 + gamma2 - gamma1 - 	&
                       &		mm1),  v1 - v2 + m1 + gamma2 - gamma1

                               summ = summ + (1-2*mod(beta1+beta2+beta,2)) * exp( &
                               LogMultinomial(gamma1 - beta1,                     &
                               gamma2 - beta2,                                    &
                               v1 - gamma1 - beta3,                               &
                               v2 - gamma2 -vv2 + beta1 + beta2 + beta3)          &
                            +  LogMultinomial(beta1, beta2, beta3, vv2 - beta1 - beta2 - beta3) &
                            +  LogBinomial(mm1, v1 - v2 + m1 + gamma2 - gamma1 - beta)          &
                            +  LogBinomial(mm2, beta)                                           &
                                                        ) * BinomialSign(mm2, beta)             &
                            *  tandelta ** (2*(beta1 - beta2 + beta));
                            ! See the comments of other cases.
                            end do
                        end do
                    end do
                end do
            end do
        end do
    else ! mm2 < 0, mm1 < 0
	!print *, 'sec 4';
    summ = 0D0;
!        do gamma1 = max(0, v1 - v2 + m1), min(v1, v1 + m1) this bound seems wrong
		 do gamma1 = 0, v1
            do gamma2 = max(0, v2 - v1 + gamma1 - m1), v2
                do beta1 = 0, min(vv2, gamma1)
                    do beta2 = 0, min(vv2, gamma2)
                       do beta3 = max(0, vv2 - v2 - beta1 - beta2 + gamma2), 	&
                       &	min(min(vv2, v1 - gamma1), vv2 - beta1 - beta2)
                            do beta = 0,  v1 - v2 + m1 + gamma2 - gamma1
                               
                               summ = summ + (1-2*mod(beta1+beta2+beta,2)) * exp( &
                               LogMultinomial(gamma1 - beta1,                     &
                               gamma2 - beta2,                                    &
                               v1 - gamma1 - beta3,                               &
                               v2 - gamma2 -vv2 + beta1 + beta2 + beta3)          &
                            +  LogMultinomial(beta1, beta2, beta3, vv2 - beta1 - beta2 - beta3) &
                            +   LogBinomial(mm1, v1 - v2 + m1 + gamma2 - gamma1 - beta)         &
                            +   LogBinomial(mm2, beta)                            & 
                                                        )                         &
                            *  BinomialSign(mm1, v1 - v2 + m1 + gamma2 - gamma1 - beta)			&
                            *  BinomialSign(mm2, beta)							  &
                            *  tandelta ** (2*(beta1 - beta2 + beta));
                            ! see the comments of other cases. 
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    end if 
!       print *, "summ = ", summ;

	!print *, '';
	!print *, 'summ = ', summ;
	TMBracket = summ * (1-2*mod(vv2+mm2+nn1+nn2+n1+n2,2)) * sindelta**e2       &
	          &      * cosdelta**e1 / tandelta**mm2  * exp( logNM(n1, m1)      &
		      &        + logNM(n2, m2) - logNM(nn1, mm1) - logNM(nn2, mm2) );
    !print *, sindelta**e2, cosdelta**e1, tandelta**mm2,  logNM(n1, m1)      &
		      !&        + logNM(n2, m2) - logNM(nn1, mm1) - logNM(nn2, mm2) ;

       !
	   ! this code differs from Xingbo's code by a factor of (-1)^mm2.
	   ! The reason is different conventions:
       ! phi(p1) phi(p2) = sum TM-bracket * phi(p1 + p2) * phi((p1-p2)/2) 
       ! vs
       ! phi(p1) phi(p2) = sum TM-bracket * phi(p1 + p2) * phi((p2-p1)/2), 
       ! the last function caused a sign difference, (-1)^M2, if the m quantum 
       ! number of phi((p1-p2)/2) is M2. 
       !
       ! First posted by Young on Jan 26, 2012, updated by Young on Feb. 3, 2012
       ! modified on July 24, 2012 for Maris basis

	!print *, ' ------------------------ ';

    End Function TMBracket

  double precision function TMC(nn, mm, n, m, n1, m1, n2, m2, x1, x2)
  implicit none
!
! patch of the Talmi-Moshinsky coefficient subroutine
! improve the numerical stability as well as the numerical efficiency 
! 
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


!																				!
!===============================================================================!
!																				!
!			PART III. Wavefunction Integrals	MarisBasis.pdf Eq. 26-28		!
!																				!
!===============================================================================!
!																				!
	Function WFI2to1q1(n1, m1, x1, n2, m2, x2, np, mp, xp)
	implicit none
!	q1 + q2 -> q, with angular momentum exchange from q1, 
!	MarisBasis.pdf Eq. 26a
	
	integer ::	n1, m1, n2, m2, np, mp;
	double precision ::	x1, x2, xp, WFI2to1q1, TMBracket, Inmk, Inmnm1;
	
	integer ::	NN, mu, nu, n;
	double precision ::	a, b, c;
	
	if( m1 + m2 + 1 == mp ) then
	
		n = n1 + n2 - np + (abs(m1) + abs(m2) - abs(mp) - 1)/2;
		nu = n1 + n2 + (abs(m1) + abs(m2) - abs(m1 + m2))/2;
		
		if( n >= 0) then
		
			WFI2to1q1 = TMBracket(np, mp, n, -1, n1, m1, n2, m2, x1, x2) 		&
			&			 * sqrt(x2)/(x1+x2)**1.5D0 * Inmk(n, -1, 1);
!			print *, 'partII of WFI2to1q1=', WFI2to1q1;
			!print '(8I4, 3F10.4)', np, mp, n, -1, n1, m1, n2, m2, real(x1), real(x2), 					&
			!&		 real(TMBracket(np, mp, n, -1, n1, m1, n2, m2, x1, x2));
		else
			
			WFI2to1q1 = 0D0;
		
		endif
				
		a = 0D0;
		do NN = max(0, np-1), min(nu, np+1)
			
				a = a +  TMBracket(NN, mp-1, nu-NN, 0, n1, m1, n2, m2, x1, x2)	&
				&			* Inmk(nu-NN,0,0) * Inmnm1(NN, mp-1, np, mp);
		
		enddo
		
		WFI2to1q1 = WFI2to1q1 + sqrt(x1)/(x1+x2)**1.5D0 * a; 
	
	else
		
		WFI2to1q1 = 0D0;
	
	endif
	
	End Function


	Function WFI2to1q2(n1, m1, x1, n2, m2, x2, np, mp, xp)
	implicit none
!	q1 + q2 -> q, with angular momentum exchange from q2, 
!	MarisBasis.pdf Eq. 26b
	
	integer ::	n1, m1, n2, m2, np, mp;
	double precision ::	x1, x2, xp, WFI2to1q2, TMBracket, Inmk, Inmnm1;
	
	integer ::	NN, mu, nu, n;
	double precision ::	a, b, c;
	
	if( m1 + m2 + 1 == mp ) then
	
		n = n1 + n2 - np + (abs(m1) + abs(m2) - abs(mp) - 1)/2;
		nu = n1 + n2 + (abs(m1) + abs(m2) - abs(m1 + m2))/2;
		
		if( n >= 0) then
		
			WFI2to1q2 = -1D0 * TMBracket(np, mp, n, -1, n1, m1, n2, m2, x1, x2) &
			&			 * sqrt(x1)/(x1+x2)**1.5D0 * Inmk(n, -1, 1);
		
		else
			
			WFI2to1q2 = 0D0;
		
		endif
		
!		print *, 'partII of WFI2to1q2=', WFI2to1q2;
!		print *, np, mp, n, -1, n1, m1, n2, m2, real(x1), real(x2), 			&
!				 real(TMBracket(np, mp, n, -1, n1, m1, n2, m2, x1, x2));
		a = 0D0;
		do NN = max(0, np-1), min(nu, np+1)
			
				a = a +  TMBracket(NN, mp-1, nu-NN, 0, n1, m1, n2, m2, x1, x2)	&
				&			* Inmk(nu-NN,0,0) * Inmnm1(NN, mp-1, np, mp);
		
		enddo
		
		WFI2to1q2 = WFI2to1q2 + sqrt(x2)/(x1+x2)**1.5D0 * a; 
	
	else
		
		WFI2to1q2 = 0D0;
	
	endif
	
	End Function

	Function WFI2to1qp(n1, m1, x1, n2, m2, x2, np, mp, xp)
	implicit none
	! q1 + q2 -> q', MarisBasis Eq. 26c

	integer ::	n1, m1, n2, m2, np, mp;
	double precision ::	x1, x2, xp, WFI2to1qp, TMBracket, Inmk, Inmnm1;
	
	integer ::	NN, mu, nu;
	double precision ::	a, b, c;	
	
	if(m1 + m2 + 1 == mp) then
		
		a = 0D0;
		nu = n1 + n2 + (abs(m1) + abs(m2) - abs(m1 + m2))/2;
		do NN = max(0, np-1), min(nu, np+1)
			
			a = a + TMBracket(NN, mp-1, nu-NN, 0, n1, m1, n2, m2, x1, x2) 		&
			&		* Inmk(nu-NN,0,0) * Inmnm1(NN, mp-1, np, mp);
		
		enddo
		WFI2to1qp = a / (x1 + x2);
	
	else
		
		WFI2to1qp = 0D0;
	
	endif
	
	End Function
	
	Function WFI2to1(n1, m1, x1, n2, m2, x2, np, mp, xp)
	implicit none
	! q1 + q2 -> q', Eq. 26d, MarisBasis.pdf
	
	integer ::	n1, m1, n2, m2, np, mp, n;
	double precision ::	x1, x2, xp, WFI2to1, TMBracket, Inmk, Inmnm1;
	
	if( m1+m2 == mp ) then
		
		n = n1 + n2 - np + (abs(m1) + abs(m2) - abs(m1+m2))/2;
		
		if( n >= 0 ) then
		
			WFI2to1 = TMBracket(np, mp, n, 0, n1, m1, n2, m2, x1, x2) 			&
			&			* Inmk(n, 0, 0) / (x1 + x2);
			
		else
			
			WFI2to1 = 0D0;
		
		endif
		
	else 
		
		WFI2to1 = 0D0;
	
	endif
	
	End Function
	
	Function WFI3to1(n1, m1, x1, n2, m2, x2, n3, m3, x3, np, mp, xp)
	implicit none
	! q1 + q2 + q3 -> q' Eq. 27a
	integer ::	n1, m1, n2, m2, n3, m3, np, mp;
	double precision ::	x1, x2, x3, xp, WFI3to1, Inmk, TMBracket;
	
	integer ::	NN, chi, xi;
	double precision ::	a, b, c;
	
	if( m1 + m2 + m3 == mp ) then
		
		a = 0D0;
		xi = n1 + n2 + (abs(m1) + abs(m2) - abs(m1 + m2))/2;
		chi= n3 - np + (abs(m1 + m2) + abs(m3) - abs(mp))/2;
		
		do NN = max(0, -chi), xi
			
			a = a + TMBracket(NN, m1+m2, xi-NN, 0, n1, m1, n2, m2, x1, x2) 		&
			&		* TMBracket(np, mp, chi+NN, 0, NN, m1+m2, n3, m3, x1+x2, x3)&
			&		* Inmk(xi-NN, 0, 0) * Inmk(chi+NN, 0, 0);
		
		enddo
		
		WFI3to1 = a / (x1 + x2 + x3);
	
	else 
		
		WFI3to1 = 0D0;
	
	endif
	
	End Function
	
	Function WFI2to2(n1, m1, x1, n2, m2, x2, n1p, m1p, x1p, n2p, m2p, x2p)
	implicit none
	! q1 + q2 -> q'1 + q'2, Eq. 27b
	
	integer ::	n1, m1, n2, m2, n1p, m1p, n2p, m2p;
	double precision ::	WFI2to2, TMBracket, Inmk, x1, x2, x1p, x2p;
	
	double precision ::	a;
	integer ::	NN, xi, chi, n, np;
	
	if( m1 + m2 == m1p + m2p) then
		
		xi = n1 + n2 + (abs(m1) + abs(m2) - abs(m1+m2))/2;
		chi= n1p+n2p + (abs(m1p)+ abs(m2p)- abs(m1p+m2p))/2;
		a = 0D0;
		do NN = 0, min(xi, chi)
			
			n = xi - NN;
			np= chi - NN;
	
			a = a + TMBracket(NN, m1+m2, n, 0, n1, m1, n2, m2, x1, x2) 			&
			&		* TMBracket(NN, m1+m2, np, 0, n1p, m1p, n2p, m2p, x1p, x2p)	&
			&		* Inmk(n, 0, 0) * Inmk(np, 0, 0);
		
			!print *, ' ';
			!print *, '---inside wfi2to2---';
			!write(*,'(8I4, 2G12.4)'), NN, m1+m2, n, 0, n1, m1, n2, m2, x1, x2;
			!print *, TMBracket(NN, m1+m2, n, 0, n1, m1, n2, m2, x1, x2);
!!
			!write(*, '(8I4,2G12.4)'), NN, m1+m2, np, 0, n1p, m1p, n2p, m2p, x1p, x2p;
			!print *, TMBracket(NN, m1+m2, np, 0, n1p, m1p, n2p, m2p, x1p, x2p);
!
			!print *, '--- --- ';
		enddo 
		
		WFI2to2 = a / (x1 + x2);
	else
		
		WFI2to2 = 0D0;
	
	end if
	
	End Function

		
!																				!
!===============================================================================!
!																				!
!			PART IV. Fractorials, Binomials, Multinomials and other factors		!
!																				!
!===============================================================================!
!																				!
        
    Function LogN(n,m)
    implicit none
!	LogN(n,k) = log( sqrt(n!/(n+k)!) )
    
        double precision :: LogN, lognm 
        integer :: n, m, ii
    
        lognm = 0D0;
        do ii = n+1, n+abs(m)
            lognm = lognm + log(dble(ii));
        end do

        LogN = -0.5D0 * lognm;

    End Function LogN    
    
    Function LogNM(n, m)
    implicit none
! 	log(sqrt((n+m)!n!))
    
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
    
    Function LogMultinomial(a, b, c, d)
    implicit none
! 4-multinomial coefficient: (a+b+c+d)!/(a!*b!*c!*d!)
! notice that if any argument is negative, the total value should be zero.
! this is specified by multinomialsign(a, b, c, d)
    
    double precision :: LogMultinomial, loga, logb, logc, logd;
    integer :: a, b, c, d, i

    loga = 0D0;
    do i = 2, a
        loga = loga + log(dble(i));
    end do

    logb = 0D0;
    do i = 2, b
        logb = logb + log(dble(i));
    end do

    logc = 0D0;
    do i = 2, c
        logc = logc + log(dble(i));
    end do
    
    logd = 0D0;
    do i = 2, d
        logd = logd + log(dble(i));
    end do

    LogMultinomial = -loga - logb - logc - logd;
    do i = 2, a+b+c+d
        LogMultinomial = LogMultinomial + log(dble(i));
    end do
    
    End Function

	Function MultinomialSign(a, b, c, d)
	implicit none
	! sign of multinomial coefficient
	! +1, if 0 <= a, b, c, d <= a+b+c+d
	! 0, otherwise;
	!TODO, generalization to negative power 
	! multinomial theorem. recurrence relation from binomials
	
	integer :: a, b, c, d, MultinomialSign;
	
	if( a >= 0 .and. b >= 0 .and. c >= 0 .and. d >= 0) then
		
		MultinomialSign = +1;
	
	else 
		
		MultinomialSign = 0;
		
	end if
	
	End Function    
