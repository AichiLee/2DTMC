(*********************************************
 * Talmi Moshinsky Transformation 
 * Mathematica package
 * Young Li, 10/05/2011, leeyoung@iastate.edu
 * version 0.0.0
 *********************************************)

 BeginPackage["TalmiMoshinsky`"]

 (*TODO*)

TalmiMoshinsky::usage =
"This package provides the calculation of Talmi-Moshinsky transformation brackets."

JacobiVariable::usage =
"JacobiVariable[varlist or number, opts->optlist] gives the Jacobi variables (mass list, coordinate list, momentum list)of the given list of variables. Argument varlist represents the list of particle masses. A number is also allowed, yet which represents the number of particles with the same mass 1. opts could be Momentum and Coordinate. optlist represents the list of symbols of particle momenta or coordinates. If not specified, Null is assumed."

TransformationBracket::usage = 
"TransformationBracket[N, M, n, m, Subscript[n, 1], Subscript[m, 1], Subscript[n, 2], Subscript[m, 2], L1, L2, l1, l2 ] gives the Talmi-Moshinsky transformation bracket of two Harmonic Oscillator wavefunctions."

AngleBracket::usage = 
"\[LeftAngleBracket] N1, M1, N2, M2, L1, L2, delimiter, n1, m1, n2, m2, l1, l2 \[RightAngleBracket] is a nicer input form of transformation bracket."

Prefactor::usage = 
"Prefactor[N1,M1,N2,M2,n1,m1,n2,m2,L1,L2,l1,l2] = (-1)^(N1+N2+n1+n2)*Sqrt[((n1+Abs[m1])! * (n2+Abs[m2])!*n1!*n2!)/((N1+Abs[M1])!*(N2+Abs[M2])!*N1!*N2!)]*L1*L2/l1/l2 gives the prefactor for Talmi-Moshinsky Transformations."

CheckBracket::usage="check brackets"

TM::usage="TM[N, M, n, m, n1, m1, n2, m2, x1, x2] calculates the 2D Talmi-Moshinsky coefficient with generalized binomial and multinomial coefficients"

TMC::usage="TMC[N, M, n, m, n1, m1, n2, m2, x1, x2] calculates the 2D Talmi-Moshinsky coefficient without using the generalized binomial and multinomial coefficients"

checkTM::usage="checkTM[nmax, x1, x2]"
checkvars::usage="checkvars[nmax, x1, x2]"

vars::usage="vars[nmax,x1,x2] all nonzero coefficients for nmax "

Begin["TalmiMoshinsky`Private`"]

JacobiVariable::ivdpar = "The argument `1` is not a valid expression. A list or a positive integer is expected.";
JacobiVariable::ivdlen = "The argument `1` is not valid in length.";

TalmiMoshinsky[] = Module[{}, Print["TODO: TalmiMoshinksy Package Welcome Information and Initialization"];];

JacobiVariable[mass_:1, OptionsPattern[{Momentum -> Null, Coordinate -> Null }]]:=
    Module[{mlist, plist, qlist, n, i, q1, q2, p1, p2, m1, m2},
	
	If[ListQ[mass], n = Length[mass]; mlist = mass,
	    If[ NumericQ[mass] && IntegerQ[mass] && Positive[mass], n = mass; mlist = Table[1,{mass}],	
	    Message[JacobiVariable::ivdpar, mass];Return[] 
	    ]
	];

	If[TrueQ[OptionValue[Momentum] =!= Null], 
	    If[ Length[ OptionValue[Momentum] ] != n,
		Message[JacobiVariable::ivdlen, OptionValue[Momentum]]; Return[]
	    ]
	];
	If[TrueQ[OptionValue[Coordinate] =!= Null], 
	    If[ Length[ OptionValue[Coordinate] ] != n,
		Message[JacobiVariable::ivdlen, OptionValue[Coordinate]]; Return[]
	    ]
	];
	
	qlist = OptionValue[Coordinate];
	plist = OptionValue[Momentum];
	    
	For[i = 1, i < n, i++,
	    
	    m1 = mlist[[i]];
	    m2 = mlist[[i+1]];
	    mlist[[i]] = (m1 * m2)/(m1 + m2);
	    mlist[[i+1]] = m1 + m2;
	    
	    If[TrueQ[OptionValue[Momentum] =!= Null],
		p1 = plist [[i]];
		p2 = plist [[i+1]];
		plist [[i]] = (m2*p1 - m1*p2)/(m1+m2);
		plist [[i+1]] = p1 + p2;
	    ];
	    
	    If[TrueQ[OptionValue[Coordinate] =!= Null],
		q1 = qlist [[i]];
		q2 = qlist [[i+1]];
		qlist [[i]] = q1 - q2;
		qlist [[i+1]] = (m1*q1 + m2*q2)/(m1 + m2);
	    ]
	];

	Return[{mlist,qlist, plist}];
    ]


Prefactor[N1_Integer?( # >= 0 &), M1_Integer, N2_Integer?(# >= 0 &), M2_Integer, n1_Integer?(# >= 0 &), m1_Integer, n2_Integer?(# >= 0 &), m2_Integer,L1_:Sqrt[2], L2_:1/Sqrt[2], l1_:1, l2_:1 ] := (-1)^(N1+N2+n1+n2)*(L1*L2/l1/l2)*Sqrt[((n1+Abs[m1])!*(n2+Abs[m2])!*n1!*n2!)/((N1+Abs[M1])!*(N2+Abs[M2])!*N1!*N2!)]; 

TransformationBracket[N1_, M1_, N2_, M2_, n1_, m1_, n2_, m2_, L1_: Sqrt[2], L2_: 1/Sqrt[2], l1_: 1, l2_: 1] := 
 Module[{v1, v2, V1, V2, SinDelta, CosDelta, TanDelta, e1, e2, E1, E2, sum, summ, detector},
(*  Print[{N1,M1,N2,M2,n1,m1,n2,m2,L1,L2,l1,l2}]; *)   
 
  If[
   2 N1 + Abs[M1] + 2 N2 + Abs[M2] != 2 n1 + Abs[m1] + 2 n2 + Abs[m2] || m1 + m2 != M1 + M2, Return[0]
   ];

  v1 = n1 + (Abs[m1] - m1)/2;
  v2 = n2 + (Abs[m2] - m2)/2; 
  V2 = N2 + (Abs[M2] - M2)/2; 
  V1 = N1 + (Abs[M1] - M1)/2;
  e1 = 2 n1 + Abs[m1];
  e2 = 2 n2 + Abs[m2]; 
  E2 = 2 N2 + Abs[M2]; 
  E1 = 2 N1 + Abs[M1];

(*  Print[{v1,v2,V1,V2,e1,e2,E1,E2}]; *)

  SinDelta = (l2/L1);
  CosDelta = (l1/L1);
  TanDelta = (l2/l1);

(*  Print[{SinDelta,CosDelta,TanDelta}]; *)
  
  sum = (-1)^(V2 + M2) * SinDelta^e2 * CosDelta^(e1)*1/TanDelta^(M2) * ( 
    If[M2 >= 0,
     If[M1 >= 0,
       Sum[
         (-1)^(beta1 - beta2 + beta)*
        Multinomial[gamma1 - beta1, gamma2 - beta2, 
         v1 - gamma1 - beta3, v2 - gamma2 - V2 + beta1 + beta2 + beta3]*
        Multinomial[beta1, beta2, beta3, V2 - beta1 - beta2 - beta3]*
        Binomial[M1, v1 - v2 + m1 + gamma2 - gamma1 - beta]*
        Binomial[M2, beta]  * TanDelta^(2*(beta1 - beta2 + beta))

       , {gamma1, 0, v1}
       , {gamma2, Max[0, v2 - v1 + gamma1 - m1], v2}
       , {beta1, 0, Min[V2, gamma1]}
       , {beta2, 0, Min[V2, gamma2]}
       , {beta3, Max[0, V2 - v2 - beta1 - beta2 + gamma2], Min[V2, v1 - gamma1]}
       , {beta, Max[0, v1 - v2 + m1 + gamma2 - gamma1 - M1 - 1 ], 
        Min[v1 - v2 + m1 + gamma2 - gamma1, M2]}
       ]
      ,
      Sum[
        (-1)^(beta1 - beta2 + beta)*
        Multinomial[gamma1 - beta1, gamma2 - beta2, 
         v1 - gamma1 - beta3, v2 - gamma2 - V2 + beta1 + beta2 + beta3]*
        Multinomial[beta1, beta2, beta3, V2 - beta1 - beta2 - beta3]*
        Binomial[M1, v1 - v2 + m1 + gamma2 - gamma1 - beta]*
        Binomial[M2, beta]  *TanDelta^(2*(beta1 - beta2 + beta))

       (* If[detector == 0, (*Print[{detector, gamma1, gamma2, beta1, beta2, beta3, beta, Multinomial[gamma1 - beta1, gamma2 - beta2, 
         v1 - gamma1 - beta3, v2 - gamma2 - V2 + beta1 + beta2 + beta3],
         Multinomial[beta1, beta2, beta3, V2 - beta1 - beta2 - beta3], 
          Binomial[M1, v1 - v2 + m1 + gamma2 - gamma1 - beta],
        Binomial[M2, beta] }]*) 1, Print[{detector, gamma1, gamma2, beta1, beta2, beta3, beta,"|",M1, v1-v2+m1+gamma2-gamma1-beta,Binomial[M1, v1-v2+m1+gamma2-gamma1-beta]}]];

        detector *)

       , {gamma1, 0, v1}
       , {gamma2, Max[0, v2 - v1 + gamma1 - m1], v2}
       , {beta1, 0, Min[V2, gamma1]}
       , {beta2, 0, Min[V2, gamma2]}
       , {beta3, Max[0, V2 - v2 - beta1 - beta2 + gamma2], Min[V2, v1 - gamma1]}
       , {beta, 0, Min[v1 - v2 + m1 + gamma2 - gamma1, M2]}
       ]
      ]
     ,
     If[M1 > 0,
      Sum[
       (-1)^(beta1 - beta2 + beta)*
        Multinomial[gamma1 - beta1, gamma2 - beta2, 
         v1 - gamma1 - beta3, v2 - gamma2 - V2 + beta1 + beta2 + beta3]*
        Multinomial[beta1, beta2, beta3, V2 - beta1 - beta2 - beta3]*
        Binomial[M1, v1 - v2 + m1 + gamma2 - gamma1 - beta]*
        Binomial[M2, beta]  *TanDelta^(2*(beta1 - beta2 + beta))
       , {gamma1, 0, v1}
       , {gamma2, Max[0, v2 - v1 + gamma1 - m1], v2}
       , {beta1, 0, Min[V2, gamma1]}
       , {beta2, 0, Min[V2, gamma2]}
       , {beta3, Max[0, V2 - v2 - beta1 - beta2 + gamma2], Min[V2, v1 - gamma1]}
       , {beta, Max[0, v1 - v2 + m1 + gamma2 - gamma1 - M1 - 1], 
        v1 - v2 + m1 + gamma2 - gamma1}
       ]
      ,
      
      Sum[
       (-1)^(beta1 - beta2 + beta)*
        Multinomial[gamma1 - beta1, gamma2 - beta2, 
         v1 - gamma1 - beta3, v2 - gamma2 - V2 + beta1 + beta2 + beta3]*
        Multinomial[beta1, beta2, beta3, V2 - beta1 - beta2 - beta3]*
        Binomial[M1, v1 - v2 + m1 + gamma2 - gamma1 - beta]*
        Binomial[M2, beta] * TanDelta^(2*(beta1 - beta2 + beta))
       , {gamma1, 0, v1}
       , {gamma2, Max[0, v2 - v1 + gamma1 - m1], v2}
       , {beta1, 0, Min[V2, gamma1]}
       , {beta2, 0, Min[V2, gamma2]}
       , {beta3, Max[0, V2 - v2 - beta1 - beta2 + gamma2], Min[V2, v1 - gamma1]}
       , {beta, 0, v1 - v2 + m1 + gamma2 - gamma1}
       ]
      
      ]
     
     ]);
(* 	Print[summ]; *)
  
  Return[sum*(-1)^(N1 + N2 - n1 - n2)*
    Sqrt[((n1 + Abs[m1])!*(n2 + Abs[m2])!*n1!* n2!)/((N1 + Abs[M1])!*(N2 + Abs[M2])!*N2!*N1!)]];
  ]/;(l1 > 0 && l2 > 0 && L1 == Sqrt[l1^2 + l2^2] && L2 == 1/Sqrt[1/l1^2 + 1/l2^2]) (*The selection rule for natural lengths*)

(* (l1 > 0 && l2 > 0 && L1 == Sqrt[l1^2 + l2^2] && L2 == 1/Sqrt[1/l1^2 + 1/l2^2]) (*The selection rule for natural lengths*) *)


AngleBracket[ {N1_, M1_, N2_, M2_, L1_:Sqrt[2], L2_:1/Sqrt[2]}, delimiter_, {n1_, m1_, n2_, m2_, l1_:1, l2_:1} ]:= TransformationBracket[N1, M1, N2, M2, n1, m1, n2, m2, L1, L2, l1, l2]

(*
<<WavefunctionIntegration.m

EvaluateBracket[N1_, M1_, N2_, M2_, n1_, m1_, n2_, m2_, L1_: Sqrt[2], L2_: 1/Sqrt[2], l1_: 1, l2_ : 1] := 
1/(2 Pi)^4 * NIntegrate[r1* r2
	* PolarWavefunction2D[N1, -M1, 
		Sqrt[r1^2 + r2^2 + 2 r1 r2 Cos[\[CurlyPhi]1 - \[CurlyPhi]2]], 
		Arg[r1 E^ (I \[CurlyPhi]1) + r2 E^(I \[CurlyPhi]2)], L1] 
	* PolarWavefunction2D[N2, -M2, 
		1/2 * Sqrt[r1^2 + r2^2 - 2 r1 r2 Cos[\[CurlyPhi]1 - \[CurlyPhi]2]], 
     		Arg[r1 E^(I \[CurlyPhi]1) - r2 E^(I \[CurlyPhi]2)], L2]
	* PolarWavefunction2D[n1, m1, r1, \[CurlyPhi]1, l1] 
	* PolarWavefunction2D[n2, m2, r2, \[CurlyPhi]2, l2],
	 {r1, 0, Max[2 N1 + Abs[M1] + 2 N2 + Abs[M2]]}, {\[CurlyPhi]1, 0, 2 Pi}, {r2, 0, Max[2 N1 + Abs[M1] + 2 N2 + Abs[M2]]}, {\[CurlyPhi]2, 0, 2 Pi}
		];

CheckBracket[maxpoints_:100, max_:10, sign_:0, zero_:0] := 
Module[{data,t1,t2,n1,m1,n2,m2,N1,M1,N2,M2},

data = Table[{n1 = RandomInteger[{0, max}], m1 = RandomInteger[{max*sign, max}], 
    n2 = RandomInteger[{0, max}], m2 = RandomInteger[{sign*max, max}],
    N1 = RandomInteger[{0, n1 + n2}], M1 = RandomInteger[{sign*max, max}], 
    n1 + n2 - N1 + RandomInteger[{0, zero}], 
    m1 + m2 - M1 + RandomInteger[{0, zero}]}, {maxpoints}];

t1 = Timing[TransformationBracket @@@ data];
t2 = Timing[EvaluateBracket @@@ data];

Print[{"timing the x'n brackets:" <> ToString[t1[[1]]], "timing the direct evaluation of x'n brackets:"<>ToString[t2[[1]]], "number of quanta that fail: " <>  Length[Select[Abs[t1[[2]] - t2[[2]]], # > 10^-3 || ! NumberQ[#] &]], 
  "number of nonzero quanta: " <> Length[Select[t1[[2]], # != 0 || ! NumberQ[#] &]] // N}];
]
*)

TM[N_, M_, n_, m_, n1_, m1_, n2_, m2_, x1_, x2_] := 
 TransformationBracket[N, M, n, m, n1, m1, n2, m2, Sqrt[x1 + x2], 
  Sqrt[(x1 x2)/(x1 + x2)], Sqrt[x1], Sqrt[x2]]
  
TMC[N_, M_, n_, m_, n1_, m1_, n2_, m2_, x1_, x2_] :=
 
 Module[{sin\[Delta] = Sqrt[x2/(x1 + x2)], 
   cos\[Delta] = Sqrt[x1/(x1 + x2)], 
   tan\[Delta] = Sqrt[x2/x1], \[Rho]},
  If[
   2 N + Abs[M] + 2 n + Abs[m] == 2 n1 + Abs[m1] + 2*n2 + Abs[m2] && 
    M + m == m1 + m2
   ,
   (-1)^(N + n + n1 + 
       n2) Sqrt[((n1 + Abs[m1])! n1! (n2 + Abs[m2])! n2!)/((N + 
        Abs[M])! N! (n + Abs[m])! n!)]*cos\[Delta]^(2 n1 + Abs[m1])*
    sin\[Delta]^(2 n2 + Abs[m2])*Sum[
     \[Rho] = 
      n1 - n2 - \[Gamma] + (Abs[m1] - Abs[m2] - \[Alpha] - \[Beta])/
        2;
     Sum[
      If[
       (* Mod[n+\[Rho]-\[Rho]3,2]==0&&*)
       (*Mod[
       n-\[Rho]-\[Rho]3,2]==0&&*)
       (*2(\[Rho]3+\[Gamma]3)>=Abs[
       Sign[m]*\[Alpha]+Sign[M]*\[Beta]+m2-m1] &&*)
       
       Mod[2 (\[Rho]3 + \[Gamma]3) + Sign[m] \[Alpha] + 
          Sign[M]*\[Beta] + m2 - m1, 4] == 0,
       (* If[!IntegerQ[(2(\[Rho]3+\[Gamma]3)+Sign[m]*\[Alpha]+Sign[
       M]*\[Beta]+m2-m1)/4],
       Print[{\[Gamma]3,\[Rho]3,\[Beta],\[Alpha],(2(\[Rho]3+\[Gamma]3)\
+Sign[m]*\[Alpha]+Sign[M]*\[Beta]+m2-m1)/
       4}]]; *)
       (-1)^(\[Rho]3 + (Abs[m] - \[Alpha])/2)*
        tan\[Delta]^(\[Alpha] + 2 \[Rho]) 
        *Multinomial[(N - \[Gamma]3 + \[Gamma])/
          2, (N - \[Gamma]3 - \[Gamma])/2, \[Gamma]3]
        *Multinomial[(n - \[Rho]3 + \[Rho])/2, (n - \[Rho]3 - \[Rho])/
          2, \[Rho]3]
        *Binomial[Abs[M], (Abs[M] + \[Beta])/2]
        *Binomial[Abs[m], (Abs[m] + \[Alpha])/2]
        *Binomial[\[Rho]3 + \[Gamma]3, (2 (\[Rho]3 + \[Gamma]3) + 
            Sign[m]*\[Alpha] + Sign[M]*\[Beta] + m2 - m1)/4]
       ,
       0
       ]
      , {\[Rho]3, 
       Max[0, Abs[Sign[m]*\[Alpha] + Sign[M]*\[Beta] + m2 - m1]/
           2 - \[Gamma]3] + 
        Mod[Max[0, 
           Abs[Sign[m]*\[Alpha] + Sign[M]*\[Beta] + m2 - m1]/
             2 - \[Gamma]3] + n + \[Rho], 2], 
       Min[n, \[Rho] + n, n - \[Rho]], 2}]
     , {\[Alpha], -Abs[m], Abs[m], 2}, {\[Beta], -Abs[M], Abs[M], 
      2}, {\[Gamma]3, 
      Max[Round@(Abs[Sign[m]*\[Alpha] + Sign[M]*\[Beta] + m2 - m1]/2 -
           n), 0], N}, {\[Gamma], \[Gamma]3 - N, N - \[Gamma]3, 2}
     ]
   ,
   0]
  ]

vars[nmax_, x1_, x2_] := Module[{n1, n2, m1, m2, N, M, n, m},
  Part[Reap[
    Do[
     n2 = (nmax - Abs[m1] - Abs[m2])/2 - n1;
     M = m1 + m2 - m;
     N = (nmax - Abs[m] - Abs[m1 + m2 - m])/2 - n;
     If[
      IntegerQ[n2] && n2 >= 0 && IntegerQ[N] && N >= 0
      ,
      Sow[{N, M, n, m, n1, m1, n2, m2, x1, x2}]
      ]
     , {m1, -nmax, nmax}
     , {m2, -(nmax - Abs[m1]), nmax - Abs[m1]}
     , {n1, 0, Floor[(nmax - Abs[m1] - Abs[m2])/2]}
     , {m, -nmax, nmax}
     , {n, 0, (nmax - Abs[m] - Abs[m1 + m2 - m])}
     ]
    ], 2, 1]
  ]
  
  checkTM[nmax_,x1_,x2_]:=Module[ {v, t1, t2, t0, tm1, tm2}, 
{t0, v} = Timing[vars[nmax, x1, x2]];
{t1, tm1} = Timing[TM @@@ v];
{t2, tm2} = Timing[TMC @@@ v];
Print[{t0, t1, t2, Length[v], Total[Abs[tm1 - tm2]], Total@(Abs[tm1] + Abs[tm2])/2}];
 ]
 

checkvars[nmax_, x1_, x2_] := Module[{v = 0, n1, n2, m1, m2, N, M, n, m, y1, y2},
	y1 = SetPrecision[x1, 16];
	y2 = SetPrecision[x2, 16];
  
    Do[
     n2 = (nmax - Abs[m1] - Abs[m2])/2 - n1;
     M = m1 + m2 - m;
     N = (nmax - Abs[m] - Abs[m1 + m2 - m])/2 - n;
     If[
      IntegerQ[n2] && n2 >= 0 && IntegerQ[N] && N >= 0
      ,
	  v = Max[ v, Abs[TM[N, M, n, m, n1, m1, n2, m2, x1, x2] 
	               - TMC[N, M, n, m, n1, m1, n2, m2, y1, y2]] ]
      ]
     , {m1, -nmax, nmax}
     , {m2, -(nmax - Abs[m1]), nmax - Abs[m1]}
     , {n1, 0, Floor[(nmax - Abs[m1] - Abs[m2])/2]}
     , {m, -nmax, nmax}
     , {n, 0, (nmax - Abs[m] - Abs[m1 + m2 - m])}
     ];
	 v
] 

End[ ]

EndPackage[ ]
