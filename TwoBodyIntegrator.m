(*
  This module has been verified to work with FeynCalc 10.
  Please feel free to use it with no accreditation.
  The formulae can be checked against "From Spinors to Supersymmetry" Appendix D,
  coauthored by Howard Haber.
*)

(* Load FeynCalc - uncomment if necessary *)
<< FeynCalc`

(* Pre-Defined Phase Space Integrations *)

Kallen\[Lambda][aa_, bb_, cc_] = aa^2 + bb^2 + cc^2 - 2 (aa bb + bb cc + cc aa);

twobody$R2[xsq_, m1_ : 0, m2_ : 0] = \[Pi]/(2 xsq) \[Sqrt]Kallen\[Lambda][xsq, m1^2, m2^2] // Simplify;
twobody$T1q\[Mu][xsq_, m1_ : 0, m2_ : 0] = ((xsq + m1^2 - m2^2)/(2 xsq)) twobody$R2[xsq, m1, m2] // Simplify;
twobody$T2q\[Mu][xsq_, m1_ : 0, m2_ : 0] = twobody$T1q\[Mu][xsq, m2, m1];
twobody$T11metric[xsq_, m1_ : 0, m2_ : 0] = -(1/(12 xsq^2)) (xsq Kallen\[Lambda][xsq, m1^2, m2^2]) twobody$R2[xsq, m1, m2] // Simplify;
twobody$T11q\[Mu]q\[Nu][xsq_, m1_ : 0, m2_ : 0] = 1/(12 xsq^2) (4 Kallen\[Lambda][xsq, m1^2, m2^2] + 12 m1^2 xsq) twobody$R2[xsq, m1, m2] // Simplify;
twobody$T22metric[xsq_, m1_ : 0, m2_ : 0] = twobody$T11metric[xsq, m1^2, m2^2];
twobody$T22q\[Mu]q\[Nu][xsq_, m1_ : 0, m2_ : 0] = 1/(12 xsq^2) (4 Kallen\[Lambda][xsq, m1^2, m2^2] + 12 m2^2 xsq) twobody$R2[xsq, m1, m2] // Simplify;
twobody$T12metric[xsq_, m1_ : 0, m2_ : 0] = -twobody$T11metric[xsq, m1^2, m2^2];
twobody$T12q\[Mu]q\[Nu][xsq_, m1_ : 0, m2_ : 0] = -(1/(12 xsq^2)) (4 Kallen\[Lambda][xsq, m1^2, m2^2] - 6 xsq (xsq - m1^2 - m2^2)) twobody$R2[xsq, m1, m2];

(* Two-Body Phase Space Integration Module *)

Clear@TwoBodyIntegratorModule
TwoBodyIntegratorModule[MsquaredExpr_, m1_ : 0, m2_ : 0, mom1_ : p1, mom2_ : p2, fourvector_ : q] := 
 Module[{t0$twobody, t1$twobody, t2$twobody, t11$twobody, t12$twobody, t22$twobody, result$twobody, twobody$mu, twobody$nu},
  
  t0$twobody = (MsquaredExpr /. {p1 -> 0, p2 -> 0}) // Simplify;
  t1$twobody = (FourDivergence[MsquaredExpr, FV[mom1, twobody$mu]] /. {p1 -> 0, p2 -> 0});
  t2$twobody = (FourDivergence[MsquaredExpr, FV[mom2, twobody$mu]] /. {p1 -> 0, p2 -> 0});
  t11$twobody = 1/2 (FourDivergence[MsquaredExpr, FV[mom1, twobody$mu], FV[mom1, twobody$nu]] /. {p1 -> 0, p2 -> 0});
  t22$twobody = 1/2 (FourDivergence[MsquaredExpr, FV[mom2, twobody$mu], FV[mom2, twobody$nu]] /. {p1 -> 0, p2 -> 0});
  t12$twobody = (FourDivergence[MsquaredExpr, FV[mom1, twobody$mu], FV[mom2, twobody$nu]] /. {p1 -> 0, p2 -> 0});
  
  result$twobody = (t0$twobody*twobody$R2[SP[fourvector], m1, m2] + 
     Contract[t1$twobody*FV[fourvector, twobody$mu]]*
      twobody$T1q\[Mu][SP[fourvector], m1, m2] + 
     Contract[t2$twobody*FV[fourvector, twobody$mu]]*
      twobody$T2q\[Mu][SP[fourvector], m1, m2] +
     Contract[t12$twobody*MT[twobody$mu, twobody$nu]]*
      twobody$T12metric[SP[fourvector], m1, m2] +
     Contract[t11$twobody*MT[twobody$mu, twobody$nu]]*
      twobody$T11metric[SP[fourvector], m1, m2] +
     Contract[t22$twobody*MT[twobody$mu, twobody$nu]]*
      twobody$T22metric[SP[fourvector], m1, m2] +
     Contract[
       t12$twobody*FV[fourvector, twobody$mu]*FV[fourvector, twobody$nu]]*
      twobody$T12q\[Mu]q\[Nu][SP[fourvector], m1, m2] +
     Contract[
       t11$twobody*FV[fourvector, twobody$mu]*FV[fourvector, twobody$nu]]*
      twobody$T11q\[Mu]q\[Nu][SP[fourvector], m1, m2] +
     Contract[
       t22$twobody*FV[fourvector, twobody$mu]*FV[fourvector, twobody$nu]]*
      twobody$T22q\[Mu]q\[Nu][SP[fourvector], m1, m2]
    );
  result$twobody
  ];

(* Example usage:
   (This block is provided for illustration. You can remove or modify it as needed.)

test$MatrixElement = (SpinorVBar[p, mp] . (ΓV*GA[μ] + Γp*FV[p, μ] + Γk*FV[k, μ] + ΓPV*GA[μ] . GA5) . 
                       SpinorU[k, mk])*
                     (SpinorUBar[p1, m1] . (ΓVf*GA[μ] + Γp1*FV[p1, μ] + Γp2*FV[p2, μ] + ΓPVf*GA[μ] . GA5) . 
                       SpinorV[p2, m2]) +
                     (SpinorVBar[p, mp] . (ΓS + ΓPS*GA5) . SpinorU[k, mk])*
                     (SpinorUBar[p1, m1] . (ΓSf + ΓPSf*GA5) . SpinorV[p2, m2]) +
                     (SpinorVBar[p, mp] . (ΓT*(GA[μ, ν] - GA[ν, μ]) + ΓPT*(GA[μ, ν] - GA[ν, μ]) . GA5) . 
                       SpinorU[k, mk])*
                     (SpinorUBar[p1, m1] . (ΓTf*(GA[μ, ν] - GA[ν, μ]) + ΓPTf*(GA[μ, ν] - GA[ν, μ]) . GA5) . 
                       SpinorV[p2, m2])
                     // Contract // Simplify;
                     
test$structurelist = {ΓV, ΓS, ΓPV, ΓPS, Γp, Γk, ΓVf, ΓSf, ΓPVf, ΓPSf, Γp1, Γp2, ΓT, ΓPT, ΓTf, ΓPTf};

test$Msquared = 1/4 FermionSpinSum[test$MatrixElement*ComplexConjugate[test$MatrixElement]];
test$Msquared = test$Msquared // DiracSimplify // Simplify;

(* Remove extra terms for computational efficiency *)
test$Msquared = test$Msquared /. 
    Thread[{Γk, Γp2, ΓS, ΓSf, ΓPSf, ΓPS, ΓT, ΓPT, ΓTf, ΓPTf} -> 0] // Simplify;

(* Perform phase space integration *)
integratedResult = Simplify[TwoBodyIntegratorModule[test$Msquared], Assumptions -> SP[q] > 0];
*)
