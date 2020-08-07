(* ::Package:: *)

(*gg\[Rule]Hh Cross Section*)
Jacobi[shat_, pt_, M1_, M2_] := (2*pt*shat)/
  Sqrt[-4*shat*(M1^2 + pt^2) + (M1^2 - M2^2 + shat)^2];
Prefactor = (GF^2*alphas^2)/(2^8*(2*Pi)^3);


(*Total differential cross section*)
dsigmaI[mhH_, pt_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  Prefactor*Jacobi[mhH^2, pt, M1, M2]*
   GeV2tofb*(Abs[
      CtriI[mhH^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhH^2] + 
       CboxI[\[Alpha], \[Beta]]*
        FboxTop[mhH^2, 
         M1^2/2 + M2^2/2 - mhH^2/2 - 
          1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
         M1^2/2 + M2^2/2 - mhH^2/2 + 
          1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
         M1, M2]]^2 + 
     Abs[CboxI[\[Alpha], \[Beta]]*
       GboxTop[mhH^2, 
        M1^2/2 + M2^2/2 - mhH^2/2 - 
         1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
        M1^2/2 + M2^2/2 - mhH^2/2 + 
         1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
        M1, M2]]^2);
dsigmaII[mhH_, pt_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  Prefactor*Jacobi[mhH^2, pt, M1, M2]*
   GeV2tofb*(Abs[
      CtriII[mhH^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhH^2] + 
       CboxII[\[Alpha], \[Beta]]*
        FboxTop[mhH^2, 
         M1^2/2 + M2^2/2 - mhH^2/2 - 
          1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
         M1^2/2 + M2^2/2 - mhH^2/2 + 
          1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
         M1, M2]]^2 + 
     Abs[CboxII[\[Alpha], \[Beta]]*
       GboxTop[mhH^2, 
        M1^2/2 + M2^2/2 - mhH^2/2 - 
         1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
        M1^2/2 + M2^2/2 - mhH^2/2 + 
         1/2 Sqrt[-4*mhH^2*(M1^2 + pt^2) + (M1^2 - M2^2 + mhH^2)^2], 
        M1, M2]]^2);


(*The PDF*)
ih = 0;
ggPDFintCenter[mhH_, S_] := 
  NIntegrate[(xf[ih, x, mhH, 0]/x)*((S*x*xf[ih, mhH^2/(S*x), mhH, 0])/
     mhH^2)*(2*mhH)/(S*x), {x, mhH^2/S, 1}];


dsigmadMhHLHC14I[mhH_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  ggPDFintCenter[mhH, 14000^2]*
   Re[NIntegrate[
     dsigmaI[mhH, pT, M1, M2, m12, \[Alpha], \[Beta]], {pT, 0, Sqrt[
      mhH^2/4 - M1^2/2 - M2^2/2 + (M1^2 - M2^2)^2/(4*mhH^2)]}]];
dsigmadMhHLHC14II[mhH_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  ggPDFintCenter[mhH, 14000^2]*
   Re[NIntegrate[
     dsigmaII[mhH, pT, M1, M2, m12, \[Alpha], \[Beta]], {pT, 0, Sqrt[
      mhH^2/4 - M1^2/2 - M2^2/2 + (M1^2 - M2^2)^2/(4*mhH^2)]}]];


ggtohHLHC14I[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  NIntegrate[
   dsigmadMhHLHC14I[mhH, M1, M2, m12, \[Alpha], \[Beta]], {mhH, 
    M1 + M2 + 0.1, 2000 + 0.1}];
ggtohHLHC14II[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  NIntegrate[
   dsigmadMhHLHC14II[mhH, M1, M2, m12, \[Alpha], \[Beta]], {mhH, 
    M1 + M2 + 0.1, 2000 + 0.1}];
