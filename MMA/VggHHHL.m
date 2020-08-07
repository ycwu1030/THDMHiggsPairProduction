(* ::Package:: *)

(*Define variables for gg\[Rule]Hh production calculation*)
(*From Appendix A of hep-ph/9603205*)


beta[M1_, shat_] := Sqrt[1 - (4*M1^2)/shat];
\[Rho]t[M1_] := M1^2/mt^2;


(*Triangle and Box functions*)
FtriTop[shat_] := 
  2/(shat/mt^2)*(2 + (4 - shat/mt^2)*mt^2*
      C0[0, 0, shat, mt^2, mt^2, mt^2]);

FboxTop[shat_, that_, uhat_, M1_, M2_] := 
  1/(shat/mt^2)^2*(4*shat/mt^2 + 
     8*shat/mt^2*mt^2*C0[0, 0, shat, mt^2, mt^2, mt^2] - 
     2*shat/mt^2*(shat/mt^2 + M1^2/mt^2 + M2^2/mt^2 - 8)*
      mt^4*(D0[0, 0, M1^2, M2^2, shat, uhat, mt^2, mt^2, mt^2, mt^2] +
         D0[0, 0, M1^2, M2^2, shat, that, mt^2, mt^2, mt^2, mt^2] + 
        D0[0, M1^2, 0, M2^2, that, uhat, mt^2, mt^2, mt^2, 
         mt^2]) + (M1^2/mt^2 + M2^2/mt^2 - 8)*
      mt^2*((that - M1^2)/mt^2*C0[0, M1^2, that, mt^2, mt^2, mt^2] + (
         that - M2^2)/mt^2*C0[0, M2^2, that, mt^2, mt^2, mt^2] + (
         uhat - M1^2)/mt^2*C0[0, M1^2, uhat, mt^2, mt^2, mt^2] + (
         uhat - M2^2)/mt^2*
         C0[0, M2^2, uhat, mt^2, mt^2, 
          mt^2] - (that/mt^2*uhat/mt^2 - (M1^2/mt^2)*(M2^2/mt^2))*
         mt^2*D0[0, M1^2, 0, M2^2, that, uhat, mt^2, mt^2, mt^2, 
          mt^2])); 

GboxTop[shat_, that_, uhat_, M1_, M2_] := 
  1/((shat/mt^2)*((that/mt^2)*(uhat/mt^2) - (M1^2/mt^2)*(M2^2/
        mt^2)))*(((that/mt^2)^2 + (M1^2/mt^2)*(M2^2/mt^2) - 
        8*that/mt^2)*
      mt^2*(shat/mt^2* C0[0, 0, shat, mt^2, mt^2, mt^2] + (
         that - M1^2)/mt^2*C0[0, M1^2, that, mt^2, mt^2, mt^2] + (
         that - M2^2)/mt^2*C0[0, M2^2, that, mt^2, mt^2, mt^2] - 
        shat/mt^2*that/mt^2* mt^2*
         D0[0, 0, M1^2, M2^2, shat, that, mt^2, mt^2, mt^2, 
          mt^2]) + ((uhat/mt^2)^2 + (M1^2/mt^2)*(M2^2/mt^2) - 
        8*uhat/mt^2)*
      mt^2*(shat/mt^2* C0[0, 0, shat, mt^2, mt^2, mt^2] + (
         uhat - M1^2)/mt^2*C0[0, M1^2, uhat, mt^2, mt^2, mt^2] + (
         uhat - M2^2)/mt^2*C0[0, M2^2, uhat, mt^2, mt^2, mt^2] - 
        shat/mt^2*uhat/mt^2* mt^2*
         D0[0, 0, M1^2, M2^2, shat, uhat, mt^2, mt^2, mt^2, 
          mt^2]) - ((that/mt^2)^2 + (uhat/mt^2)^2 - 
        2*(M1^2/mt^2)*(M2^2/mt^2))*(that/mt^2 + uhat/mt^2 - 8)*mt^2*
      C0[M1^2, M2^2, shat, mt^2, mt^2, mt^2] - 
     2*(that/mt^2 + uhat/mt^2 - 
        8)*(that/mt^2*uhat/mt^2 - (M1^2/mt^2)*(M2^2/mt^2))*
      mt^4*(D0[0, 0, M1^2, M2^2, shat, that, mt^2, mt^2, mt^2, mt^2] +
         D0[0, 0, M1^2, M2^2, shat, uhat, mt^2, mt^2, mt^2, mt^2] + 
        D0[0, M1^2, 0, M2^2, that, uhat, mt^2, mt^2, mt^2, mt^2]));


(* box coefficients *)
(*CboxI[\[Alpha]_,\[Beta]_]:=ghuuI[\[Alpha],\
\[Beta]]^2;*)

CboxI[\[Alpha]_, \[Beta]_] := 
  ghuuI[\[Alpha], \[Beta]]*gHuuI[\[Alpha], \[Beta]];
(* triangle coefficients *)

CtriI[shat_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := (
   ghuuI[\[Alpha], \[Beta]]*vev*
    ghhH[M1, M2, m12, \[Alpha], \[Beta]])/(
   shat - M1^2 + I*M1*hWidI[M1, \[Alpha], \[Beta]]) + (
   gHuuI[\[Alpha], \[Beta]]*vev*
    ghHH[M1, M2, m12, \[Alpha], \[Beta]])/(
   shat - M2^2 + I*M2*HWidI[M2, M1, m12, \[Alpha], \[Beta]]);
(*CtriIres[shat_,M1_,M2_,m12_,\[Alpha]_,\[Beta]_]:=(gHuuI[\[Alpha],\
\[Beta]]*vev*ghhH[M1,M2,m12,\[Alpha],\[Beta]])/(shat-M2^2+I*M2*HWidI[\
M2,M1,m12,\[Alpha],\[Beta]]);*)


(* box coefficients *)
(*CboxII[\[Alpha]_,\[Beta]_]:=ghuuII[\[Alpha],\
\[Beta]]^2;*)

CboxII[\[Alpha]_, \[Beta]_] := 
  ghuuII[\[Alpha], \[Beta]]*gHuuII[\[Alpha], \[Beta]];
(* triangle coefficients *)

CtriII[shat_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := (
   ghuuII[\[Alpha], \[Beta]]*vev*
    ghhH[M1, M2, m12, \[Alpha], \[Beta]])/(
   shat - M1^2 + I*M1*hWidII[M1, \[Alpha], \[Beta]]) + (
   gHuuII[\[Alpha], \[Beta]]*vev*
    ghHH[M1, M2, m12, \[Alpha], \[Beta]])/(
   shat - M2^2 + I*M2*HWidII[M2, M1, m12, \[Alpha], \[Beta]]);
(*CtriIIres[shat_,M1_,M2_,m12_,\[Alpha]_,\[Beta]_]:=(gHuuII[\[Alpha],\
\[Beta]]*vev*ghhH[M1,M2,m12,\[Alpha],\[Beta]])/(shat-M2^2+I*M2*HWidII[\
M2,M1,m12,\[Alpha],\[Beta]]);*)
