(* ::Package:: *)

(*gg\[Rule]hh Cross Section*)
Jacobi[shat_,pt_,M1_]:=(2*pt*shat)/Sqrt[shat (-4*( M1^2+ pt^2)+shat)];
Prefactor=(GF^2*alphas^2)/(2^8*(2*Pi)^3)*1/2;


(*Total differential cross section*)
dsigmaI[mhh_, pt_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  Prefactor*Jacobi[mhh^2, pt, M1]*
   GeV2tofb*(Abs[
      CtriI[mhh^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhh^2] + 
       CboxI[\[Alpha], \[Beta]]*
        FboxTop[mhh^2, 
         M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1]]^2 + 
     Abs[CboxI[\[Alpha], \[Beta]]*
       GboxTop[mhh^2, 
        M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1]]^2);
dsigmaII[mhh_, pt_, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  Prefactor*Jacobi[mhh^2, pt, M1]*
   GeV2tofb*(Abs[
      CtriII[mhh^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhh^2] + 
       CboxII[\[Alpha], \[Beta]]*
        FboxTop[mhh^2, 
         M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1]]^2 + 
     Abs[CboxII[\[Alpha], \[Beta]]*
       GboxTop[mhh^2, 
        M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1]]^2);


(*Prepare PDFs*)
(*We only consider the gluon PDF*)
ih = 0;
ggPDFintCenter[mhh_, S_] := 
  NIntegrate[(xf[ih, x, mhh, 0]/x)*((S*x*xf[ih, mhh^2/(S*x), mhh, 0])/
     mhh^2)*(2*mhh)/(S*x), {x, mhh^2/S, 1}];


(* p p differential cross section*)
dsigmadMhhLHC14I[mhh_,M1_,M2_,m12_,\[Alpha]_,\[Beta]_]:=ggPDFintCenter[mhh,14000^2]*Re[NIntegrate[dsigmaI[mhh,pT,M1,M2,m12,\[Alpha],\[Beta]],{pT,0,Sqrt[mhh^2/4-M1^2]}]];
dsigmadMhhLHC14II[mhh_,M1_,M2_,m12_,\[Alpha]_,\[Beta]_]:=ggPDFintCenter[mhh,14000^2]*Re[NIntegrate[dsigmaII[mhh,pT,M1,M2,m12,\[Alpha],\[Beta]],{pT,0,Sqrt[mhh^2/4-M1^2]}]];


(* Total p p cross section*)
ggtohhLHC14I[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  NIntegrate[
   dsigmadMhhLHC14I[mhh, M1, M2, m12, \[Alpha], \[Beta]], {mhh, 
    2*M1 + 0.1, 2000 + 0.1}];
ggtohhLHC14II[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  NIntegrate[
   dsigmadMhhLHC14II[mhh, M1, M2, m12, \[Alpha], \[Beta]], {mhh, 
    2*M1 + 0.1, 2000 + 0.1}];


dsigmaIN[mhh_?NumericQ, pt_?NumericQ, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := Block[{tmp},
  tmp=Prefactor*Jacobi[mhh^2, pt, M1]*
   GeV2tofb*(Abs[
      CtriI[mhh^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhh^2] + 
       CboxI[\[Alpha], \[Beta]]*
        FboxTop[mhh^2, 
         M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1]]^2 + 
     Abs[CboxI[\[Alpha], \[Beta]]*
       GboxTop[mhh^2, 
        M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1]]^2);
    ClearCache[];
    tmp
];
dsigmaIIN[mhh_?NumericQ, pt_?NumericQ, M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := Block[{tmp},
tmp=Prefactor*Jacobi[mhh^2, pt, M1]*
   GeV2tofb*(Abs[
      CtriII[mhh^2, M1, M2, m12, \[Alpha], \[Beta]]*FtriTop[mhh^2] + 
       CboxII[\[Alpha], \[Beta]]*
        FboxTop[mhh^2, 
         M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
         M1]]^2 + 
     Abs[CboxII[\[Alpha], \[Beta]]*
       GboxTop[mhh^2, 
        M1^2 - mhh^2/2 - 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1^2 - mhh^2/2 + 1/2 Sqrt[mhh^2*(mhh^2 - 4*(M1^2 + pt^2))], 
        M1]]^2);ClearCache[];
 tmp
];


ggtohhLHC14Itogether[M1_, M2_, m12_, \[Alpha]_, \[Beta]_,S_:14000^2]:=Block[{r,mhhMax=2000.0},
r=ImplicitRegion[
 mhh0 > (2*M1+1.0)/mhhMax && mhh0 < 1.0 && pT0 > 0 && 
  pT0 < Sqrt[mhh0^2/4 - (M1/mhhMax)^2] && x0 > mhh0^2*mhhMax^2/S && x0 <= 1, {x0, mhh0, 
  pT0}];
  NIntegrate[
 Re[dsigmaIN[mhhS*mhhMax, pTS*mhhMax, M1, M2, 
    m12, \[Alpha], \[Beta]]]*(2 xf[ih, xS, mhhS*mhhMax, 0] xf[ih, mhhS^2*mhhMax^2/(S*xS), mhhS*mhhMax, 0])/(mhhS*mhhMax*xS), {xS, mhhS, pTS}\[Element]r,PrecisionGoal->2,WorkingPrecision->5]
]
