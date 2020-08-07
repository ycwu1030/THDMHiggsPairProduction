(* ::Package:: *)

(*2HDM Couplings*)


(*Loop integral functions*)
f[x_] := If[x <= 1, 
   ArcSin[x^(
    1/2)]^2, -(Log[(1 + Sqrt[1 - 1/x])/(1 - Sqrt[1 - 1/x])] - I*Pi)^2/
    4];
Hfermion[x_] := (2*(x + (x - 1)*f[x]))/x^2;
Hvector[x_] := -((2*x^2 + 3*x + 3*(2*x - 1)*f[x])/x^2);
Afermion[x_] := (2*f[x])/x^2;


(*Couplings in Type-I 2HDM*)
ghuuI[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] + Cos[\[Beta] - \[Alpha]]/Tan[\[Beta]];
ghddI[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] + Cos[\[Beta] - \[Alpha]]/Tan[\[Beta]];
ghllI[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] + Cos[\[Beta] - \[Alpha]]/Tan[\[Beta]];
ghVV[\[Alpha]_, \[Beta]_] := Sin[\[Beta] - \[Alpha]];

gHuuI[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] - Sin[\[Beta] - \[Alpha]]/Tan[\[Beta]];
gHddI[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] - Sin[\[Beta] - \[Alpha]]/Tan[\[Beta]];
gHllI[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] - Sin[\[Beta] - \[Alpha]]/Tan[\[Beta]];

gHVV[\[Alpha]_, \[Beta]_] := Cos[\[Beta] - \[Alpha]];

gAuuI[\[Alpha]_, \[Beta]_] := 1/Tan[\[Beta]];
gAddI[\[Alpha]_, \[Beta]_] := Tan[\[Beta]];
gAllI[\[Alpha]_, \[Beta]_] := Tan[\[Beta]];


(*Couplings in Type-II 2HDM*)
ghuuII[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] + Cos[\[Beta] - \[Alpha]]/Tan[\[Beta]];
ghddII[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] - Cos[\[Beta] - \[Alpha]]*Tan[\[Beta]];
ghllII[\[Alpha]_, \[Beta]_] := 
  Sin[\[Beta] - \[Alpha]] - Cos[\[Beta] - \[Alpha]]*Tan[\[Beta]];
ghVV[\[Alpha]_, \[Beta]_] := Sin[\[Beta] - \[Alpha]];

gHuuII[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] - Sin[\[Beta] - \[Alpha]]/Tan[\[Beta]];
gHddII[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] + Sin[\[Beta] - \[Alpha]]*Tan[\[Beta]];
gHllII[\[Alpha]_, \[Beta]_] := 
  Cos[\[Beta] - \[Alpha]] + Sin[\[Beta] - \[Alpha]]*Tan[\[Beta]];

gHVV[\[Alpha]_, \[Beta]_] := Cos[\[Beta] - \[Alpha]];

gAuuII[\[Alpha]_, \[Beta]_] := 1/Tan[\[Beta]];
gAddII[\[Alpha]_, \[Beta]_] := Tan[\[Beta]];
gAllII[\[Alpha]_, \[Beta]_] := Tan[\[Beta]];


(*Cubic Couplings, following 1407.0281*)
ghhh[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := 
  3/vev*(M1^2*(Sin[\[Beta] - \[Alpha]] + 
        2*Sin[\[Beta] - \[Alpha]]*Cos[\[Beta] - \[Alpha]]^2 + (
        2*Cos[\[Beta] - \[Alpha]]^3)/Tan[2*\[Beta]]) - (2*m12^2)/(
      Sin[\[Beta]]*
       Cos[\[Beta]])*(Sin[\[Beta] - \[Alpha]] + 
        Cos[\[Beta] - \[Alpha]]/Tan[2*\[Beta]])*
      Cos[\[Beta] - \[Alpha]]^2);

ghhH[M1_, M2_, 
   m12_, \[Alpha]_, \[Beta]_] := -(Cos[\[Beta] - \[Alpha]]/
    vev)*((M2^2 + 2*M1^2)*(1 + (
        2*Sin[\[Beta] - \[Alpha]]*Cos[\[Beta] - \[Alpha]])/
        Tan[2*\[Beta]] - 2*Cos[\[Beta] - \[Alpha]]^2) - (2*m12^2)/(
      Sin[\[Beta]]*
       Cos[\[Beta]])*(2 + (
        3*Sin[\[Beta] - \[Alpha]]*Cos[\[Beta] - \[Alpha]])/
        Tan[2*\[Beta]] - 3*Cos[\[Beta] - \[Alpha]]^2));

ghHH[M1_, M2_, m12_, \[Alpha]_, \[Beta]_] := -Sin[\[Beta] - \[Alpha]]/
   vev*((M1^2 + 2*M2^2)*Sin[2*\[Alpha]]/Sin[2*\[Beta]] - (2*m12^2)/(
      Sin[\[Beta]]*
       Cos[\[Beta]])*(2 - (
        3*Sin[\[Beta] - \[Alpha]]*Cos[\[Beta] - \[Alpha]])/
        Tan[2*\[Beta]] - 3*Sin[\[Beta] - \[Alpha]]^2));


(*Tree-level estimation of the h/H total width*)

(*For h*)
hWidI[M1_, \[Alpha]_, \[Beta]_] := (Nc*GF*mb^2)/(4 Sqrt[2]*Pi)*
    ghddI[\[Alpha], \[Beta]]^2*M1*(1 - (4*mb^2)/M1^2)^(3/2) + (
    Nc*GF*mc^2)/(4 Sqrt[2]*Pi)*ghuuI[\[Alpha], \[Beta]]^2*
    M1*(1 - (4*mc^2)/M1^2)^(3/2) + (GF*mtau^2)/(4 Sqrt[2]*Pi)*
    ghllI[\[Alpha], \[Beta]]^2*M1*(1 - (4*mtau^2)/M1^2)^(3/2) + (
    GF*alphaEM^2*M1^3)/(128 Sqrt[2]*Pi^3)*
    Abs[Nc*Qu^2*ghuuI[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mt^2)] + 
      Nc*Qd^2*ghddI[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mb^2)] + 
      Ql^2*ghllI[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mtau^2)] + 
      ghVV[\[Alpha], \[Beta]]*Hvector[M1^2/(4*mw^2)]]^2 + (
    GF*alphas^2*M1^3)/(36 Sqrt[2]*Pi^3)*
    Abs[3/4*ghuuI[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mt^2)] + 
      3/4*ghddI[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mb^2)]]^2;
hWidII[M1_, \[Alpha]_, \[Beta]_] := (Nc*GF*mb^2)/(4 Sqrt[2]*Pi)*
    ghddII[\[Alpha], \[Beta]]^2*M1*(1 - (4*mb^2)/M1^2)^(3/2) + (
    Nc*GF*mc^2)/(4 Sqrt[2]*Pi)*ghuuII[\[Alpha], \[Beta]]^2*
    M1*(1 - (4*mc^2)/M1^2)^(3/2) + (GF*mtau^2)/(4 Sqrt[2]*Pi)*
    ghllII[\[Alpha], \[Beta]]^2*M1*(1 - (4*mtau^2)/M1^2)^(3/2) + (
    GF*alphaEM^2*M1^3)/(128 Sqrt[2]*Pi^3)*
    Abs[Nc*Qu^2*ghuuII[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mt^2)] + 
      Nc*Qd^2*ghddII[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mb^2)] + 
      Ql^2*ghllII[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mtau^2)] + 
      ghVV[\[Alpha], \[Beta]]*Hvector[M1^2/(4*mw^2)]]^2 + (
    GF*alphas^2*M1^3)/(36 Sqrt[2]*Pi^3)*
    Abs[3/4*ghuuII[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mt^2)] + 
      3/4*ghddII[\[Alpha], \[Beta]]*Hfermion[M1^2/(4*mb^2)]]^2;


(*For H*)
HWidI[M2_, M1_, 
   m12_, \[Alpha]_, \[Beta]_] := (Nc*GF*mb^2)/(4 Sqrt[2]*Pi)*
    gHddI[\[Alpha], \[Beta]]^2*M2*(1 - (4*mb^2)/M2^2)^(3/2) + (
    Nc*GF*mc^2)/(4 Sqrt[2]*Pi)*gHuuI[\[Alpha], \[Beta]]^2*
    M2*(1 - (4*mc^2)/M2^2)^(3/2) + (GF*mtau^2)/(4 Sqrt[2]*Pi)*
    gHllI[\[Alpha], \[Beta]]^2*M2*(1 - (4*mtau^2)/M2^2)^(3/2) + (
    GF*alphaEM^2*M2^3)/(128 Sqrt[2]*Pi^3)*
    Abs[Nc*Qu^2*gHuuI[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mt^2)] + 
      Nc*Qd^2*gHddI[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mb^2)] + 
      Ql^2*gHllI[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mtau^2)] + 
      gHVV[\[Alpha], \[Beta]]*Hvector[M2^2/(4*mw^2)]]^2 + (
    GF*alphas^2*M2^3)/(36 Sqrt[2]*Pi^3)*
    Abs[3/4*gHuuI[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mt^2)] + 
      3/4*gHddI[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mb^2)]]^2 + 
   Which[M2 <= M1 + 2*me, 0.,
    M1 + 2*me < M2 <= M1 + 2*ms, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, m12, \[Alpha], \[Beta]])^2*
     me^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*ms < M2 <= M1 + 2*mmu, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mmu < M2 <= M1 + 2*mc, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 
          1.)*(2. - 1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mc < M2 <= M1 + 2*mtau, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuI[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mtau < M2 <= M1 + 2*mb, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuI[\[Alpha], \[Beta]])^2 + 
       mtau^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 
          1.)*(2. - 1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mb < M2 <= 2*M1, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuI[\[Alpha], \[Beta]])^2 + 
       mtau^2/vev^2*(ghllI[\[Alpha], \[Beta]])^2 + (3*mb^2)/
         vev^2*(ghddI[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    2*M1 < M2, 
    1/(32*Pi*M2)*(ghhH[M1, M2, m12, \[Alpha], \[Beta]])^2*
     Sqrt[1 - (4*M1^2)/M2^2]
    ];
HWidII[M2_, M1_, 
   m12_, \[Alpha]_, \[Beta]_] := (Nc*GF*mb^2)/(4 Sqrt[2]*Pi)*
    gHddII[\[Alpha], \[Beta]]^2*M2*(1 - (4*mb^2)/M2^2)^(3/2) + (
    Nc*GF*mc^2)/(4 Sqrt[2]*Pi)*gHuuII[\[Alpha], \[Beta]]^2*
    M2*(1 - (4*mc^2)/M2^2)^(3/2) + (GF*mtau^2)/(4 Sqrt[2]*Pi)*
    gHllII[\[Alpha], \[Beta]]^2*M2*(1 - (4*mtau^2)/M2^2)^(3/2) + (
    GF*alphaEM^2*M2^3)/(128 Sqrt[2]*Pi^3)*
    Abs[Nc*Qu^2*gHuuII[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mt^2)] + 
      Nc*Qd^2*gHddII[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mb^2)] + 
      Ql^2*gHllII[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mtau^2)] + 
      gHVV[\[Alpha], \[Beta]]*Hvector[M2^2/(4*mw^2)]]^2 + (
    GF*alphas^2*M2^3)/(36 Sqrt[2]*Pi^3)*
    Abs[3/4*gHuuII[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mt^2)] + 
      3/4*gHddII[\[Alpha], \[Beta]]*Hfermion[M2^2/(4*mb^2)]]^2 + 
   Which[M2 <= M1 + 2*me, 0.,
    M1 + 2*me < M2 <= M1 + 2*ms, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, m12, \[Alpha], \[Beta]])^2*
     me^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*ms < M2 <= M1 + 2*mmu, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mmu < M2 <= M1 + 2*mc, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 
          1.)*(2. - 1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mc < M2 <= M1 + 2*mtau, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuII[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mtau < M2 <= M1 + 2*mb, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuII[\[Alpha], \[Beta]])^2 + 
       mtau^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 
          1.)*(2. - 1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    M1 + 2*mb < M2 <= 2*M1, 
    1/(32*Pi^3*M2)*(ghhH[M1, M2, 
        m12, \[Alpha], \[Beta]])^2*(me^2/
         vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*ms^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2 + 
       mmu^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*mc^2)/
         vev^2*(ghuuII[\[Alpha], \[Beta]])^2 + 
       mtau^2/vev^2*(ghllII[\[Alpha], \[Beta]])^2 + (3*mb^2)/
         vev^2*(ghddII[\[Alpha], \[Beta]])^2)*((M1^2/M2^2 - 1.)*(2. - 
          1/2*Log[M1^2/M2^2]) + (1 - 5*M1^2/M2^2)/
         Sqrt[4*M1^2/M2^2 - 
           1]*(ArcTan[(2*M1^2/M2^2 - 1)/Sqrt[4*M1^2/M2^2 - 1]] - 
          ArcTan[1/Sqrt[4*M1^2/M2^2 - 1]])),
    2*M1 < M2, 
    1/(32*Pi*M2)*(ghhH[M1, M2, m12, \[Alpha], \[Beta]])^2*
     Sqrt[1 - (4*M1^2)/M2^2]
    ];
