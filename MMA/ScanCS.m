(* ::Package:: *)

Install["LoopTools"];

pdfpath = "";
SetDirectory[pdfpath];
<< mstwpdf.m
prefix=pdfpath<>"Grids/mstw2008nlo";
ReadPDFGrid[prefix,0];


type=ToExpression[$ScriptCommandLine[[2]]];
proc=ToExpression[$ScriptCommandLine[[3]]];
MHL=ToExpression[$ScriptCommandLine[[4]]];
TYPENAME={"I","II"};
PROCNAME={"HLHL","HHHL"};
PROCSRCNAME={{"VggHLHL`","CSggHLHL`"},{"VggHHHL`","CSggHHHL`"}};


SourceDir="";
SetDirectory[SourceDir];
Get["PhysicsConstants`"];
Get["THDMCouplings`"];
Get[PROCSRCNAME[[proc]][[1]]];
Get[PROCSRCNAME[[proc]][[2]]];


FuncPtr={{ggtohhLHC14I, ggtohhLHC14II}, {ggtohHLHC14I, ggtohHLHC14II}};
Func2BeRun=FuncPtr[[proc]][[type]];
Print[ToString[Func2BeRun]];
Print["Type: ",type];
Print["Process:",proc];


str=OpenWrite["THDM_"<>TYPENAME[[type]]<>"_gg2"<>PROCNAME[[proc]]<>"_CS_MHL"<>ToString[MHL]<>".dat"];
WriteString[str,"MHL  m12  tb  cs\n"];
tbList={0.5,1.0,5.0,10.0};
For[itb=1,itb<=Length[tbList],itb++,
	tb=tbList[[itb]];
	For[m12=0,m12<=200,m12+=10,
	Print["Calculating CS for: ","MHL=",MHL," tb=",tb," m12=",m12];
	cs=Func2BeRun[MHL,125.0,m12,ArcTan[tb],ArcTan[tb]]//Quiet;
	tmp=ToString[MHL]<>"  "<>ToString[m12]<>"  "<>ToString[tb]<>"  "<>ToString[cs]<>"\n";
	WriteString[str,tmp];
	]
]
Close[str];
