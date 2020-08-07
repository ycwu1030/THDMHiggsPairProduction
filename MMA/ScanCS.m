(* ::Package:: *)

Install["LoopTools"];

pdfpath = "/Users/ycwu/Workingspace/Utilities/MSTW2008/";
SetDirectory[pdfpath];
<< mstwpdf.m
prefix=pdfpath<>"Grids/mstw2008nlo";
ReadPDFGrid[prefix,0];


type=ToExpression[$ScriptCommandLine[[2]]];
proc=ToExpression[$ScriptCommandLine[[3]]];
MHL=ToExpression[$ScriptCommandLine[[4]]];
PROCSRCNAME={{"VggHLHL`","CSggHLHL`"},{"VggHHHL`","CSggHHHL`"}};


SourceDir="/Users/ycwu/Desktop/TempFiles/2HDM_HiggsPair/MMA";
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


str=OpenWrite["CSdata.dat"];
tbList={0.5,1.0,5.0,10.0};
For[itb=1,itb<=Length[tbList],itb++,
	tb=tbList[[itb]];
	For[m12=0,m12<=200,m12+=10,
	Print["Calculating CS for: ","MHL=",MHL," tb=",tb," m12=",m12];
	cs=Func2BeRun[MHL,125.0,m12,ArcTan[tb],ArcTan[tb]]//Quiet;
	Write[str,StringRiffle[ToString/@{MHL,tb,m12,cs},{"","\t","\n"}]];
	]
]
