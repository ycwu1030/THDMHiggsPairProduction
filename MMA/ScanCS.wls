#!/usr/bin/env wolframscript
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


(*FuncPtr={{ggtohhLHC14I, ggtohhLHC14II}, {ggtohHLHC14I, ggtohHLHC14II}};
Func2BeRun=FuncPtr[[proc]][[type]];*)
Print[type,proc];
