(* ::Package:: *)

BeginPackage["PerfectGeometry`"]

SavePerfectGeometry::usage="Saves the arrays needed to define magnetic geometry in PERFECT.
SavePerfectGeometry[BHat,dBHatdpsi,dBHatdtheta,JHat,IHat,dIHatdpsi]"

Begin["Private`"]

Options[SavePerfectGeometry]={Filename->"input.geometry.h5",Append->False}

SavePerfectGeometry[BHat_,dBHatdpsi_,dBHatdtheta_,JHat_,IHat_,dIHatdpsi_,OptionsPattern[]]:=
Module[{Npsi,Ntheta,filenameval,appendval},
filenameval=OptionValue[Filename];
appendval=OptionValue[Append];
If[FileExtension[filenameval]!="h5",Throw["Filename suffix must be .h5 for Mathematica to write to and HDF5 file."]];
{Npsi,Ntheta}=Dimensions[BHat];
If[Dimensions[dBHatdpsi]!={Npsi,Ntheta},Throw["Dimensions of dBHatdpsi do not match dimensions of BHat"]];
If[Dimensions[dBHatdtheta]!={Npsi,Ntheta},Throw["Dimensions of dBHatdtheta do not match dimensions of BHat"]];
If[Dimensions[JHat]!={Npsi,Ntheta},Throw["Dimensions of JHat do not match dimensions of BHat"]];
If[Dimensions[IHat]!={Npsi},Throw["Dimensions of IHat do not match dimensions of BHat"]];
If[Dimensions[dIHatdpsi]!={Npsi},Throw["Dimensions of dIHatdpsi do not match dimensions of BHat"]];
Export[filenameval,
{
BHat,
dBHatdpsi,
dBHatdtheta,
JHat,
IHat,
dIHatdpsi
},
{"Datasets",
{
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/BHat"],
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/dBHatdpsi"],
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/dBHatdtheta"],
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/JHat"],
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/IHat"],
StringJoin["/Npsi",ToString[Npsi],"Ntheta",ToString[Ntheta],"/dIHatdpsiB"]
}
},"Append"->appendval];
]

End[]
EndPackage[]
