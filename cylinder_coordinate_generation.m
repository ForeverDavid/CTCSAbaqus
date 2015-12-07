(* ::Package:: *)

(* ::Input:: *)
(*p1.p1+p2.p2-2 p1.p2*)
(**)


(* ::Input:: *)
(*((p1.p1+p2.p2-2 p1.p2)/.{p1->p1i+dp1 t1,p2->p2i+dp2 t2})//tf//Expand*)


(* ::Input:: *)
(*tf[e_]:=e/.Dot[z1_+z2_,z3_+z4_]->Dot[z1,z3]+Dot[z2,z3]+Dot[z1,z4]+Dot[z2,z4]*)


(* ::Input:: *)
(*Map[#/.z_:>(z/.t1->1) t1^Count[z,t1,4]&,Map[#/.z_:>(z/.t2->1) t2^Count[z,t2,4]&,%]]*)


(* ::Input:: *)
(*FullSimplify[Solve[{D[%,t1]==0,D[%,t2]==0},{t1,t2}],TransformationFunctions->{Automatic,tf1}]*)


(* ::Input:: *)
(*tf1[e_]:=e/.Dot[z1_,z2_]->Dot[z2,z1]*)


(* ::Input:: *)
(*int[cyl1_,cyl2_]:=Module[{p1i=cyl1[[1,1]],dp1=cyl1[[1,2]]-cyl1[[1,1]],r1=cyl1[[2]],p2i=cyl2[[1,1]],dp2=cyl2[[1,2]]-cyl2[[1,1]],r2=cyl2[[2]],loc,t1,t2},loc={t1->(dp2.dp1 (-dp2.p1i+dp2.p2i)+dp2.dp2 (p1i.dp1-p2i.dp1))/((dp1.dp2)^2-dp1.dp1 dp2.dp2),t2->(dp1.dp1 (-dp2.p1i+dp2.p2i)+dp2.dp1 (p1i.dp1-p2i.dp1))/((dp1.dp2)^2-dp1.dp1 dp2.dp2)};*)
(*(-r2/Norm[dp2]<(t1/.loc)<1+r2/Norm[dp2])&&(-r1/Norm[dp1]<(t2/.loc)<1+r1/Norm[dp1])&&((t1^2 dp1.dp1-2 t1 t2 dp1.dp2+t1 dp1.p1i-2 t1 dp1.p2i+t2^2 dp2.dp2+t2 dp2.p2i+t1 p1i.dp1-2 t2 p1i.dp2+p1i.p1i-2 p1i.p2i+t2 p2i.dp2+p2i.p2i)/.loc)<(r1+r2)^2]*)


(* ::Input:: *)
(*cylinders=Table[{RandomReal[{-100,100},{2,3}],RandomReal[5]},{100}];*)
(*nint=ParallelTable[Or@@Table[int[cylinders[[i]],cylinders[[j]]],{j,i+1,100}],{i,100}];*)
(*c=Cases[{cylinders,nint}//Transpose,{z_,False}->z,{1}];*)
(*Graphics3D[{EdgeForm[None],Directive[Opacity@RandomReal[{.4,.9}],Hue[RandomReal[]]],Cylinder[First@#,Last@#]}&/@c,Boxed->False,ImageSize->800]*)
