(* ::Package:: *)

(* ::Input::Initialization:: *)
(* :Title: NavierStokes *)

(* :Authors: Diego Andrade *)

(* :Context: Navier-Stokes` *)

(* :Package Version: 12.0 *)

(* :Copyright: Copyright 2020 Intelestial Research Inc.  *)

(* :History:
	Version 0.1 by Diego Andrade 2020.
*)

(* :Keywords: *)

(* :Source: *)

(* :Warning: *)

(* :Wolfram Language Version: 12.0 *)

(* :Limitation: *)



(* ::Input::Initialization:: *)
BeginPackage["NavierStokes`"]
Unprotect@@Names["NavierStokes`*"];
ClearAll@@Names["NavierStokes`*"];

fe::usage="f[x]"
ge::usage="g[x]"
deriv2ux::usage="This function computes \!\(\*FractionBox[\(\*SuperscriptBox[\(\[PartialD]\), \(2\)]u\), \(\[PartialD]\*SuperscriptBox[\(x\), \(2\)]\)]\)"
deriv2uy::usage="This function computes  \!\(\*FractionBox[\(\*SuperscriptBox[\(\[PartialD]\), \(2\)]u\), \(\[PartialD]\*SuperscriptBox[\(y\), \(2\)]\)]\) "
deriv1u2x::usage="This function computes \!\(\*FormBox[\(TraditionalForm\`\*FractionBox[\(\[PartialD]\*SuperscriptBox[\(u\), \(2\)]\), \(\[PartialD]x\)]\),
TraditionalForm]\)"

deriv2vx::usage="This function computes \!\(\*FormBox[FractionBox[\(\*SuperscriptBox[\(\[PartialD]\), \(2\)]v\), \(\[PartialD]\*SuperscriptBox[\(x\), \(2\)]\)],
TraditionalForm]\)"
deriv1uvy::usage="This function computes \!\(\*FractionBox[\(\[PartialD]\((uv)\)\), \(\[PartialD]y\)]\)"
deriv2vy::usage="This function computes \!\(\*FormBox[FractionBox[\(\*SuperscriptBox[\(\[PartialD]\), \(2\)]v\), \(\[PartialD]\*SuperscriptBox[\(y\), \(2\)]\)],
TraditionalForm]\)"
deriv1v2y::usage="This function computes \!\(\*FormBox[FractionBox[\(\[PartialD]\*SuperscriptBox[\(v\), \(2\)]\), \(\[PartialD]y\)],
TraditionalForm]\)"
deriv1uvx::usage="This function computes \!\(\*FormBox[\(TraditionalForm\`\*FractionBox[\(\[PartialD]\((uv)\)\), \(\[PartialD]x\)]\),
TraditionalForm]\)"
compF::usage="This function computes the F projection"

compG::usage="This function computes the G projection"

computeRHS::usage="This function computes RHS"

computepNew::usage="This function computes p of the next step, p_new"

epsE::usage="Epsilon w.r.t. right wall,  ret=0, if i=imax+1, ret=1, if i<imax+1"

epsN::usage="Epsilon w.r.t. top wall,  ret=0, if j=jmax+1, ret=1, if j<jmax+1"

epsS::usage="Epsilon w.r.t. bottom wall,  ret=0, if j=2, ret=1, if j>2"

epsW::usage="Epsilon w.r.t. left wall, ret=0, if i=2, ret=1, if i>2"

computeRit::usage="Compute residual for pressure Equation"

computeU::usage="Compute U for the next step"

computeV::usage="Compute V for the next step"

selectDelta::usage="Select delta t for the next time step"

boundaryValuesU::usage="set boundary values for velocity U"

boundaryValuesV::usage="set boundary values for velocity V"

boundaryValuesP::usage="set boundary values for pressure P"

boundaryValuesF::usage="set boundary values for pressure F"

boundaryValuesG::usage="set boundary values for pressure G"

Begin["`Private`"]

fe[x_]:=Module[{},x^2];

ge[x_]:=Module[{},x^8];

(* Second derivative of u with respect to (w.r.t.) x \[PartialD]^2u/\[PartialD]x^2*)
deriv2ux[MatU_,imax_,jmax_,dx_]:=
Module[{U=MatU,im=imax,jm=jmax,delx=dx,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\)-2\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\))/delx^2;
	];];
Return[solution];
];

(* First derivative of u^2 w.r.t.x \[PartialD]^2u/\[PartialD]y^2*)
deriv2uy[MatU_,imax_,jmax_,dy_]:=
Module[{U=MatU,im=imax,jm=jmax,dely=dy,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\)-2\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\))/dely^2;
		(*Print["S[",i,",",j,"]=",solution[[i,j]]];*)
	];];
Return[solution];
];

(* First derivative of u^2 w.r.t.x \[PartialD]u^2/\[PartialD]x*)
deriv1u2x[MatU_,imax_,jmax_,dx_,gamma_]:=
Module[{U=MatU,im=imax,jm=jmax,delx=dx,gm=gamma,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=1/delx (((\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\))/2)^2-((\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2)^2)+
gm/delx (Abs[\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\)]/2 (\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\))/2-Abs[\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)]/2 (\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2);
	];];
Return[solution];
];


(* First derivative of uv w.r.t.y \[PartialD](uv)/\[PartialD]y *)
deriv1uvy[MatU_,MatV_,imax_,jmax_,dy_,gamma_]:=
Module[{U=MatU,V=MatV,im=imax,jm=jmax,dely=dy,gm=gamma,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=1/dely*((\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\))/2 *(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2-(\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j - 1\)\(\[RightDoubleBracket]\)\)]\))/2 *(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2)+
gamma/dely (Abs[\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\)]/2 *(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2-Abs[\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j - 1\)\(\[RightDoubleBracket]\)\)]\)]/2 *(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2);
	];];
Return[solution];
];

(* Second derivative of v w.r.t.x \[PartialD]^2v/\[PartialD]x^2 *)
deriv2vx[MatV_,imax_,jmax_,dx_]:=
Module[{V=MatV,im=imax,jm=jmax,delx=dx,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
For[i=2,i<im+2,i++,
	For[j=2,j<jm+1,j++,
		solution[[i,j]]=(\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\)-2\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\))/delx^2;
	];];
Return[solution];
];

(* Second derivative of v w.r.t.y \[PartialD]^2v/\[PartialD]y^2 *)
deriv2vy[MatV_,imax_,jmax_,dy_]:=
Module[{V=MatV,im=imax,jm=jmax,dely=dy,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
For[i=2,i<im+2,i++,
	For[j=2,j<jm+1,j++,
		solution[[i,j]]=(\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\)-2\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\))/dely^2;
	];];
Return[solution];
];

(* First derivative of v^2 w.r.t.y \[PartialD]v^2/\[PartialD]y *)
deriv1v2y[MatV_,imax_,jmax_,dy_,gamma_]:=
Module[{V=MatV,im=imax,jm=jmax,dely=dy,gm=gamma,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=1/dely*(((\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2)^2-((\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2)^2)+
gamma/dely (Abs[\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\)]/2*(\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2-Abs[\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)]/2*(\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j - 1\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2);
	];];
Return[solution];
];

(* First derivative of uv w.r.t.x \[PartialD](uv)/\[PartialD]x *)
deriv1uvx[MatU_,MatV_,imax_,jmax_,dx_,gamma_]:=
Module[{U=MatU,V=MatV,im=imax,jm=jmax,delx=dx,gm=gamma,solution,i,j},
solution= Table[0,{im+2},{jm+2}];
	For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		solution[[i,j]]=1/delx ((\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2 (\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\))/2-(\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j + 1\)\(\[RightDoubleBracket]\)\)]\))/2 (\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2)+
gamma/delx (Abs[\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i, j + 1\)\(\[RightDoubleBracket]\)\)]\)]/2 (\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i + 1, j\)\(\[RightDoubleBracket]\)\)]\))/2-Abs[\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)+\!\(\*SubscriptBox[\(U\), \(\(\[LeftDoubleBracket]\)\(i - 1, j + 1\)\(\[RightDoubleBracket]\)\)]\)]/2 (\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i - 1, j\)\(\[RightDoubleBracket]\)\)]\)-\!\(\*SubscriptBox[\(V\), \(\(\[LeftDoubleBracket]\)\(i, j\)\(\[RightDoubleBracket]\)\)]\))/2);
	];];
Return[solution];
];

(* Computes F projection according to Equation 3.29 *)
compF[MatU_,imax_,jmax_,dx_,dt_,Re_,gx_,d2ux_,d2uy_,d1u2x_,d1uvy_]:=
Module[{U=MatU,im=imax,jm=jmax,delx=dx,delt=dt,i,j},
F= Table[0,{imax+2},{jmax+2}];(* F *)
For[i=2,i<im+1,i++,
	For[j=2,j<jm+2,j++,
		F[[i,j]] =U[[i,j]]+delt*(1/Re*(d2ux[[i,j]]+d2uy[[i,j]])-d1u2x[[i,j]]-d1uvy[[i,j]]+gx);

	];];
Return[F];
];


(* Computes G projection according to Equation 3.29 *)
compG[MatV_,imax_,jmax_,dy_,dt_,Re_,gy_,d2vx_,d2vy_,d1uvx_,d1v2y_]:=
Module[{V=MatV,im=imax,jm=jmax,dely=dy,delt=dt,i,j},
G= Table[0,{imax+2},{jmax+2}];(* G *)
For[i=2,i<im+2,i++,
	For[j=2,j<jm+1,j++,
		G[[i,j]] =V[[i,j]]+delt*(1/Re*(d2vx[[i,j]]+d2vy[[i,j]])-d1uvx[[i,j]]-d1v2y[[i,j]]+gy);

	];];
Return[G];
];

(* Compute right-hand size of the pressure Eq.(3.38) *)
computeRHS[imax_,jmax_,delt_,delx_,dely_,F_,G_]:=
Module[{i,j,rhs},
rhs= Table[0,{imax+2},{jmax+2}];(* RHS *)
For[i=2,i<imax+2,i++,
	For[j=2,j<jmax+2,j++,
		rhs[[i,j]] = 1/delt*((F[[i,j]]-F[[i-1,j]])/delx+(G[[i,j]]-G[[i,j-1]])/dely);

	];];
Return[rhs];
];

epsE[i_,imax_]:=
Module[{ret},
ret=0.0;
If[i==imax+1,ret=0.0,ret=1.0];
Return[ret];
];

epsN[j_,jmax_]:=
Module[{ret},
ret=0.0;
If[j==jmax+1,ret=0.0,ret=1.0];
Return[ret];
];

epsS[j_]:=
Module[{ret},
ret=0.0;
If[j==2,ret=0.0,ret=1.0];
Return[ret];
];

epsW[i_]:=
Module[{ret},
ret=0.0;
If[i==2,ret=0.0,ret=1.0];
Return[ret];
];

(* compute p of the next step,p_new Eq. 3.44*) 
computepNew[imax_,jmax_,omega_,delx_,dely_,p_,rhs_]:=
Module[{pnew,i,j},
pnew= Table[0,{imax+2},{jmax+2}];(* pressure at time it+1 *)
For[j=2,j<jmax+2,j++,
	pnew[[1,j]]=p[[2,j]];
pnew[[imax+2,j]]=p[[imax+1,j]];
	];
For[i=2,i<imax+2,i++,
	pnew[[i,1]]=p[[i,2]];
pnew[[i,jmax+2]]=p[[i,jmax+1]];
	];
For[i=2,i<imax+2,i++,
	For[j=2,j<jmax+2,j++,
		pnew[[i,j]]=(1-omega)*p[[i,j]]+omega/((epsE[i,imax]+epsW[i])/delx^2+(epsN[j,jmax]+epsS[j])/dely^2)*((epsE[i,imax]*p[[i+1,j]]+epsW[i]*pnew[[i-1,j]])/delx^2+(epsN[j,jmax]*p[[i,j+1]]+epsS[j]*pnew[[i,j-1]])/dely^2-rhs[[i,j]]);

];];

Return[pnew];
];

(* compute residual for pressure Eq.according to Eq.(3.45) *)
computeRit[imax_,jmax_,delx_,dely_,p_,rhs_]:=
Module[{i,j},
rit= Table[0,{imax+2},{jmax+2}];
For[i=2,i<imax+2,i++,
	For[j=2,j<jmax+2,j++,
		rit[[i,j]]=(epsE[i,imax]*(p[[i+1,j]]-p[[i,j]])-epsW[i]*(p[[i,j]]-p[[i-1,j]]))/delx^2+(epsN[j,jmax]*(p[[i,j+1]]-p[[i,j]])- epsS[j]*(p[[i,j]]-p[[i,j-1]]))/dely^2-rhs[[i,j]];

];];
Return[rit];
];

(* compute u and v of the next step according to Eqs.(3.34) and (3.35) *)
computeU[MatU_,imax_,jmax_,delt_,delx_,F_,p_]:=
Module[{U=MatU,i,j},
For[i=2,i<imax+1,i++,
	For[j=2,j<jmax+2,j++,
		U[[i,j]]=F[[i,j]]-delt/delx*(p[[i+1,j]]-p[[i,j]]);
];];
Return[U];
];

computeV[MatV_,imax_,jmax_,delt_,dely_,G_,p_]:=
Module[{V=MatV,i,j},
For[i=2,i<imax+2,i++,
	For[j=2,j<jmax+1,j++,
		V[[i,j]]=G[[i,j]]-delt/dely*(p[[i,j+1]]-p[[i,j]]);
];];
Return[V];
];

(* Select delta t for the next time step using Eq.(3.50) *)
selectDelta[tau_,Re_,delx_,dely_,MatU_,MatV_]:=
Module[{U=MatU,V=MatV,delta,i,j},

		delta=tau*Min[Re/2*(1/delx^2+1/dely^2)^-1,delx/If[Abs[Max[U]]!=0,Abs[Max[U]],1],dely/If[Abs[Max[V]]!=0,Abs[Max[V]],1]];

Return[delta];
];

(* set boundary values for velocity u and v *)
boundaryValuesU[MatU_,ustar_,imax_,jmax_,wN_,wE_,wW_,wS_]:=
Module[{U=MatU,i,j},
(* w.r.t.u *)
(*left wall *)
	If[wW==2,
For[j=2,j<jmax+2,j++,
		U[[1,j]]=0.0;
];
];
(*right wall *)
	If[wE==2,
For[j=2,j<jmax+2,j++,
		U[[imax+2,j]]=0.0;
];
];
(*bottom wall *)
	If[wS==2,
For[i=2,i<imax+2,i++,
		U[[i,1]]=-U[[i,2]];
];
];
(*top wall *)
	If[wN==2,
For[i=2,i<imax+2,i++,
		U[[i,jmax+2]]=2*ustar-U[[i,jmax+1]];
];
];
	
Return[U];
];

boundaryValuesV[MatV_,imax_,jmax_,wN_,wE_,wW_,wS_]:=
Module[{V=MatV,i,j},
(* w.r.t.u *)
(*top wall *)
	If[wN==2,
For[i=2,i<imax+2,i++,
		V[[i,jmax+2]]=0.0;
];
];
(*bottom wall *)
	If[wS==2,
For[i=2,i<imax+2,i++,
		V[[i,1]]=0.0;
];
];
(*left wall *)
If[wW==2,
For[j=2,j<jmax+2,j++,
		V[[1,j]]=-V[[2,j]];
];
];
(*right wall *)
	If[wE==2,
For[j=2,j<jmax+2,j++,
		V[[imax+2,j]]=-V[[imax+1,j]];
];
];

Return[V];
];

boundaryValuesP[MatP_,imax_,jmax_,wN_,wE_,wW_,wS_]:=
Module[{P=MatP,i,j},
(* w.r.t.u *)
(* left wall *)
	If[wW==2,
For[j=2,j<jmax+2,j++,
		P[[1,j]]=P[[2,j]];
];
];
(* right wall *)
	If[wE==2,
For[j=2,j<jmax+2,j++,
		P[[imax+2,j]]=P[[imax+1,j]];
];
];
(* bottom wall *)
If[wS==2,
For[i=2,i<imax+2,i++,
		P[[i,1]]=P[[i,2]];
];
];
(* top wall *)
	If[wN==2,
For[i=2,i<imax+2,i++,
		P[[i,jmax+2]]=P[[i,jmax+1]];
];
];

Return[P];
];

boundaryValuesF[MatU_,imax_,jmax_,wN_,wE_,wW_,wS_]:=
Module[{U=MatU,F,i,j},
F= Table[0,{imax+2},{jmax+2}];
(* w.r.t.u *)
(* bottom wall *)
If[wW==2,
For[j=2,j<jmax+2,j++,
		F[[1,j]]=U[[1,j]];
];
];
(* right wall *)
	If[wE==2,
For[j=2,j<jmax+2,j++,
		F[[imax+1,j]]=U[[imax+1,j]];
];
];

Return[F];
];

boundaryValuesG[MatV_,imax_,jmax_,wN_,wE_,wW_,wS_]:=
Module[{V=MatV,G,i,j},
G= Table[0,{imax+2},{jmax+2}];
(* bottom wall *)
If[wS==2,
For[i=2,i<imax+2,i++,
		G[[i,1]]=V[[i,1]];
];
];
	(* top wall *)
	If[wN==2,
For[i=2,i<imax+2,i++,
		G[[i,jmax+1]]=V[[i,jmax+1]];
];
];

Return[G];
];

End[]
Protect@@Names["NavierStokes`*"];
EndPackage[]






