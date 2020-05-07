//
//  chorin.hpp
//  NS
//
//  Created by Diego Andrade on 4/29/20.
//  Copyright © 2020 Diego Andrade. All rights reserved.
//

#ifndef chorin_hpp
#define chorin_hpp

#include <stdio.h>
#include <iostream>
#include <math.h>

#include "datadef.h"
#include "matrix.hpp"

using namespace std;

class chorin{
    
public:
    
    chorin();
    void chorinI(int imax, int jmax, REAL** m);
    ~chorin();
    
    REAL** solution;
    //REAL** F;
    //REAL** G;
    REAL** RHS;
    REAL** RIT;
    REAL** PNEW;
    REAL** P;
    REAL** U;
    REAL** V;
    
    /*REAL ** d2ux;
      // compute (d^2 u/d y^2)
    REAL ** d2uy;
      // compute (d u^2/d x)
    REAL ** d1u2x;
      // compute (d uv/d y)
    REAL ** d1uvy;
      // compute (d^2 v/d x^2)
    REAL ** d2vx ;
      // compute (d^2 v/d y^2)
    REAL ** d2vy;
      // compute (d v^2/d y)
    REAL ** d1v2y;
      // compute (d uv/d x)
    REAL ** d1uvx;*/
    
    
    REAL** m;

    
    REAL** boundaryValuesU(REAL** U, REAL ustar,
                           int imax, int jmax,
                           int wN, int wE, int wW, int wS);
    
    REAL** boundaryValuesV(REAL** V, REAL ustar,
                            int imax, int jmax,
                            int wN, int wE, int wW, int wS);
    
    REAL** deriv2ux(REAL** U, REAL** solution , int imax, int jmax, REAL delx);
    
    REAL** deriv2uy(REAL** U, REAL** d2uy, int imax, int jmax, REAL dely);
    
    REAL** deriv1u2x(REAL** U, REAL** solution, int imax, int jmax, REAL delx, REAL gamma);
    
    REAL** deriv1uvy(REAL** U, REAL** V, REAL** solution, int imax, int jmax, REAL dely, REAL gamma);
    
    REAL** deriv2vx(REAL** V, REAL** solution, int imax, int jmax, REAL delx);
    
    REAL** deriv2vy(REAL** V, REAL** solution, int imax, int jmax, REAL dely);
    
    REAL** deriv1v2y(REAL** V, REAL** solution, int imax, int jmax, REAL dely, REAL gamma);
    
    REAL** deriv1uvx(REAL ** U, REAL ** V, REAL** solution, int imax, int jmax, REAL delx, REAL gamma);
   
    REAL** compF(REAL** F, REAL** U,int imax,int  jmax,REAL delx,REAL delt, REAL Reynolds, REAL GX,
                 REAL** d2ux, REAL** d2uy, REAL** d1u2x, REAL** d1uvy);
    
    REAL** compG(REAL** G, REAL** V,int imax,int  jmax,REAL dely,REAL delt, REAL Reynolds, REAL GY,
                 REAL** d2vx, REAL** d2vy, REAL** d1uvx, REAL** d1v2y);
    
    REAL** computeRHS(REAL** RHS, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL ** F, REAL ** G);
    
    REAL epsE(int i, int imax);
    
    REAL epsW(int i);
    
    REAL epsN(int j, int jmax);
    
    REAL epsS(int j);
    
    REAL** computepNew(int imax, int jmax, REAL omega, REAL delx, REAL dely, REAL** P, REAL ** RHS);
    
    REAL** computeRit(REAL** RIT, int imax, int jmax, REAL delx, REAL dely, REAL** P, REAL** RHS);
    
    REAL normFrobeius(int ima, int jmax, REAL** M);
    
    REAL** computeU(REAL** U, int imax, int jmax, REAL delt, REAL delx, REAL** F, REAL** P);
    
    REAL** computeV(REAL** V, int imax, int jmax, REAL delt, REAL dely, REAL** G, REAL** P);
    
    REAL select_delt(int imax, int jmax, REAL tau,REAL Re,REAL delx,REAL dely,REAL** U,REAL** v);
    
    REAL maxMat(REAL** M, int imax, int jmax);
    
    REAL** boundaryValuesP (REAL** P, int imax, int jmax, int wN, int wS, int wE, int wW);
    
    REAL** boundaryValuesF (REAL** U, REAL**F, int imax, int jmax, int wN, int wS, int wE, int wW);
    
    REAL** boundaryValuesG (REAL** V, REAL**G, int imax, int jmax, int wN, int wS, int wE, int wW);
    
    REAL** RMATRIX_ZERO(REAL zero, int imax, int jmax);

    
};
#endif /* chorin_hpp */
