//
//  main.cpp
//  NS
//
//  Created by Diego Andrade on 4/28/20.
//  Copyright © 2020 Diego Andrade. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fcntl.h>
#include "datadef.h"
#include "matrix.hpp"
#include "chorin.hpp"
#include <iomanip>

using namespace std;


REAL ** RMATRIX (REAL v, int imax, int jmax)
{
    
    int i,j;
    
    REAL** m = new REAL*[imax];
    for(i = 0; i < imax; ++i){
        m[i] = new REAL[jmax];
        for(j=0;j<jmax; ++j){
            m[i][j]=v;
        }
    }

    return(m);
}


#pragma mark <F> MAIN FUNCTION
int main(int argc, char * argv[]) {
    
    // insert code here...
    std::cout << "Computational Fluid Dynamics Solver \n";
    std::cout << "Diego Andrade \n";
    std::cout << "2020 COVID edition \n";
    
    matrix obj_m;
    chorin obj_c;
    
    
    //Interior Points
    int imax = 125;
    int jmax = 125;
    
    //obj_c.chorinI(imax, jmax, m);
    
    int iter = 0;
    
    REAL xlength = 1;
    REAL ylenght = 1;
    REAL delx = xlength/((REAL)imax); // Spatial step interval in X direction
  
    //cout << "delx : " << delx << "value " << "\t";
    REAL dely = ylenght/((REAL)jmax); // Spatial step interval in Y direction
    REAL ustar = 1; // normalized velocity of lid
    
    // Time conditions
    // Initialize time value and time index
    REAL timet = 0;
    REAL timen = 0; // timet:time value,timen:time index
    REAL delt = 0.02;  //Time Step
    REAL tau = 0.5; //Safety factor for time step size control \[Tau]
    REAL tend = 40;
  
 
    
    //Pressure Itertion Data
    //int it = 0;  // SOR counter
    //REAL res = 0; // Norm of pressure equation residual
    REAL eps = 0.0001; // Stopping tolerance for pressure iteration
    REAL gamma = 0.9; //Upwind differencing factor \[Gamma]
    
   //Problem Dependent quatities
   //Rey= 1000;//Reynolds Number

    REAL Reynolds = 10000;
    REAL GX = 0; // Body forces
    REAL GY = 0; // Body forces

    //REAL UI = 0; // Initial Velocity
    //REAL VI = 0; // Initial Velocity

    //REAL PI = 0; // Initial Pressure

    //REAL Umed = 1; // Mean velocity Value of the lid

    int itermax = 100; // Max iteration in SOR

    REAL omega = 1.7; // SOR Relaxation
    //REAL TOL = 0.001; // SOR Tolerance
    REAL ritnorm = 1;

    //int rindex = 1;
    
    //Data Array Initialization *)
    //initialize and assign initial values to u,v,and p*)
   // REAL ** m = RMATRIX(0.,  imax+2,  jmax+2);
    REAL ** U = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** V = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** P = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** PNEW = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** F = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** G = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** RHS = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** RIT = RMATRIX(0.0,  imax+2,  jmax+2);
    //REAL ** solution = RMATRIX(0.,  imax+2,  jmax+2);
    
    REAL ** d2ux = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d2uy = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d1u2x = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d1uvy = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d2vx = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d2vy = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d1v2y = RMATRIX(0.0,  imax+2,  jmax+2);
    REAL ** d1uvx = RMATRIX(0.0,  imax+2,  jmax+2);
    
    
    
    
    int wW = 2;
    int wE = 2;
    int wN = 2;
    int wS = 2;  // flags to indicate no slip condition on a wall
    
    U = obj_c.boundaryValuesU(U, ustar, imax, jmax, wN, wE, wW, wS);
    V = obj_c.boundaryValuesV(V, ustar, imax, jmax, wN, wE, wW, wS);
    
    // compute (d^2 u/d x^2)
    d2ux  = obj_c.deriv2ux(U, d2ux, imax, jmax, delx);
    // compute (d^2 u/d y^2)
    d2uy  = obj_c.deriv2uy(U, d2uy, imax, jmax, delx);
    // compute (d u^2/d x)
    d1u2x = obj_c.deriv1u2x(U, d1u2x,imax, jmax, delx, gamma);
    // compute (d uv/d y)
    d1uvy = obj_c.deriv1uvy(U, V, d1uvy, imax, jmax, dely, gamma);
    // compute (d^2 v/d x^2)
    d2vx  = obj_c.deriv2vx(V, d2vx, imax, jmax, delx);
    // compute (d^2 v/d y^2)
    d2vy = obj_c.deriv2vy(V, d2vy, imax, jmax, dely);
    // compute (d v^2/d y)
    d1v2y = obj_c.deriv1v2y(V, d1v2y, imax, jmax, dely, gamma);
    // compute (d uv/d x)
    d1uvx = obj_c.deriv1uvx(U, V, d1uvx ,imax, jmax, delx, gamma);
    
    F = obj_c.compF(F, U, imax, jmax, delx, delt, Reynolds, GX, d2ux, d2uy, d1u2x, d1uvy);
   // obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"F");
    G = obj_c.compG(G, V, imax, jmax, dely, delt, Reynolds, GY, d2vx, d2vy, d1uvx, d1v2y);
   // obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"G");
    RHS = obj_c.computeRHS(RHS, imax, jmax, delt, delx, dely, F, G);
  //  obj_m.PRINT_MATRIX(RHS, imax+2, jmax+2,"RHS");
  //  obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"P");
  
    iter=0;
    ritnorm =10;
    while(iter < itermax && ritnorm > eps)
    {
       ////obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.ant");
       PNEW = obj_c.computepNew(imax,  jmax,  omega,  delx,  dely, P, PNEW, RHS);
        //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.pos");
        //obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"PNEW.pos");
        RIT = obj_c.computeRit(RIT,imax, jmax, delx, dely, P, RHS);
        ritnorm = obj_c.normFrobeius(imax, jmax, RIT)/((imax + 2)*(jmax + 2));
        P = PNEW;
        iter = iter + 1;
        cout << setw(10) << "iter# : " << iter << " \t ritnorm: " << ritnorm << "\n";
   }
  
    
    U = obj_c.computeU(U, imax, jmax, delt, delx, F, P);
    V = obj_c.computeV(V, imax, jmax, delt, dely, G, P);
    
  //  obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"U---");
  //  obj_m.PRINT_MATRIX(V, imax+2, jmax+2,"V---");
  //  obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P---");
    //
   
  //  obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"PNEW.");
//    obj_m.PRINT_MATRIX(RIT, imax+2, jmax+2,"RIT.");
   // obj_m.PRINT_MATRIX(d2uy, imax+2, jmax+2,"d2uy");
   // obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P");
   // obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"PNEW");
  //  obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"U.");
  //  obj_m.PRINT_MATRIX(V, imax+2, jmax+2,"V.");
  //  obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.");
    //  obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"F");
    //  obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"G");
    
   
    int k=1;
  
    while(timet<tend)
    {
        
        delt = obj_c.select_delt(imax, imax, tau, Reynolds,  delx, dely, U, V);
        //cout << " delt : " <<delt << endl;
        //cout << " ------------------------------------------------------------------ " << endl;
        
        U = obj_c.boundaryValuesU(U, ustar, imax, jmax, wN, wE, wW, wS);
        V = obj_c.boundaryValuesV(V, ustar, imax, jmax, wN, wE, wW, wS);
        
        //obj_m.PRINT_MATRIX(V, imax+2, jmax+2,"Vbound");
        
        
        // compute (d^2 u/d x^2)
        d2ux  = obj_c.deriv2ux(U, d2ux, imax, jmax, delx);
        // compute (d^2 u/d y^2)
        d2uy  = obj_c.deriv2uy(U, d2uy, imax, jmax, delx);
        // compute (d u^2/d x)
        d1u2x = obj_c.deriv1u2x(U, d1u2x,imax, jmax, delx, gamma);
        // compute (d uv/d y)
        d1uvy = obj_c.deriv1uvy(U, V, d1uvy, imax, jmax, dely, gamma);
        // compute (d^2 v/d x^2)
        d2vx  = obj_c.deriv2vx(V, d2vx, imax, jmax, delx);
        // compute (d^2 v/d y^2)
        d2vy = obj_c.deriv2vy(V, d2vy, imax, jmax, dely);
        // compute (d v^2/d y)
        d1v2y = obj_c.deriv1v2y(V, d1v2y, imax, jmax, dely, gamma);
        // compute (d uv/d x)
        d1uvx = obj_c.deriv1uvx(U, V, d1uvx ,imax, jmax, delx, gamma);
           
    
        /*obj_m.PRINT_MATRIX(d2ux, imax+2, jmax+2,"d2ux");
        obj_m.PRINT_MATRIX(d2uy, imax+2, jmax+2,"d2uy");
        obj_m.PRINT_MATRIX(d1u2x, imax+2, jmax+2,"d1u2x");
        obj_m.PRINT_MATRIX(d1uvy, imax+2, jmax+2,"d1uvy");
    
        obj_m.PRINT_MATRIX(d2vx, imax+2, jmax+2,"d2vx");
        obj_m.PRINT_MATRIX(d2vy, imax+2, jmax+2,"d2vy");
        obj_m.PRINT_MATRIX(d1v2y, imax+2, jmax+2,"d1v2y");
        obj_m.PRINT_MATRIX(d1uvx, imax+2, jmax+2,"d1uvx");*/
        
        //cout << "iter: " << k << " ------------------------------------------------------------------ " << endl;
        //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P");
        //obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"F");
        //obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"G");
   
        P = obj_c.boundaryValuesP(P, imax, jmax, wN, wS, wE, wW);
        F = obj_c.boundaryValuesF(U, F, imax, jmax, wN, wS, wE, wW);
        G = obj_c.boundaryValuesG(V, G, imax, jmax, wN, wS, wE, wW);
        
        //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"-P");
        //obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"-F");
       // obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"-U");
        //obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"-G");
    
        F = obj_c.compF(F, U, imax, jmax, delx, delt, Reynolds, GX, d2ux, d2uy, d1u2x, d1uvy);
        G = obj_c.compG(G, V, imax, jmax, dely, delt, Reynolds, GY, d2vx, d2vy, d1uvx, d1v2y);
    
        
        // obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.");
        //obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"F.");
        //
        //obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"G.");
        
        RHS = obj_c.computeRHS(RHS,imax, jmax, delt, delx, dely, F, G);
        //obj_m.PRINT_MATRIX(RHS, imax+2, jmax+2,"RHS.");
       
       
        iter = 0;
        ritnorm = 1;
        while(iter < itermax && ritnorm > eps)
           {
              //cout << "\niteration# : " << iter << " ====================================================================================================================================================================="<< "\n";
              //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.ant");
              PNEW = obj_c.computepNew(imax,  jmax,  omega,  delx,  dely, P, PNEW, RHS);
              //obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"PNEW.pos");
              //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.in");
              RIT = obj_c.computeRit(RIT,imax, jmax, delx, dely, P, RHS);
              //obj_m.PRINT_MATRIX(RIT, imax+2, jmax+2,"RIT.in");
              ritnorm = obj_c.normFrobeius(imax, jmax, RIT)/((imax + 2)*(jmax + 2));
              P = PNEW;
              //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P.in");
              iter = iter + 1;
              //cout << "iter# : " << iter << " ritnorm: " << ritnorm << " =====================================================================================================================================================================" << "\n";
          }
    
        //obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"-U");
        U = obj_c.computeU(U, imax, jmax, delt, delx, F, P);
        //obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"U.");
        V = obj_c.computeV(V, imax, jmax, delt, dely, G, P);
        
        timet=timet+delt;
        timen=timen+1;
       if(k%100==0)
       {
        cout << "iteration k: " << k << " time: " << timet << endl;
       }
        if(k%25000==0)
        {
            cout << " ------------------------------------------------------ ------------------------------------------------------" <<endl;
            cout << "iteration k: " << k << " time: " << timet << endl;
            //obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"U");
           // obj_m.PRINT_MATRIX(V, imax+2, jmax+2,"V");
        }
        k++;
    }


    obj_m.RMATRIX_TO_FILE(U,"/Users/diegoandrade/Documents/NS/U.txt", imax+2, jmax+2);
    obj_m.RMATRIX_TO_FILE(V,"/Users/diegoandrade/Documents/NS/V.txt", imax+2, jmax+2);

    //obj_m.PRINT_MATRIX(P, imax+2, jmax+2,"P");
    //obj_m.PRINT_MATRIX(PNEW, imax+2, jmax+2,"PNEW");
    obj_m.PRINT_MATRIX(U, imax+2, jmax+2,"U");
    obj_m.PRINT_MATRIX(V, imax+2, jmax+2,"V");
    //obj_m.PRINT_MATRIX(F, imax+2, jmax+2,"F");
    //obj_m.PRINT_MATRIX(G, imax+2, jmax+2,"G");
    //obj_m.PRINT_MATRIX(rit, imax+2, jmax+2,"rit");
  

    obj_m.FREE_RMATRIX(U, imax);
    obj_m.FREE_RMATRIX(V, imax);
    obj_m.FREE_RMATRIX(PNEW, imax);
    obj_m.FREE_RMATRIX(F, imax);
    obj_m.FREE_RMATRIX(G, imax);
    //obj_m.FREE_RMATRIX(rit, imax);*/
    
    return 0;
}
