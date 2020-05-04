//
//  chorin.cpp
//  NS
//
//  Created by Diego Andrade on 4/29/20.
//  Copyright © 2020 Diego Andrade. All rights reserved.
//

#include "chorin.hpp"

REAL smallest(REAL x, REAL y, REAL z){
    return std::min({x, y, z});
}

matrix obj_m;
matrix * q = new matrix();

chorin::chorin()
{
    q->m;
}

chorin::~chorin()
{
    
}

void chorin::chorinI(int imax, int jmax, REAL** mo)
{
    
    //solution = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
    //F = obj_m.RMATRIX(m , 0.0, imax+2, jmax+2); //delete object??
    //G = obj_m.RMATRIX(m , 0.0, imax+2, jmax+2); //delete object??
    //RHS = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
   // RIT = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
   // PNEW = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
    //P = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
   // U = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
   // V = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
    m = obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
    
    /*d2ux=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2); //delete object??
    // compute (d^2 u/d y^2)
    d2uy=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // compute (d u^2/d x)
    d1u2x=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // compute (d uv/d y)
    d1uvy=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // compute (d^2 v/d x^2)
    d2vx=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // ^2 v/d y^2)
    d2vy=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // compute (d v^2/d y)
    d1v2y=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);
    // compute (d uv/d x)
    d1uvx=obj_m.RMATRIX(m, 0.0, imax+2, jmax+2);*/
      
}

REAL** chorin::boundaryValuesU(REAL** U, REAL ustar,
                               int imax, int jmax,
                               int wN, int wE, int wW, int wS)
{
    int i,j;
    
    //left wall
    if(wW==2)
    {
        for(j=1;j<jmax+1;j++)
            U[0][j]=0.0;
    }
    //right wall
    if(wE==2)
    {
        for(j=1;j<jmax+1;j++)
            U[imax+1][j]=0.0;
    }
    //bottom wall
    if(wS==2)
       {
        for(i=1;i<imax+1;i++)
            U[i][0]= -U[i][1];
       }
    //top wall
    if(wS==2)
       {
        for(i=1;i<imax+1;i++)
            U[i][jmax+1]= 2*ustar-U[i][jmax];
        }
    return (U);
}


REAL** chorin::boundaryValuesV(REAL** V, REAL vstar,
                               int imax, int jmax,
                               int wN, int wE, int wW, int wS)
{
    int i,j;
       //top wall
        if(wS==2)
        {
            for(i=1;i<imax+1;i++)
                V[i][jmax+1]= 0.0;
        }
       //bottom wall
       if(wS==2)
       {
           for(i=1;i<imax+1;i++)
               V[i][0]= 0.0;
       }
       //left wall
       if(wW==2)
       {
           for(j=1;j<jmax+1;j++)
               V[0][j]=-V[1][j];
       }
       //right wall
       if(wE==2)
       {
           for(j=1;j<jmax+1;j++)
               V[imax+1][j]=-V[imax][j];
       }

     
    return (V);
}

REAL** chorin::deriv2ux(REAL** U, REAL** d2ux , int imax, int jmax, REAL delx)
{
    int i,j;
    //p =RMATRIX_ZERO(solution, 0.0, imax+2, jmax+2); //delete object??
    
    
    //obj_m.PRINT_MATRIX(solution, imax+2, jmax+2,"solution-");
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            d2ux[i][j]=(U[i+1][j]-2*U[i][j]+U[i-1][j])/pow(delx,2);
            //cout << delx << "\t";
        }
    }
    //obj_m.PRINT_MATRIX(solution, imax+2, jmax+2,"solution+");
    return(d2ux);
}

REAL** chorin::deriv2uy(REAL** U, REAL** d2uy, int imax, int jmax, REAL dely)
{
    int i,j;
    // d2uy = RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            d2uy[i][j]=(U[i][j+1]-2*U[i][j]+U[i][j-1])/pow(dely,2);
            //cout << "d2uy[" << i <<"]["<<j<<"]: " <<solution[i][j] << endl;
        }
    }
    
    return(d2uy);
    
}

REAL** chorin::deriv1u2x(REAL** U, REAL** d1u2x, int imax, int jmax, REAL delx, REAL gamma)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            d1u2x[i][j]=1/delx*(pow((U[i][j]+U[i+1][j])/2,2)-pow((U[i-1][j]+U[i][j])/2,2))+
            (gamma/delx)*((fabs(U[i][j]+U[i+1][j])/2)*(U[i][j]-U[i+1][j])/2-(fabs(U[i-1][j]+U[i][j])/2)*(U[i-1][j]-U[i][j])/2);
        }
    }
    
    return(d1u2x);
}

REAL** chorin::deriv1uvy(REAL** U,  REAL** V, REAL** solution, int imax, int jmax, REAL dely, REAL gamma)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            solution[i][j]=
            1/(4*dely)*
            ((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])
             -(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])) +
            (gamma/(4*dely))*
            (fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])
             -fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]));
        }
    }
    
    return(solution);
}


REAL** chorin::deriv2vx(REAL** V, REAL** d2vx, int imax, int jmax, REAL delx)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
       
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax;j++){
               d2vx[i][j]= (V[i+1][j]-2*V[i][j]+V[i-1][j])/pow(delx,2);
           }
       }
       
    return(d2vx);
}

REAL** chorin::deriv2vy(REAL** V, REAL** d2vy, int imax, int jmax, REAL dely)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
       for(i=1;i<imax+1;i++){
           for(j=1;j<jmax;j++){
               d2vy[i][j]= (V[i][j+1]-2*V[i][j]+V[i][j+1])/pow(dely,2);
           }
       }
       
       return(d2vy);
}


REAL** chorin::deriv1v2y(REAL** V, REAL** solution, int imax, int jmax, REAL dely, REAL gamma)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            solution[i][j]= 1/dely*(pow((V[i][j]+V[i][j+1])/2,2)-pow((V[i][j-1]+V[i][j])/2,2)) +
            (gamma/dely)*(fabs(V[i][j]+V[i][j+1])/2*(V[i][j]-V[i][j+1])/2-fabs(V[i][j-1]+V[i][j])/2*(V[i][j-1]-V[i][j])/2);
        }
    }
    
    return(solution);
    
    
}


REAL** chorin::deriv1uvx(REAL ** U, REAL ** V, REAL** d1uvx, int imax, int jmax, REAL delx, REAL gamma)
{
    int i,j;
    //solution =obj_m.RMATRIX_ZERO(m, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            d1uvx[i][j]= 1/delx*((U[i][j]+U[i][j+1])/2*(V[i][j]+V[i+1][j])/2-(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]+V[i][j])/2) +
            (gamma/delx)*(fabs(U[i][j]+U[i][j+1])/2*(V[i][j]-V[i+1][j])/2-fabs(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]-V[i][j])/2);
        }
    }
    
    return(d1uvx);
}

REAL** chorin::compF(REAL** F, REAL** U,int imax,int  jmax,REAL delx,REAL delt, REAL Reynolds, REAL GX,
                     REAL** d2ux, REAL** d2uy, REAL** d1u2x, REAL** d1uvy)
{
    
    int i,j;

   // F =RMATRIX_ZERO(F, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax+2;j++){
            F[i][j] = U[i][j]+delt*(((1/Reynolds)*(d2ux[i][j]+d2uy[i][j]))-d1u2x[i][j]-d1uvy[i][j]+GX);
           
        }
    }
    
    return(F);
}

REAL** chorin::compG(REAL** G, REAL** V,int imax,int  jmax,REAL dely,REAL delt, REAL Reynolds, REAL GY,
                     REAL** d2vx, REAL** d2vy, REAL** d1uvx, REAL** d1v2y)
{
    
    int i,j;
    
   // G = RMATRIX_ZERO(G, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax;j++){
            G[i][j] = V[i][j]+delt*((1/Reynolds)*(d2vx[i][j]+d2vy[i][j])-d1uvx[i][j]-d1v2y[i][j]+GY);
        }
    }
    
    return(G);
}


 REAL** chorin::computeRHS(REAL** RHS, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL ** F, REAL ** G)
{
    int i,j;
    //RHS =RMATRIX_ZERO(RHS, 0.0, imax+2, jmax+2); //delete object??
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax+1;j++){
            RHS[i][j] = 1/delt*(F[i][j]-F[i-1][j])/delx+(G[i][j]-G[i][j-1])/dely;
        }
    }
    
    return(RHS);
    
}

REAL chorin::epsE(int i, int imax)
{
    REAL ret = 0.0;
    if(i==imax)
    {
        ret=0.0;
    }else{
        ret=1.0;
    }
    
    return (ret);
}

REAL chorin::epsN(int j, int jmax)
{
    REAL ret = 0.0;
    if(j==jmax)
    {
        ret=0.0;
    }else{
        ret=1.0;
    }
    
    return (ret);
}

REAL chorin::epsS(int j)
{
    REAL ret = 0.0;
    if(j==1)
    {
        ret=0.0;
    }else{
        ret=1.0;
    }
    
    return (ret);
}

REAL chorin::epsW(int i)
{
    REAL ret = 0.0;
    if(i==1)
    {
        ret=0.0;
    }else{
        ret=1.0;
    }
    
    return (ret);
}

REAL** chorin::computepNew(REAL** PNEW, int imax, int jmax, REAL omega, REAL delx, REAL dely, REAL** P, REAL ** RHS)
{
    int i,j;
    REAL rdx2,rdy2;

    //PNEW =RMATRIX_ZERO(PNEW, 0.0, imax+2, jmax+2); //delete object??
    
    for(j=1;j<jmax+1;j++){
        PNEW[0][j]=P[1][j];
        PNEW[imax+1][j]=P[imax][j];
    }
    
    for(i=1;i<imax+1;i++){
         PNEW[i][0]=P[i][1];
         PNEW[i][jmax+1]=P[i][jmax];
     }
    
    rdx2 = 1./pow(delx,2);
    rdy2 = 1./pow(dely,2);
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax+1;j++){
            PNEW[i][j]=(1-omega)*P[i][j]+
                        (omega/(((epsE(i,imax) + epsW(i))*rdx2)+((epsN(j,jmax) + epsS(j))*rdy2))) *
                        ((epsE(i,imax)*P[i+1][j]+epsW(i)*PNEW[i-1][j])*rdx2 +
                        (epsN(j,jmax)*P[i][j+1]+epsS(j)*PNEW[i][j-1])*rdy2 -
                         RHS[i][j]);
        }
    }
    
    return(PNEW);
    
}


REAL** chorin::computeRit(REAL** RIT, int imax, int jmax, REAL delx, REAL dely, REAL** P, REAL** RHS)
{
    int i,j;
    REAL rdx2,rdy2;
    //matrix o_matrix;
    //RIT =RMATRIX_ZERO(RIT, 0.0, imax+2, jmax+2); //delete object??
    
    rdx2 = 1./delx/delx;
    rdy2 = 1./dely/dely;
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax+1;j++){
            RIT[i][j]= (epsE(i,imax)*(P[i+1][j]-P[i][j])-epsW(i)*(P[i][j]-P[i-1][j]))*rdx2 +
                        (epsN(j,jmax)*(P[i][j+1]-P[i][j])-epsS(j)*(P[i][j]-P[i][j-1]))*rdy2 - RHS[i][j];
                        
        }
    }
    
    return(RIT);
    
    
}

REAL chorin::normFrobeius(int imax, int jmax, REAL** M)
{
    int i,j;
    REAL norm=0;
    
    for(i=0;i<imax+2;i++){
        for(j=0;j<jmax+2;j++){
            norm += pow(fabs(M[i][j]),2);
        }
    }
    
    norm = sqrt(norm);
    
    return(norm);
}

REAL** chorin::computeU(REAL** U, int imax, int jmax, REAL delt, REAL delx, REAL** F, REAL** P)
{
    int i,j;
    
    for(i=1;i<imax;i++){
        for(j=1;j<jmax+1;j++){
            U[i][j]=F[i][j]-delt/delx*(P[i+1][j]-P[i][j]);
            /*
            cout << "U[" << i <<"]["<<j<<"]: "<<U[i][j]<<endl;
            cout << "F[" << i <<"]["<<j<<"]: "<<F[i][j]<<endl;
            cout << "P[" << i+1<<"]["<<j<<"]: "<<P[i+1][j]<<endl;
            cout << "P[" << i <<"]["<<j<<"]: "<<P[i][j]<<endl;
            cout << "delx :  " <<delx<<endl;
            cout << "delt :  " <<delx<<endl;
             */
        }
    }
    
    return(U);
}


REAL** chorin::computeV(REAL** V, int imax, int jmax, REAL delt, REAL dely, REAL** G, REAL** P)
{
    int i,j;
    
    for(i=1;i<imax+1;i++){
        for(j=1;j<jmax;j++){
            V[i][j]=G[i][j]-delt/dely*(P[i][j+1]-P[i][j]);
        }
    }
    
    return(V);
    
}


REAL chorin::select_delt(int imax, int jmax, REAL tau, REAL Re, REAL delx, REAL dely, REAL **U, REAL **V)
{
    REAL delta = 0;
    REAL dx2,dy2,maxU,maxV;
    dx2=1./delx/delx;
    dy2=1./dely/dely;
    
    maxU = maxMat(U, imax, jmax);
    maxV = maxMat(V, imax, jmax);
    
    delta = min({Re/2*pow((dx2+dy2),-1), pow(delx,2)/abs(maxU),pow(dely,2)/abs(maxV)});
    
    return(delta);
}

                
REAL chorin::maxMat(REAL** M, int imax, int jmax)
{
    REAL maxVal = 0.;
    int i,j;
    
    for(i=1;i<imax+2;i++){
        for(j=1;j<jmax+2;j++){
            maxVal = max(M[i][j],maxVal);
        }
    }
    
    return(maxVal);
}


REAL** chorin::boundaryValuesP (REAL** P, int imax, int jmax, int wN, int wS, int wE, int wW)
{
    int i,j;
    //left wall
    if(wW==2)
    {
        for(j=1;j<jmax;j++)
            P[0][j] = P[1][j];
    }
    //right wall
    if(wE==2)
    {
        for(j=1;j<jmax;j++)
            P[imax+1][j] = P[imax][j];
    }
    //bottom wall
    if(wS==2)
    {
        for(i=1;i<imax;i++)
            P[i][0] = P[i][1];
    }
    //top wall
    if(wN==2)
    {
        for(i=1;i<imax;i++)
            P[i][jmax+1] = P[i][jmax];
    }
    
    return(P);
}


REAL** chorin::boundaryValuesF (REAL** U, REAL**F, int imax, int jmax, int wN, int wS, int wE, int wW)
{
    int j;
    //matrix o_matrix;
    //F = RMATRIX_ZERO(F,0.0, imax+2, jmax+2); //delete object??
    
    
    //left wall
    if(wW==2)
    {
        for(j=1;j<jmax+2;j++)
            F[0][j] = U[0][j];
    }
    //right wall
    if(wE==2)
    {
        for(j=1;j<jmax+2;j++)
            F[imax+1][j] = U[imax+1][j];
    }
    
    
    return(F);
}

REAL** chorin::boundaryValuesG (REAL** V, REAL** G, int imax, int jmax, int wN, int wS, int wE, int wW)
{
    int i;
    
    // bottom  wall
    if(wS==2)
    {
        for(i=1;i<imax+2;i++)
            G[i][0] = V[i][0];
    }
    if(wN==2)
    {
        for(i=1;i<imax+2;i++)
            G[i][jmax+1] = V[i][jmax];
    }
    return(G);
}

REAL ** chorin::RMATRIX_ZERO(REAL **m, REAL zero, int imax, int jmax)
{
       int i,j;
    
       m = new REAL*[imax];
       for(i = 0; i < imax; ++i){
           m[i] = new REAL[jmax];
           for(j=0;j<jmax; ++j){
               m[i][j]=zero;
           }
       }
       
       /*for(i = 0; i < imax; ++i){
           for(j=0; j< jmax; ++j){
               m[i][j]=zero;
           }
       }*/
       
       return(m);
       
}
