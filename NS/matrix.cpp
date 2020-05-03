//
//  matrix.cpp
//  NS
//
//  Created by Diego Andrade on 4/28/20.
//  Copyright © 2020 Diego Andrade. All rights reserved.
//

#include "matrix.hpp"


matrix::matrix()
{
    
}

matrix::~matrix()
{
    
}

REAL ** matrix::RMATRIX (REAL** m, REAL v, int imax, int jmax)
{
    
    int i,j;
    
    m = new REAL*[imax];
    for(i = 0; i < imax; ++i){
        m[i] = new REAL[jmax];
        for(j=0;j<jmax; ++j){
            m[i][j]=v;
        }
    }

    return(m);
}

REAL ** matrix::RMATRIX_ZERO(double **m, double zero, int imax, int jmax)
{
    int i,j;
    
    for(i = 0; i < imax; ++i){
        for(j=0; j< jmax; ++j){
            m[i][j]=zero;
        }
    }
    
    
    return(m);
    
}

void matrix::FREE_RMATRIX(REAL** m, int imax)
{
    int i;
    for(i = 0; i < imax; ++i) {
        delete [] m[i];
    }
    
    delete [] m;
}

void matrix::PRINT_MATRIX(REAL** m, int imax, int jmax, char* name)
{
   
    int i,j;
    if(imax>=20)
    {
        cout << "\nFirst 20 Rows and COlumns of the Matrix "<< name << endl << endl;
           for( i=0;i<20;i++){
               for( j=0;j<20;j++){;
               cout << m[i][j] << "\t";
               }
           cout<<"\n";
           }
    }else{
        cout << "\nMatrix " << name << endl << endl;
        for( i=0;i<imax;i++){
            for( j=0;j<jmax;j++){;
            cout << m[i][j] << "\t";
            }
        cout<<"\n";
        
    }
    }
    
}
