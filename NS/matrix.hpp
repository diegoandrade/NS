//
//  matrix.hpp
//  NS
//
//  Created by Diego Andrade on 4/28/20.
//  Copyright © 2020 Diego Andrade. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>
#include <iostream>
#include"datadef.h"
#include <fstream>

using namespace std;

class matrix{
    
public:
    
    matrix();
    ~matrix();
    
    REAL ** m;
    
    REAL** RMATRIX (REAL** m, REAL v, int imax, int jmax);
    REAL** RMATRIX_ZERO (REAL** m, REAL zero, int imax, int jmax);
    void FREE_RMATRIX(REAL** m, int imax);
    void PRINT_MATRIX(REAL** m, int imax, int jmax, char* name);
    
    void RMATRIX_TO_FILE (REAL **m, char* file, int imax, int jmax);
   
};



#endif /* matrix_hpp */
