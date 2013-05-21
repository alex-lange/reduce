#ifndef _LLL_H
#define _LLL_H

// functions for Linear Algebra
// makes use of FLENS library

#include <flens/flens.cxx>
#include <cmath>
#include <iostream>

using namespace flens;

typedef double T;
typedef DenseVector<Array<T>> DEVector;
typedef GeMatrix<FullStorage<T, ColMajor>> GEMatrix;

typedef DEVector::IndexType IndexType;

//typedef GEMatrix::ConstView            GEMatrixConstView;
//typedef GEMatrix::View                 GEMatrixView;
typedef GEMatrix::VectorView GEMatrixVectorView;
typedef GEMatrix::ConstVectorView GEMatrixConstVectorView;

const Underscore<GEMatrix::IndexType> _;

//const double y = .75;

/*
Runs Gram-Schmidt process on basis B, makes Bstar the orthogonal basis
Returns 'alpha', the coefficients used during the process, as a matrix
*/
GEMatrix gram_schmidt( GEMatrix * B, GEMatrix * Bstar );

double norm( DEVector * v, bool taxi = false );

int lll( GEMatrix * B, double y = 0.75, bool taxi = false );
int wr( GEMatrix * B, GEMatrix * delta );
int wr_taxi( GEMatrix * B, double * delta ); 
void fill_delta( GEMatrix * B, GEMatrix * delta );
void print_matrix( GEMatrix * A, std::ostream * o = &std::cout );

#endif
