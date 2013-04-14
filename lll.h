#ifndef _LLL_H
#define _LLL_H

// functions for Linear Algebra
// makes use of FLENS library

#include <flens/flens.cxx>
#include <cmath>

using namespace flens;

typedef double T;
typedef DenseVector<Array<T>> DEVector;
typedef GeMatrix<FullStorage<T, ColMajor>> GEMatrix;

//typedef GEMatrix::ConstView            GEMatrixConstView;
//typedef GEMatrix::View                 GEMatrixView;
typedef GEMatrix::VectorView GEMatrixVectorView;
typedef GEMatrix::ConstVectorView GEMatrixConstVectorView;

const Underscore<GEMatrix::IndexType> _;

const double y = .75;

/*
Runs Gram-Schmidt process on basis B, makes Bstar the orthogonal basis
Returns 'alpha', the coefficients used during the process, as a matrix
*/
GEMatrix gram_schmidt( GEMatrix * B, GEMatrix * Bstar );

void lll( GEMatrix * B );

void print_matrix( GEMatrix * A );

#endif
