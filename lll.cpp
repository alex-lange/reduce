#include "lll.h"

using namespace std;
using namespace flens;

/*
Runs Gram-Schmidt process on basis B, makes Bstar the orthogonal basis
Returns 'alpha', the coefficients used during the process, as a matrix
*/
GEMatrix gram_schmidt( GEMatrix * B, GEMatrix * Bstar ){
  int m = B->numRows();
  int n = B->numCols();
  GEMatrix alpha(n,n);

  (*Bstar)(_,_(1,1)) = (*B)(_,_(1,1));

  for( int j = 2; j <= n; j++ ){
    (*Bstar)(_,_(j,j)) = (*B)(_,_(j,j));
    for( int i = 1; i < j; i++ ){
      // computes coefficient of projection
      double alph = (*Bstar)(_,_(i,i)).vectorView()*(*B)(_,_(j,j)).vectorView()
	/  ((*Bstar)(_,_(i,i)).vectorView() * (*Bstar)(_,_(i,i)).vectorView());

      alpha(i,j) = alph;
      
      // replaces vector with orthogonal vector
      (*Bstar)(_,_(j,j)) = (*Bstar)(_,_(j,j)) - alph * (*Bstar)(_,_(i,i));
    }
  }
  return alpha;
}


void lll( GEMatrix * B ){
  int m = B->numRows();
  int n = B->numCols();

  GEMatrix Bstar(m,n);
  GEMatrix coef = gram_schmidt( B, &Bstar );

  bool done = false;

  int rounds = 0;

  while( !done ){
    
    // step 1
    for( int j = 2; j <= n; j++ ){
      for( int i = j-1; i >= 1; i-- ){
	if( abs(coef(i,j)) > .5 ){
	  (*B)(_,_(j,j)) = (*B)(_,_(j,j)) - 
	    floor(coef(i,j)+.5) * (*B)(_,_(i,i));
	}
      }
    }

    // step 2
    coef = gram_schmidt( B, &Bstar );
    bool found = false;

    for( int j = 1; j < n; j++ ){
      DEVector suma =Bstar(_,_(j+1,j+1)).vectorView() + 
	coef(j,j+1)*Bstar(_,_(j,j)).vectorView();
      double norma = blas::dot( suma, suma );
      double normb = blas::dot(Bstar(_,_(j,j)).vectorView(),
			       Bstar(_,_(j,j)).vectorView());
      if( norma < y * normb ){
	found = true;
	blas::swap( (*B)(_,_(j,j)).vectorView(), 
		    (*B)(_,_(j+1,j+1)).vectorView());
	break;
      }
    }
    if( !found ){
      done = true;
    }
    else{
      coef = gram_schmidt( B, &Bstar );
      rounds++;
    }
    if( rounds % 100 == 0 ){
      cout << rounds << endl;
      //      print_matrix( B );
    }
  }
}


void print_matrix( GEMatrix * A ){
  int m = A->numRows();
  int n = A->numCols();
  for( int i = 1; i <= m; i++ ){
    for( int j = 1; j <= n; j++ ){
      cout << " ";
      if( (*A)(i,j) >= 0 ) cout << " ";
      cout << (*A)(i,j);
    }
    cout << endl;
  }

}
