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

double norm( DEVector * v, bool taxi ){
  if( !taxi )
    return blas::dot( *v, *v );
  else{
    double sum = 0;
    for(int i = 1; i <= v->length(); i++){
      double val = (*v)(i);
      if( val >= 0 ) sum += val;
      else sum += (-1*val);
    }

    return sum;
    //    for( int i = 0; i < 
  }
    
}


int lll( GEMatrix * B, double y, bool taxi ){
  int m = B->numRows();
  int n = B->numCols();

  GEMatrix Bstar(m,n);
  GEMatrix coef = gram_schmidt( B, &Bstar );

  bool done = false;

  int rounds = 1;

  while( !done ){
    //  this was put here
    //  GEMatrix coef = gram_schmidt( B, &Bstar );
    
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
      DEVector sumb = Bstar(_,_(j,j)).vectorView();
      double norma = norm( &suma, taxi );
      double normb = norm( &sumb, taxi );
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
  return rounds;
}


int wr( GEMatrix * B, GEMatrix * delta ){
  bool taxi = false;
  int m = B->numRows();
  int n = B->numCols();
  int num_replaced = 0;

  // check all (n choose 2) combinations of vertices
  for( int i = 1; i < n; i++ ){
    for( int j = i+1; j <= n; j++ ){
      
      int k;
      // e in [-1,1]
      for( int e = -1; e <= 1; e += 2 ){
	if( (*delta)(i,i) < (*delta)(j,j) )
	  k = j;
	else
	  k = i;

	DEVector V = (*B)(_,_(i,i)).vectorView()
	  + e * (*B)(_,_(j,j)).vectorView();


	// if the norm of V is different replace it!
	if( norm( &V, taxi ) < (*delta)(k,k) ){

	  // make sure no bugs
	  double dotprod = 0;
	  for( int x = 1; x <= m; x++ )  dotprod += V(x)*V(x);
	  (*delta)(k,k) = (*delta)(i,i) + (*delta)(j,j) + 2*e*(*delta)(i,j);
	  if( (*delta)(k,k) != dotprod ) return -1;

	  
	  for( int h = 1; h <= n; h++ ){
	    if (h != i && h != j){
	      (*delta)(k,h) = (*delta)(i,h)+e*(*delta)(j,h);
	      (*delta)(h,k) = (*delta)(k,h);
	    }
	  }

	  if( k != i ){
	    (*delta)(k,i) = (*delta)(i,i) + e*(*delta)(j,i);
	    (*delta)(i,k) = (*delta)(k,i);
	  }
	  else{
	    (*delta)(k,j) = (*delta)(i,j) + e*(*delta)(j,j);
	    (*delta)(j,k) = (*delta)(k,j);
	  }
	  
	  for( int l = 1; l <= m; l++ ){
	    (*B)(l,k) = V(l);
	  }

	  num_replaced++;
	}

      }

    }
  }
  return num_replaced;
}


int wr_taxi( GEMatrix * B, double * delta ){
  int n = B->numCols();
  int m = B->numRows();

  int num_replaced = 0;

  for( int i = 1; i < n; i++ ){
    for( int j = i+1; j <= n; j++ ){
      
      int k;
      // e in [-1,1]
      for( int e = -1; e <= 1; e += 2 ){
	if( delta[i-1] < delta[j-1] )
	  k = j;
	else
	  k = i;
	
	DEVector V = (*B)(_,_(i,i)).vectorView()
	  + e * (*B)(_,_(j,j)).vectorView();
	
	// if the norm of V is better replace it!
	double new_norm = norm( &V, true );
	if( new_norm < delta[k-1] ){
	  for( int l = 1; l <= m; l++ ){
	    (*B)(l,k) = V(l);
	  }
	  delta[k-1] = new_norm;
	  num_replaced++;
	}
	else{
	  //	  cout << "New = " << new_norm << " and old = " 

	}
      }
    }
  }
  return num_replaced;
}


void fill_delta( GEMatrix * B, GEMatrix * delta ){
  int m = B->numRows();
  int n = B->numCols();
  int n1 = delta->numRows();
  int n2 = delta->numCols();

  if( n == n1 && n == n2 ){
    for( int i = 1; i <= n; i++ ){
      for( int j = i; j <= n; j++ ){
	(*delta)(i,j) = blas::dot( (*B)(_,_(i,i)).vectorView(), 
				   (*B)(_,_(j,j)).vectorView() ); 
	if( i != j ) (*delta)(j,i) = (*delta)(i,j);
      }
    }
  }
  else{
    cerr << "Error: Matrices indices don't match" << endl;
  }
}


void print_matrix( GEMatrix * A, ostream * o ){
  int m = A->numRows();
  int n = A->numCols();
  for( int i = 1; i <= m; i++ ){
    for( int j = 1; j <= n; j++ ){
      *o << " ";
      if( (*A)(i,j) >= 0 ) *o << " ";
      *o << (*A)(i,j);
    }
    *o << endl;
  }

}
