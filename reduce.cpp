#include <iostream>
//#include <flens/flens.cxx>
//#include <cmath>

#include "lll.h"
#include "/home/alex/research/archer/g.h"

using namespace std;
//using namespace flens;


/*
Runs Gram-Schmidt process on basis B, makes Bstar the orthogonal basis
Returns 'alpha', the coefficients used during the process, as a matrix

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
  }*/


void get_matrix( g * gr, GEMatrix * gr_A ){
  int n = gr->order();
  gr_A->resize(n,n);
  
  for( int i = 0; i < n-1; i++ ){
    for( int j = i+1; j < n; j++ ){
      if( gr->is_edge(i,j) ){
	(*gr_A)(i+1,j+1)=1;
	(*gr_A)(j+1,i+1)=1;
      }
    }
  }
}

void get_dom_basis( GEMatrix * gr_A, GEMatrix * gr_B, int scale = 2 ){
  int n = gr_A->numRows();
  gr_B->resize( n + n + n, n + n + 1);

  
  for( int i = 1; i <= 2*n; i++ ){
    (*gr_B)(i,i) = 1;
  }

  for( int i = 2*n + 1; i <= 3*n; i++ ){
    for( int j = 1; j <= n; j++ ){
      if( i-2*n == j ){
	(*gr_B)(i,j) = 1 * scale;
	(*gr_B)(i,j+n) = -1 * scale;
      }
      else
	(*gr_B)(i,j) = (*gr_A)(i-2*n,j) * scale;
    }
    (*gr_B)(i,n+n+1) = -1 * scale;
  }
}


int main( int argc, char * argv[] ){
  char opt;

  if( argc == 2 ){
    if( argv[1][0] == '-' ){
      opt = argv[1][1];

    }
  }

  if( opt == 'r' ){
    int n;
    float p;
    int num_tests;
    cout << "Testing LLL domination with random graphs..." << endl;
    cout << "n? ";
    cin >> n;
    cout << "p? ";
    cin >> p;
    cout << "Number of tests? ";
    cin >> num_tests;

    cout << n << " " << p << " " << num_tests << endl;

  }
}


void test( ){

  GEMatrix B(4,4);
  GEMatrix Bstar(4,4);

  DEVector x(3), y(3);

  B = 1, 2, 3,  4,
      5, 6, 7,  8,
      9, 8, 7,  6,
      5, 4, 3, 20;

  cout << B << endl;
  cout << endl;

  /*  GEMatrix alpha = gram_schmidt( &B, &Bstar );

  cout << Bstar << endl;
  cout << endl;
  cout << alpha << endl;

  blas::swap( Bstar(_,_(1,1)).vectorView(), Bstar(_,_(2,2)).vectorView() );
  cout << endl;
  cout << Bstar << endl; */

  lll(&B);
  
  cout << B << endl;

  g gr(13);
  gr.add_edge(0,3);
  gr.add_edge(1,3);
  gr.add_edge(2,3);
  gr.add_edge(0,4);
  gr.add_edge(0,5);
  gr.add_edge(0,6);
  gr.add_edge(1,7);
  gr.add_edge(1,8);
  gr.add_edge(1,9);
  gr.add_edge(2,10);
  gr.add_edge(2,11);
  gr.add_edge(2,12);
  
  GEMatrix gr_A;
  GEMatrix gr_B;

  get_matrix( &gr, &gr_A );

  cout << gr_A << endl;

  get_dom_basis( &gr_A, &gr_B );

  print_matrix( &gr_B );

  lll( &gr_B );

  print_matrix( &gr_B );
}


