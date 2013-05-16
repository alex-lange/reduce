#include <iostream>
//#include <flens/flens.cxx>
//#include <cmath>
#include <unistd.h>
#include "lll.h"
#include "/home/alex/research/archer/g.h"

using namespace std;
//using namespace flens;


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

bool is_dom_set( g * gr, int * dom_set ){
  int n = gr->order();
  bool dominated[gr->order()];
  
  for( int i = 0; i < n; i++ ){
    dominated[i] = false;
  }

  for( int i = 0; i < n; i++ ){
    if( dom_set[i] != 0 ){
      dominated[i] = true;
      for( int j = 0; j < n; j++ ){
	if( gr->is_edge(i,j) ){
	  dominated[j] = true;
	}
      }
    }
  }

  // for( int i = 0; i < n; i++ ) cout << dominated[i] << " ";
  //cout << endl;
  
  for( int i = 0; i < n; i++ )
    if( !dominated[i] ) return false;

  return true;
  
}

bool check_vec( GEMatrix::View gr_B, int row_size, int n ){
  /*cout << "CHECK" << endl;
    cout << col << endl;*/
  bool good_one = true;
  bool still_good_one = true;
      
  for( int r = row_size; r >= row_size-n; r-- ){
    if( gr_B(r,1) != 0 ){
      good_one = false;			
      break;
    }
  }
  if( good_one ){
    bool posi = true;
    bool found_first = false;
    for( int r = 1; r <= 2*n; r++ ){
      int val = gr_B(r,1);
      if( val != 0 ){
	if( !found_first ){
	  posi = val > 0;
	  found_first = true;
	}
	else{
	  if( (posi && val < 0) || (!posi && val > 0 ) ){
	    still_good_one = false;
	    break;
	  }
	}
      }
    }
  }
  return still_good_one && good_one;
}



int main( int argc, char * argv[] ){
  char opt;
  string out_file;
  int n, row_size, col_size, found, sum_found, num_tests;
  // istream test;

  if( argc >= 3 ){
    if( argv[1][0] == '-' ){
      opt = argv[1][1];
      out_file = argv[2];
    }
  }
  else if( argc == 2 ){
    opt = 'g';
    out_file = argv[1];
  }
  else{
    cerr << "Error: Invalid args" << endl;
    return 1;
  }

  vector<string> graph6s;
  string g_string;

  string graph_file, log_file;

  graph_file = out_file + ".graphs";
  ofstream g_file( graph_file.c_str() );
  if( !g_file.is_open() ){
    cerr << "Error opening " << graph_file << endl;
    return 0;
  }
  
  log_file = out_file + ".log";
  ofstream log( log_file.c_str() );
  if( !log.is_open() ){
    cerr << "Error opening " << log_file << endl;
    return 0;
  }


  if( opt == 'g' ){
    while( getline( cin, g_string ) ){
      graph6s.push_back( g_string );
    }

    num_tests = graph6s.size();
  }
  // Hamming graph
  else if( opt == 'h' ){
    cout << "Testing LLL domination with a Hamming graph" << endl;
    int l;
    if( argc == 4 ){
      l = atoi( argv[3] );
    }
    else{
      cout << "l? ";
      cin >> l;
    }
    n = pow( 3, l );
    g ham(n);
    ham.make_hamming( l );
    graph6s.push_back( ham.to_g6() );

  }
  // Generate random Erdos-Renyi graphs
  else if( opt == 'r' ){
    float p;
        
    cout << "Testing LLL domination with random graphs..." << endl;
    cout << "n? ";
    cin >> n; 
    cout << "p? ";
    cin >> p;
    cout << "Number of tests? ";
    cin >> num_tests;
    //    cout << "Filename? ";
    //cin >> out_file;

    log << "Random graphs" << endl;
    log << n << " " << p << " " << num_tests << endl;
    

    cout << "Generating random graphs..." << endl;
    for( int i = 0; i < num_tests; i++ ){
      g gr( n );
      int added = gr.make_rand_er( p );
      graph6s.push_back( gr.to_g6() );
      sleep(2);
    }
  }

  found = 0; sum_found = 0;

  int test_num = 1;
  for( vector<string>::iterator it = graph6s.begin();
       it != graph6s.end(); it++ ){
    g gr( (*it)[0] - 63 );
    gr.read_g6( *it );
    n = gr.order();
    row_size = 3*n;
    col_size = 2*n + 1;

    gr.print_g6( &g_file );

    GEMatrix gr_A;
    GEMatrix gr_B;
    
    get_matrix( &gr, &gr_A );
    
    cout << gr_A << endl;
    
    get_dom_basis( &gr_A, &gr_B, 10 );
    
    print_matrix( &gr_B );
    
    lll( &gr_B, 0.95, false );

    cout << endl;
    
    print_matrix( &gr_B );
    
    bool success = false; bool sum_success = false;
    log << "**** " << test_num << " ****" << endl;
    log << "Max Degree = " << gr.max_degree() << endl;
    cout << "Max Degree = " << gr.max_degree() << endl;
    
    for( int c = 1; c <= col_size && !success; c++ ){

      bool still_good_one = check_vec( gr_B(_,_(c,c) ), row_size, n );

      
      if( still_good_one ){
	//	log << test_num << " FOUND ONE, column " << c << endl;;
	cout << "FOUND ONE " << c << endl;

	int dom_set[n];
	for( int r = 1; r <= n; r++ ){
	  dom_set[r-1] = gr_B(r,c);
	}
	  
	cout << "FOUND ONE " << c << endl;
	  
	// if the vector is a domination set
	if( is_dom_set( &gr, dom_set ) ){
	  success = true;
	    
	  cout << "... AND it's a dominating set" << endl;
	  log << "... AND it's a dominating set:" << endl;
	  cout << "x = ";
	  log << "x = ";
	  int dom_size = 0;
	  for( int i = 0; i < n; i++ ){
	    log << dom_set[i] << " ";
	    cout << dom_set[i] << " ";
	    if( dom_set[i] != 0 )
	      dom_size++;
	  }
	  log << endl;
	  log << "Cardinality: " << dom_size << endl;
	  
	}
	
	cout << endl;
      }

    }
    
    if( !success ){
      log << "Did not find one, attempting to sum columns..." <<endl;
      for( int c1 = 1; c1 < col_size; c1++ ){
	for( int c2 = c1+1; c2 <=col_size; c2++ ){

	  bool go = true;
	  int scale = 1;

	  while(go){

	    GEMatrix sum = gr_B(_,_(c1,c1) ) + scale * gr_B(_,_(c2,c2) );
	    bool still_good_one = check_vec( sum, row_size, n );
      
	    if( still_good_one ){
	      //	    log << test_num << " FOUND ONE, columns " << c1 << " + " << c2 << endl;;
	    
	      int dom_set[n];
	      for( int r = 1; r <= n; r++ ){
		dom_set[r-1] = sum(r,1);
	      }
	    
	      cout << "FOUND ONE " << c1 << " + " << c2 << endl;
	  
	      // if the vector is a domination set
	      if( is_dom_set( &gr, dom_set ) ){
		sum_success = true;
		
		cout << "... AND it's a dominating set" << endl;
		log << "... AND it's a dominating set:" << endl;
		cout << "x = ";
		log << "x = ";
		int dom_size = 0;
		for( int i = 0; i < n; i++ ){
		  log << dom_set[i] << " ";
		  cout << dom_set[i] << " ";
		  if( dom_set[i] != 0 )
		    dom_size++;
		}
		log << endl;
		log << "Cardinality: " << dom_size << endl;

		cout << endl;
	      }
	    }
	    if( scale == -1 )
	      go = false;
	    else
	      scale = -1;
	  }
	}
      }
    }
    
    
    if( success ){ 
      found++;
    }
    else if( sum_success ){
      sum_found++;
    }
    test_num++;
    log << endl;

  }
  log << "******************" << endl;
  log << "Original found: " << found << endl;
  log << "Sum found: " << sum_found << endl;
  log << "Total found: " << found+sum_found << " / " << num_tests << endl;
  g_file.close();
  log.close();
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

}


