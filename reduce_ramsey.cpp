#include <iostream>
//#include <flens/flens.cxx>
//#include <cmath>

/*#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif*/

#include <unistd.h>
#include "lll.h"
#include "/home/alex/research/archer/g.h"

using namespace std;
//using namespace flens;


void get_matrix( g * gr, GEMatrix * gr_A ){
  int e = gr->num_edges();
  int t = gr->num_tris();
  
  int ** tris = gr->get_tris_array();
  gr_A->resize(t,e);
  for( int r = 1; r <= t; r++ ){
    (*gr_A)(r,tris[r-1][0]) = 1;
    (*gr_A)(r,tris[r-1][1]) = 1;
    (*gr_A)(r,tris[r-1][2]) = 1;
  }
}

void get_ram_basis( GEMatrix * gr_A, GEMatrix * gr_B, int scale = 2 ){
  int t = gr_A->numRows();
  int e = gr_A->numCols();

  gr_B->resize( e + e + t, e + e + 1);

  
  for( int i = 1; i <= 2*e; i++ ){
    (*gr_B)(i,i) = 1;
  }

  for( int i = 2*e + 1; i <= 2*e+t; i++ ){
    for( int j = 1; j <= e; j++ ){
      if( i-2*e == j ){
	(*gr_B)(i,j) = 1 * scale;
	(*gr_B)(i,j+e) = -1 * scale;
      }
      else
	(*gr_B)(i,j) = (*gr_A)(i-2*e,j) * scale;
    }
    (*gr_B)(i,e+e+1) = -1 * scale;
  }
}



bool check_vec( GEMatrix::View gr_B, int row_size, int e, int t ){
  /*cout << "CHECK" << endl;
    cout << col << endl;*/
  bool good_one = true;
  bool still_good_one = true;
      
  for( int r = row_size; r > row_size-t; r-- ){
    if( gr_B(r,1) != 0 ){
      good_one = false;			
      break;
    }
  }
  if( good_one ){
    bool posi = true;
    bool found_first = false;
    for( int r = 1; r <= 2*e; r++ ){
      int val = gr_B(r,1);
      if( val != 1 && val != -1 && val != 0 ){
	still_good_one = false;
	break;
      }
      if( val != 0 ){
	if( !found_first ){
	  posi = val > 0;
	  found_first = true;
	}
	else{
	  if( (posi && val != 1) || (!posi && val != -1 ) ){
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
  int n, row_size, col_size, found, sum_found, sum3_found;
  int no_tri, num_tests;
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

  string graph_file, log_file, fail_file;

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

  fail_file = out_file + ".fail";
  ofstream failed( fail_file.c_str() );
  if( !failed.is_open() ){
    cerr << "Error opening " << failed << endl;
    return 0;
  }


  if( opt == 'g' ){
    while( getline( cin, g_string ) ){
      graph6s.push_back( g_string );
    }

    num_tests = graph6s.size();
    cout << "Number of graphs = " << num_tests << endl;
  }
  

  found = 0; sum_found = 0, sum3_found = 0; no_tri = 0; 

  int test_num = 1;
  for( vector<string>::iterator it = graph6s.begin();
       it != graph6s.end(); it++ ){
    g gr( (*it)[0] - 63 );
    gr.read_g6( *it );
    n = gr.order();
    int e = gr.num_edges();
    int t = gr.num_tris();
    gr.print_g6( &g_file );
    bool success = false; bool sum_success = false; bool sum3_success = false;
    bool no_tris = false;
    log << "**** " << test_num << " ****" << endl;
    if( t > 0 ){
      row_size = e+e+t;
      col_size = e+e+1;

      GEMatrix gr_A;
      GEMatrix gr_B;
    
      get_matrix( &gr, &gr_A );
      get_ram_basis( &gr_A, &gr_B, 2 );
    
      print_matrix( &gr_A );
      print_matrix( &gr_B );

      lll( &gr_B, 0.85, false );
    
      print_matrix( &gr_B );
        
      for( int c = 1; c <= col_size && !success; c++ ){

	bool still_good_one = check_vec( gr_B(_,_(c,c) ), row_size, e, t );

      
	if( still_good_one ){
	  //	  log << test_num << " FOUND ONE, column " << c << endl;;
	  cout << "FOUND ONE " << c << endl;
		
	  int color[e];
	  for( int r = 1; r <= e; r++ ){
	    color[r-1] = gr_B(r,c);
	  }
	
	  if( gr.check_coloring( color ) ){
	    success = true;
	    log << test_num << " FOUND ONE, column " << c << endl;
	    log << "... AND it's a valid coloring" << endl;
	    cout << "... AND it's a valid coloring" << endl;
		    
	    for( int r = 0; r < e; r++ ){
	      log << color[r] << " ";
	      cout << color[r] << " ";
	    }
	    cout << endl;
	    log << endl;

	  }
	}
      }
    
      if( !success ){
	log << "Did not find one, attempting to sum columns..." <<endl;
	for( int c1 = 1; c1 < col_size && !sum_success; c1++ ){
	  for( int c2 = c1+1; c2 <=col_size && !sum_success; c2++ ){
	    bool go = true;
	    int scale = 1;
	    while( go ){
	      GEMatrix sum = gr_B(_,_(c1,c1) ) + scale * gr_B(_,_(c2,c2) );
	      bool still_good_one = check_vec( sum, row_size, e, t );
	  
	      if( still_good_one ){
		//		log << test_num << " FOUND ONE, columns " << c1 
		//    << " + " << scale << " * " << c2 << endl;
		cout << test_num << " FOUND ONE, columns " << c1 
		     << " + " << scale << " * " << c2 << endl;
		int color[e];
		for( int r = 1; r <= e; r++ ){
		  color[r-1] = sum(r,1)*sum(r,1);
		}
	    
		if( gr.check_coloring( color ) ){
		  sum_success = true;
		  log << test_num << " FOUND ONE, columns " << c1 
		    << " + " << scale << " * " << c2 << endl;
		  cout << "... AND it's a valid coloring" << endl;
		  log << "... AND it's a valid coloring" << endl;
		  for( int r = 0; r < e; r++ ){
		    log << color[r] << " ";
		    cout << color[r] << " ";
		  }
		  cout << endl;
		  log << endl;
	
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
      if( !success && !sum_success ){
	log << "Did not find with sum, attempting to sum three columns..." <<endl;
	cout << "Did not find with sum, attempting to sum three columns..." <<endl;
	for( int c1 = 1; c1 < col_size-1 && !sum3_success; c1++ ){
	  for( int c2 = c1+1; c2 < col_size && !sum3_success; c2++ ){
	    for( int c3 = c2+1; c3 <= col_size && !sum3_success; c3++ ){
	      bool go = true;
	      int scale1 = 1;
	      int scale2 = 1;
	      while( go ){
		GEMatrix sum = gr_B(_,_(c1,c1) ) + scale1 * gr_B(_,_(c2,c2) )
		  + scale2 * gr_B(_,_(c3,c3) );
		bool still_good_one = check_vec( sum, row_size, e, t );
	  
		if( still_good_one ){
		  /*	  log << test_num << " FOUND ONE, columns "  << c1
		      << " + " << scale1 << " * " << c2
		      << " + " << scale2 << " * " << c3 << endl;*/
		  cout<< test_num << " FOUND ONE, columns "  << c1
		      << " + " << scale1 << " * " << c2
		      << " + " << scale2 << " * " << c3 << endl;
		  int color[e];
		  for( int r = 1; r <= e; r++ ){
		    color[r-1] = sum(r,1)*sum(r,1);
		  }
	    
		  if( gr.check_coloring( color ) ){
		    log << test_num << " FOUND ONE, columns "  << c1
		      << " + " << scale1 << " * " << c2
		      << " + " << scale2 << " * " << c3 << endl;
		    sum3_success = true;
		    cout << "... AND it's a valid coloring" << endl;
		    log << "... AND it's a valid coloring" << endl;
		    for( int r = 0; r < e; r++ ){
		      log << color[r] << " ";
		      cout << color[r] << " ";
		    }
		    cout << endl;
		    log << endl;
	
		  }
		}
		if( scale1 == -1 && scale2 == -1)
		  go = false;
		else if( scale1 == -1 ){
		  scale2 = -1;
		  scale1 = 1;
		}
		else if( scale2 == -1 ){
		  scale1 = -1;
		}
		else if( scale1 == 1 ){
		  scale1 = -1;
		}
	      }
	    }
	  }
	}
      }
    }
    else{
      log << "No triangles in graph, success?" << endl;
      no_tris = true;
      no_tri++;
    }
  
    if( success ){ 
      found++;
    }
    else if( sum_success ){
      sum_found++;
    }
    else if( sum3_success ){
      sum3_found++;
    }
    else if( !no_tris ){
      failed << *it << endl;
    }
    test_num++;
    log << endl;
    
  }
  log << "******************" << endl;
  log << "Original found: " << found << endl;
  log << "Sum2 found: " << sum_found << endl;
  log << "Sum3 found: " << sum3_found << endl;
  log << "With no triangles " << no_tri << endl;
  int total = found + sum_found + sum3_found;
  log << "Total found (actual): " << total << " / " << num_tests - no_tri << endl;
  log << "Total found: " << total << " / " << num_tests << endl;
  g_file.close();
  log.close();
  failed.close();
}




void test( ){

  GEMatrix B(4,4);
  GEMatrix Bstar(4,4);
  flens::transpose(B);
  //lapack::svd();

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


