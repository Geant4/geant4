// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//tools_build_use tools inlib zlib

//  This program produces a out.root file.
//
//  See rroot.C for an example of how to manipulate
// (and check !) the content of this file with CERN-ROOT.


#include <tools/wroot/file>
#include <tools/wroot/to>

#include <tools/histo/h1d>
#include <tools/histo/h2d>
#include <tools/histo/p1d>
#include <tools/wroot/ntuple>

#include <tools/randd>
#include <tools/randf>

#ifdef TOOLS_DONT_HAVE_ZLIB
#else
#include <tools/gzip_buffer>
#endif

#include <tools/args>
#include <iostream>

#include <ctime>
#include <cstdlib>


int main(int argc,char** argv) {

  tools::args args(argc,argv);

  //////////////////////////////////////////////////////////
  /// create a .root file : ////////////////////////////////
  //////////////////////////////////////////////////////////
  bool verbose = args.is_arg("-verbose");
  std::string file = "out.root";
  tools::wroot::file rfile(std::cout,file,verbose);
#ifdef TOOLS_DONT_HAVE_ZLIB
#else
  if(args.is_arg("-noz")){
  } else {
    rfile.add_ziper('Z',tools::gzip_buffer);
    rfile.set_compression(1);
  }
#endif

  tools::wroot::directory* dir = rfile.dir().mkdir("histo");
  if(!dir) {
    std::cout << "can't create diectory." << std::endl;
    return EXIT_FAILURE;
  }

  //tools::wroot::directory* dxxx = dir->mkdir("xxx");
  //if(!dxxx) {
  //  std::cout << "can't create diectory." << std::endl;
  //  return EXIT_FAILURE;
  //}

  unsigned int entries = 10000;

  int nofColumns = 1000;
  std::cout << "nofColumns = " << nofColumns << std::endl;

  //////////////////////////////////////////////////////////
  /// create and fill a ntuple : ///////////////////////////
  //////////////////////////////////////////////////////////
 {//WARNING : the ntuple can't be on the stack. It is owned
  //          by the directory.
  tools::wroot::ntuple* ntu = 
    new tools::wroot::ntuple(rfile.dir(),"rg_rbw","Randoms");

  std::vector<tools::wroot::ntuple::column<double>* > columns;
  for (int i=0; i<nofColumns; ++i) {
    std::stringstream ss;
    ss << "rgauss_";
    ss << i;
    std::string str = ss.str();
    columns.push_back(ntu->create_column<double>(str));
  }  

  if(args.is_arg("-large")){
    entries = 300000000; //to test >2Gbytes file.
  }

  ntu->set_basket_size(1000000);

  std::cout << "Fill ntuple start ..." << std::endl;
  time_t start1 = time(NULL);
  tools::rgaussd rg(1,2);
  for(unsigned int count=0;count<entries;count++) {    
    double value = rg.shoot();
    for (int i=0; i<nofColumns; ++i) {
      if(!columns[i]->fill(value)) {
        std::cout << "col_rgauss fill failed." << std::endl;
        break;
      }
    }  
    if(!ntu->add_row()) {
      std::cout << "ntuple fill failed." << std::endl;
      break;
    }
  }
  std::cout << "Fill ntuple finished ..." << std::endl;
  printf("%.2f\n", (double)(time(NULL) - start1));
  }

  //////////////////////////////////////////////////////////
  /// write and close file : ///////////////////////////////
  //////////////////////////////////////////////////////////
 {unsigned int n;  
  std::cout << "Write ntuple start ..." << std::endl;
  time_t start2 = time(NULL);
  if(!rfile.write(n)) {
    std::cout << "file write failed." << std::endl;
  }
  std::cout << "Write ntuple finished ..." << std::endl;
  printf("%.2f\n", (double)(time(NULL) - start2));
  }
  
  rfile.close();

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


  return EXIT_SUCCESS;
}
