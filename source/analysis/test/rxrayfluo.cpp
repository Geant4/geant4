// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//tools_build_use tools inlib csz zlib


#include <tools/args>
#include <tools/rroot/file>
#include <tools/rroot/streamers>
#include <tools/rroot/fac>
#include <tools/rroot/tree>
#include <tools/rroot/ntuple>
#include <tools/ntuple_binding>

#ifdef TOOLS_DONT_HAVE_ZLIB
#else
#include <tools/gzip_buffer>
#endif

#include <iostream>
#include <cstdlib>

int main(int argc,char** argv) {


  tools::args args(argc,argv);

  bool verbose = args.is_arg("-verbose");

  std::string file = "xrayfluo.root";
  tools::rroot::file rfile(std::cout,file,verbose);
#ifdef TOOLS_DONT_HAVE_ZLIB
#else
  rfile.add_unziper('Z',tools::gunzip_buffer);
#endif

  //////////////////////////////////////////////////////
  /// read the histogram : /////////////////////////////
  //////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("h0");
  if(!key) {
    std::cout << "key for h0 not found." << std::endl;
    return EXIT_FAILURE;
  }
  unsigned int sz;
  char* buf = key->get_object_buffer(sz);
  if(!buf) {
    std::cout << "can't get data buffer for h0." << std::endl;
    return EXIT_FAILURE;
  }
  //std::cout << "size of h0 : " << sz << std::endl;

  tools::rroot::buffer b(std::cout,rfile.byte_swap(),sz,buf,key->key_length(),verbose);
  tools::histo::h1d* h = tools::rroot::TH1D_stream(b); //we get ownership on h.
  if(!h) {
    std::cout << "streaming failed for h0." << std::endl;
    return EXIT_FAILURE;
  }
  h->hprint(std::cout);
  delete h;}
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  /// read the ntuple TTree : //////////////////////////
  //////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("ntuple");
  if(!key) {
    std::cout << "key for ntuple not found." << std::endl;
    return EXIT_FAILURE;
  }
  unsigned int sz;
  char* buf = key->get_object_buffer(sz);
  if(!buf) {
    std::cout << "can't get data buffer for ntuple." << std::endl;
    return EXIT_FAILURE;
  }
  tools::rroot::buffer b(std::cout,rfile.byte_swap(),sz,buf,key->key_length(),verbose);
  tools::rroot::fac fac(rfile);
  tools::rroot::tree tree(rfile,fac);
  if(!tree.stream(b)) {
    std::cout << "TTree streaming failed." << std::endl;
    return EXIT_FAILURE;
  }
  tree.dump(std::cout,"","  ");
  tools::uint64 entries = tree.entries();  
 {for(tools::uint32 i=0;i<5;i++){
    if(!tree.show(std::cout,i)) {
      std::cout << "show failed for entry " << i << std::endl;
      return EXIT_FAILURE;
    }       
  }}
 {for(tools::uint64 i=tools::mx<tools::int64>(5,entries-5);i<entries;i++){
    if(!tree.show(std::cout,(tools::uint32)i)) {
      std::cout << "show failed for entry " << i << std::endl;
      return EXIT_FAILURE;
    }
  }}

  }
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////
  /// read "ntuple" with tools::rroot::ntuple : ////////
  //////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("ntuple");
  if(!key) {
    std::cout << "key for ntuple not found." << std::endl;
    return EXIT_FAILURE;
  }
  unsigned int sz;
  char* buf = key->get_object_buffer(sz);
  if(!buf) {
    std::cout << "can't get data buffer for ntuple." << std::endl;
    return EXIT_FAILURE;
  }
  tools::rroot::buffer b(std::cout,rfile.byte_swap(),sz,buf,key->key_length(),verbose);
  tools::rroot::fac fac(rfile);
  tools::rroot::tree tree(rfile,fac);
  if(!tree.stream(b)) {
    std::cout << "TTree streaming failed." << std::endl;
    return EXIT_FAILURE;
  }
  tools::rroot::ntuple ntu(tree); //use the flat ntuple API.
  if(!ntu.initialize(std::cout)) {
    std::cout << "can't initialize ntuple." << std::endl;
    return EXIT_FAILURE;
  }

  tools::read::icolumn<double>* col = ntu.find_column<double>("energies");
  if(!col) {
    std::cout << "can't find column energies." << std::endl;
    return EXIT_FAILURE;
  }

  tools::histo::h1d he("energies",100,-1,3);
  ntu.start();
  while(ntu.next()){
    double v;
    if(!col->get_entry(v)) {
      std::cout << "get_entry(double) failed." << std::endl;
      return EXIT_FAILURE;
    }
    //std::cout << v << std::endl;
    he.fill(v);
  }
  // Should print :
  //   ENTRIES = 22243
  //   MEAN VALUE = 0.731019
  //   R . M . S = 0.943105
  he.hprint(std::cout);}

  /////////////////////////////////////////////////////////////////////////
  /// read "ntuple" with tools::rroot::ntuple and ntuple_binding : ////////
  /////////////////////////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("ntuple");
  if(!key) {
    std::cout << "key for ntuple not found." << std::endl;
    return EXIT_FAILURE;
  }
  unsigned int sz;
  char* buf = key->get_object_buffer(sz);
  if(!buf) {
    std::cout << "can't get data buffer for ntuple." << std::endl;
    return EXIT_FAILURE;
  }
  tools::rroot::buffer b(std::cout,rfile.byte_swap(),sz,buf,key->key_length(),verbose);
  tools::rroot::fac fac(rfile);
  tools::rroot::tree tree(rfile,fac);
  if(!tree.stream(b)) {
    std::cout << "TTree streaming failed." << std::endl;
    return EXIT_FAILURE;
  }
  tools::rroot::ntuple ntu(tree); //use the flat ntuple API.
  tools::ntuple_binding nbd;
  double energy;
  nbd.add_column("energies",energy);
  if(!ntu.initialize(std::cout,nbd)) {
    std::cout << "can't initialize ntuple." << std::endl;
    return EXIT_FAILURE;
  }

  tools::histo::h1d he("energies",100,-1,3);
  ntu.start();
  while(ntu.next()){
    if(!ntu.get_row()) {
      std::cout << "get_row() failed." << std::endl;
      return EXIT_FAILURE;
    }
    //std::cout << energy << std::endl;
    he.fill(energy);
  }
  // Should print :
  //   ENTRIES = 22243
  //   MEAN VALUE = 0.731019
  //   R . M . S = 0.943105
  he.hprint(std::cout);}
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  

  return EXIT_SUCCESS;
}
