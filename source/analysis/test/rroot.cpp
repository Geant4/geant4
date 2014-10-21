// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//tools_build_use tools inlib csz zlib


#include <tools/args>
#include <tools/fileis>

#include <tools/rroot/file>
#include <tools/rroot/rall>

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

  std::string file;
  if(!args.file(file)) {
    std::cout << " give a root file." << std::endl;
    return EXIT_FAILURE;
  }

  bool verbose = args.is_arg("-verbose");
  bool ls = args.is_arg("-ls");
  bool dump = args.is_arg("-dump");

 {bool is;
  tools::file::is_root(file,is);
  if(!is) {
    std::cout << " file is not a root file." << std::endl;
    return EXIT_FAILURE;
  }}

  tools::rroot::file rfile(std::cout,file,verbose);
#ifdef TOOLS_DONT_HAVE_ZLIB
#else
  rfile.add_unziper('Z',tools::gunzip_buffer);
#endif

  if(ls) {
    std::cout << "format version " << rfile.version() << std::endl;
  }
      
  const std::vector<tools::rroot::key*>& keys = rfile.dir().keys();
  tools::rroot::read(std::cout,rfile,keys,true,ls,dump,0);

  ///////////////////////////////////////////////////////////////
  /// if reading the out.root produced with wroot.cpp : /////////
  ///////////////////////////////////////////////////////////////
 {tools::rroot::TDirectory* dir = tools::rroot::find_dir(rfile.dir(),"histo");
  if(dir) {
   {tools::rroot::key* key = dir->find_key("rg");
    if(key) {
      tools::histo::h1d* h = tools::rroot::key_to_h1d(*key);
      if(h) {
        std::cout << "h1d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean " << h->mean() << ", rms " << h->rms()
                  << std::endl;
        delete h;
      }
    }}
   {tools::rroot::key* key = dir->find_key("rf");
    if(key) {
      tools::histo::h1d* h = tools::rroot::key_to_h1d(*key);
      if(h) {
        std::cout << "h1d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean " << h->mean() << ", rms " << h->rms()
                  << std::endl;
        delete h;
      }
    }}
   {tools::rroot::key* key = dir->find_key("rgbw");
    if(key) {
      tools::histo::h2d* h = tools::rroot::key_to_h2d(*key);
      if(h) {
        std::cout << "h2d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean_x " << h->mean_x() << ", rms_x " << h->rms_x()
                  << ", mean_y " << h->mean_y() << ", rms_y " << h->rms_y()
                  << std::endl;
        delete h;
      }
    }}
   {tools::rroot::key* key = dir->find_key("prof");
    if(key) {
      tools::histo::p1d* h = tools::rroot::key_to_p1d(*key);
      if(h) {
        std::cout << "p1d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean " << h->mean() << ", rms " << h->rms()
                  << std::endl;
        delete h;
      }
    }}
   {tools::rroot::key* key = dir->find_key("prof2D");
    if(key) {
      tools::histo::p2d* h = tools::rroot::key_to_p2d(*key);
      if(h) {
        std::cout << "p2d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean_x " << h->mean_x() << ", rms_x " << h->rms_x()
                  << ", mean_y " << h->mean_y() << ", rms_y " << h->rms_y()
                  << std::endl;
        delete h;
      }
    }}
   {tools::rroot::key* key = dir->find_key("rggbw");
    if(key) {
      tools::histo::h3d* h = tools::rroot::key_to_h3d(*key);
      if(h) {
        std::cout << "h3d : " << h->title()
                  << ", all_entries " << h->all_entries()
                  << ", entries " << h->entries()
                  << ", mean_x " << h->mean_x() << ", rms_x " << h->rms_x()
                  << ", mean_y " << h->mean_y() << ", rms_y " << h->rms_y()
                  << ", mean_z " << h->mean_z() << ", rms_z " << h->rms_z()
                  << std::endl;
        delete h;
      }
    }}
  }}
  // read an ntuple :
 {tools::rroot::key* key = rfile.dir().find_key("rg_rbw");
  if(key) {
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
    //tools::uint64 entries = tree.entries();
   {for(tools::uint32 i=0;i<5;i++){
      if(!tree.show(std::cout,i)) {
        std::cout << "show failed for entry " << i << std::endl;
        return EXIT_FAILURE;
      }
    }}

    // read with the flat ntuple API :
   {tools::rroot::ntuple ntu(tree); //use the flat ntuple API.
    tools::ntuple_binding nbd;
    double v_rgauss;
    nbd.add_column("rgauss",v_rgauss);
    std::string v_string;
    nbd.add_column("strings",v_string);
    if(!ntu.initialize(std::cout,nbd)) {
      std::cout << "can't initialize ntuple with ntuple_binding." << std::endl;
      return EXIT_FAILURE;
    }
    tools::histo::h1d hg("rgauss",100,-5,5);
    ntu.start();
    unsigned int count = 0;
    while(ntu.next()){
      if(!ntu.get_row()) {
        std::cout << "get_row() failed." << std::endl;
        return EXIT_FAILURE;
      }
      hg.fill(v_rgauss);
      if(count<5) std::cout << "v_string " << v_string << std::endl;
      count++;
    }
    std::cout << "ntuple_binding(rgauss) : " << hg.mean() << " " << hg.rms() << std::endl;}

  }}

  ///////////////////////////////////////////////////////////////
  /// if reading the pawdemo.root : /////////////////////////////
  ///////////////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("h10");
  if(key) {
    tools::histo::h1d* h = tools::rroot::key_to_h1d(*key);
    if(h) {
      std::cout << "h1d : h10"
                << ", all_entries " << h->all_entries()
                << ", entries " << h->entries()
                << ", mean " << h->mean() << ", rms " << h->rms()
                << std::endl;
      delete h;
    }
  }}

  /////////////////////////////////////////////////////////////////////
  /// if reading the prof.root produced with croot_TProfile.cpp : /////
  /////////////////////////////////////////////////////////////////////
 {tools::rroot::key* key = rfile.dir().find_key("prof");
  if(key) {
    tools::histo::p1d* h = tools::rroot::key_to_p1d(*key);
    if(h) {
      std::cout << "p1d : prof"
                << ", all_entries " << h->all_entries()
                << ", entries " << h->entries()
                << ", mean " << h->mean() << ", rms " << h->rms()
                << std::endl;
      delete h;
    }
  }}
 {tools::rroot::key* key = rfile.dir().find_key("prof2D");
  if(key) {
    tools::histo::p2d* h = tools::rroot::key_to_p2d(*key);
    if(h) {
      std::cout << "p2d : prof"
                << ", all_entries " << h->all_entries()
                << ", entries " << h->entries()
                << ", mean_x " << h->mean_x() << ", rms_x " << h->rms_x()
                << ", mean_y " << h->mean_y() << ", rms_y " << h->rms_y()
                << std::endl;
      delete h;
    }
  }}


  return EXIT_SUCCESS;
}
