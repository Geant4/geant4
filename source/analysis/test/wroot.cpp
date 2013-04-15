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

#include <tools/random>
#include <tools/randf>

#ifdef TOOLS_DONT_HAVE_ZLIB
#else
#include <tools/gzip_buffer>
#endif

#include <tools/args>
#include <iostream>

int main(int argc,char** argv) {

  tools::args args(argc,argv);

  //////////////////////////////////////////////////////////
  /// create a .root file : ////////////////////////////////
  //////////////////////////////////////////////////////////
  std::string file = "out.root";
  tools::wroot::file rfile(std::cout,file);
#ifdef TOOLS_DONT_HAVE_ZLIB
#else
  if(args.is_arg("-noz")){
  } else {
    rfile.add_ziper('Z',tools::gzip_buffer);
    rfile.set_compression(9);
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

  //////////////////////////////////////////////////////////
  /// create some histos : /////////////////////////////////
  //////////////////////////////////////////////////////////

  unsigned int entries = 1000000;

  tools::random::gauss rg(1,2);
  tools::random::bw rbw(0,1);

 {tools::histo::h1d h("Gauss",100,-5,5);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),1.4);
  }
  // plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"1.4*entries");
  // write :
  if(!tools::wroot::to(*dir,h,"rg")) return EXIT_FAILURE;}

 {tools::histo::p1d h("Profile",100,-5,5,-2,2);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),rbw.shoot(),1);
  }
  if(!tools::wroot::to(*dir,h,"prof")) return EXIT_FAILURE;}

 {tools::histo::h2d h("Gauss_BW",20,-5,5,20,-2,2);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),rbw.shoot(),0.8);
  }
  //plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"rand bw");
  h.add_annotation(tools::histo::key_axis_z_title(),"0.8*entries");
  // write :
  if(!tools::wroot::to(*dir,h,"rgbw")) return EXIT_FAILURE;}

  //////////////////////////////////////////////////////////
  /// create and fill a ntuple : ///////////////////////////
  //////////////////////////////////////////////////////////
 {//WARNING : the ntuple can't be on the stack. It is owned
  //          by the directory.
  tools::wroot::ntuple* ntu = 
    new tools::wroot::ntuple(rfile.dir(),"rg_rbw","Randoms");
  tools::wroot::ntuple::column<int>* col_index =
    ntu->create_column<int>("index");
  tools::wroot::ntuple::column<double>* col_rgauss =
    ntu->create_column<double>("rgauss");
  tools::wroot::ntuple::column<float>* col_rbw =
    ntu->create_column<float>("rbw");

  if(args.is_arg("-large")){
    entries = 300000000; //to test >2Gbytes file.
    ntu->set_basket_size(1000000);
  }

  tools::randf::bw rbwf(0,1);
  for(unsigned int count=0;count<entries;count++) {    
    if(!col_index->fill(count)) {
      std::cout << "col_index fill failed." << std::endl;
      break;
    }
    if(!col_rgauss->fill(rg.shoot())) {
      std::cout << "col_rgauss fill failed." << std::endl;
      break;
    }
    if(!col_rbw->fill(rbwf.shoot())) {
      std::cout << "col_rbw fill failed." << std::endl;
      break;
    }
    if(!ntu->add_row()) {
      std::cout << "ntuple fill failed." << std::endl;
      break;
    }
  }}

  //////////////////////////////////////////////////////////
  /// create a ntuple from a ntuple_booking object. ////////
  //////////////////////////////////////////////////////////
 {tools::ntuple_booking nbk;
  nbk.m_name = "rg_rbw_2";
  nbk.m_title = "Randoms";
  nbk.add_column<double>("rgauss");
  nbk.add_column<float>("rbw");
  //nbk.add_column<bool>("not_handled");

  tools::wroot::ntuple* ntu = new tools::wroot::ntuple(rfile.dir(),nbk);
  if(ntu->columns().size()) {

    tools::wroot::ntuple::column<double>* col_rgauss =
      ntu->find_column<double>("rgauss");
    tools::wroot::ntuple::column<float>* col_rbw =
      ntu->find_column<float>("rbw");

    tools::randf::bw rbwf(0,1);
    for(unsigned int count=0;count<1000;count++) {    
      if(!col_rgauss->fill(rg.shoot())) {
        std::cout << "col_rgauss fill failed." << std::endl;
        break;
      }
      if(!col_rbw->fill(rbwf.shoot())) {
        std::cout << "col_rbw fill failed." << std::endl;
        break;
      }
      if(!ntu->add_row()) {
        std::cout << "ntuple fill failed." << std::endl;
        break;
      }
    }
  }}

  //////////////////////////////////////////////////////////
  /// consistency check : create an empty ntuple : /////////
  //////////////////////////////////////////////////////////
 {tools::wroot::ntuple* ntu = 
    new tools::wroot::ntuple(rfile.dir(),"empty","empty");
  ntu->create_column<int>("empty");}

  //////////////////////////////////////////////////////////
  /// write and close file : ///////////////////////////////
  //////////////////////////////////////////////////////////
 {unsigned int n;
  if(!rfile.write(n)) {
    std::cout << "file write failed." << std::endl;
  }}
  
  rfile.close();

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


  return EXIT_SUCCESS;
}
