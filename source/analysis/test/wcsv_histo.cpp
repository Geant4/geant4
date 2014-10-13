// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program produces a out.root file.
//
//  See rroot.C for an example of how to manipulate
// (and check !) the content of this file with CERN-ROOT.

#ifdef TOOLS_MEM
#include <tools/mem>
#endif //TOOLS_MEM

#include <tools/wcsv_histo>

#include <tools/histo/h1d>
#include <tools/histo/h2d>
#include <tools/histo/h3d>
#include <tools/histo/p1d>
#include <tools/histo/p2d>

#include <tools/histo/h1df>
#include <tools/histo/h2df>

#include <tools/randd>
#include <tools/randf>

#include <tools/args>

#include <fstream>

#include <iostream>
#include <cstdlib>

int wcsv_histo(int argc,char** argv) {

#ifdef TOOLS_MEM
  tools::mem::set_check_by_class(true);{
#endif //TOOLS_MEM
  tools::args args(argc,argv);

  bool verbose = args.is_arg("-verbose");

  //////////////////////////////////////////////////////////
  /// create some histos : /////////////////////////////////
  //////////////////////////////////////////////////////////

  unsigned int entries = 1000000;

  tools::rgaussd rg(1,2);
  tools::rbwd rbw(0,1);

  char hc = '#';
  char sep = ',';

 {std::ofstream writer("out_h1d.csv");
  if(writer.fail()) {std::cout << "can't open out_h1d.csv." << std::endl;return EXIT_FAILURE;}
  tools::histo::h1d h("Rand gauss",100,-5,5);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),1.4);
  // plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"1.4*entries");
  h.add_annotation("empty","");
  if(verbose) {
    std::cout << "h1d : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean " << h.mean() << ", rms " << h.rms()
              << std::endl;
  }
  // write :
  if(!tools::wcsv::hto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_h1df.csv");
  if(writer.fail()) {std::cout << "can't open out_h1df.csv." << std::endl;return EXIT_FAILURE;}
  tools::rgaussf rf(1,2);
  tools::histo::h1df h("GaussF",100,-5,5);
  for(unsigned int count=0;count<entries;count++) h.fill(rf.shoot(),1.4F);
  // plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"1.4*entries");
  if(verbose) {
    std::cout << "h1df : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean " << h.mean() << ", rms " << h.rms()
              << std::endl;
  }
  // write :
  if(!tools::wcsv::hto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_p1d.csv");
  if(writer.fail()) {std::cout << "can't open out_p1d.csv." << std::endl;return EXIT_FAILURE;}
  tools::histo::p1d h("Profile",100,-5,5,-2,2);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),rbw.shoot(),1);
  if(verbose) {
    std::cout << "p1d : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean " << h.mean() << ", rms " << h.rms()
              << std::endl;
  }
  if(!tools::wcsv::pto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_h2d.csv");
  if(writer.fail()) {std::cout << "can't open out_h2d.csv." << std::endl;return EXIT_FAILURE;}
  tools::histo::h2d h("Gauss_BW",20,-5,5,20,-2,2);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),rbw.shoot(),0.8);
  //plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"rand bw");
  h.add_annotation(tools::histo::key_axis_z_title(),"0.8*entries");
  if(verbose) {
    std::cout << "h2d : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean_x " << h.mean_x() << ", rms_x " << h.rms_x()
              << ", mean_y " << h.mean_y() << ", rms_y " << h.rms_y()
              << std::endl;
  }
  if(!tools::wcsv::hto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_p2d.csv");
  if(writer.fail()) {std::cout << "can't open out_p2d.csv." << std::endl;return EXIT_FAILURE;}
  tools::histo::p2d h("Profile2D",100,-5,5,100,-5,5,-2,2);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),rg.shoot(),rbw.shoot(),1);
  if(verbose) {
    std::cout << "p2d : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean_x " << h.mean_x() << ", rms_x " << h.rms_x()
              << ", mean_y " << h.mean_y() << ", rms_y " << h.rms_y()
              << std::endl;
  }
  if(!tools::wcsv::pto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_h3d.csv");
  if(writer.fail()) {std::cout << "can't open out_h3d.csv." << std::endl;return EXIT_FAILURE;}
  tools::histo::h3d h("Gauss_Gauss_BW",20,-5,5,20,-5,5,20,-2,2);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),rg.shoot(),rbw.shoot(),0.8);
  //plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_z_title(),"rand bw");
  if(verbose) {
    std::cout << "h3d : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean_x " << h.mean_x() << ", rms_x " << h.rms_x()
              << ", mean_y " << h.mean_y() << ", rms_y " << h.rms_y()
              << ", mean_z " << h.mean_z() << ", rms_z " << h.rms_z()
              << std::endl;
  }
  if(!tools::wcsv::hto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

 {std::ofstream writer("out_h1d_edges.csv");
  if(writer.fail()) {std::cout << "can't open out_h1d.csv." << std::endl;return EXIT_FAILURE;}
  std::vector<double> edges;
  double width = (5.0-(-5.0))/100;
  for(unsigned int index=0;index<=100;index++) edges.push_back(-5+index*width);
  tools::histo::h1d h("Gauss",edges);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),1.4);
  // plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"1.4*entries");
  if(verbose) {
    std::cout << "h1d : edges : " << h.title()
              << ", all entries " << h.all_entries()
              << ", entries " << h.entries()
              << ", mean " << h.mean() << ", rms " << h.rms()
              << std::endl;
  }
  // write :
  if(!tools::wcsv::hto(writer,h.s_cls(),h,sep,hc)) return EXIT_FAILURE;
  writer.close();}

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif //TOOLS_MEM

  return EXIT_SUCCESS;
}

#ifndef __CLING__
int main(int argc,char** argv) {return wcsv_histo(argc,argv);}
#endif
