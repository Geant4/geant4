// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#include <tools/waxml/begend>
#include <tools/waxml/histos>
#include <tools/waxml/ntuple>

#include <tools/random>

#include <iostream>
#include <cstdlib>

int main(int,char**) {

#ifdef TOOLS_MEM
  tools::mem::set_check_by_class(true);{
#endif

  //////////////////////////////////////////////////////////
  /// create a .aida file : ////////////////////////////////
  //////////////////////////////////////////////////////////
  std::ofstream writer("out.aida");
  if(writer.fail()) {
    std::cout << "can't open out.aida." << std::endl;
    return EXIT_FAILURE;
  }
  tools::waxml::begin(writer);

  //////////////////////////////////////////////////////////
  /// create and write some histos : ///////////////////////
  //////////////////////////////////////////////////////////
  unsigned int entries = 1000000;
  tools::random::gauss rg(1,2);
  tools::random::bw rbw(0,1);

 {tools::histo::h1d h("Gauss",100,-5,5);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),1.4);
  }
  if(!tools::waxml::write(writer,h,"/histo","rg")) {
    std::cout << "can't write h1d." << std::endl;
    return EXIT_FAILURE;
  }}

 {tools::histo::p1d h("Profile",100,-5,5,-2,2);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),rbw.shoot(),1);
  }
  if(!tools::waxml::write(writer,h,"/histo","prof")) {
    std::cout << "can't write prof." << std::endl;
    return EXIT_FAILURE;
  }}

 {tools::histo::h2d h("Gauss_BW",20,-5,5,20,-2,2);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),rbw.shoot(),0.8);
  }
  if(!tools::waxml::write(writer,h,"/histo","rgbw")) {
    std::cout << "can't write h2d." << std::endl;
    return EXIT_FAILURE;
  }}

  //////////////////////////////////////////////////////////
  /// create and write a flat ntuple : /////////////////////
  //////////////////////////////////////////////////////////
 {tools::waxml::ntuple ntu(writer);

  // create some columns with basic types :
  tools::waxml::ntuple::column<double>* col_rgauss =
    ntu.create_column<double>("rgauss");
  tools::waxml::ntuple::column<double>* col_rbw =
    ntu.create_column<double>("rbw");

  ntu.write_header("/tuple","rg_rbw","Randoms");

  // fill :
  for(unsigned int count=0;count<10000;count++) {    
    col_rgauss->fill(rg.shoot());
    col_rbw->fill(rbw.shoot());
    ntu.add_row(); // it will write columns data as a <row> in the file.
  }
  ntu.write_trailer();}

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  tools::waxml::end(writer);
  writer.close();

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif

  return EXIT_SUCCESS;
}
