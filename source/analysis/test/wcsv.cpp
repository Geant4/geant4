// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program produces a out.csv file.
//
//  See rcsv.kumac, rcsv.C for an example of how to manipulate
// the content of this file with CERN-PAW and CERN-ROOT.

#include <tools/wcsv_ntuple>

#include <tools/randd>
#include <tools/randf>
#include <tools/sto>

#include <fstream>

#include <iostream>
#include <cstdlib>

int main(int,char**) {

  //////////////////////////////////////////////////////////
  /// create a .csv file : /////////////////////////////////
  //////////////////////////////////////////////////////////
  std::ofstream writer("out.csv");
  if(writer.fail()) {
    std::cout << "can't open out.csv." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int entries = 1000;
  tools::rgaussd rg(1,2);
  tools::rbwf rbw(0,1);

  //////////////////////////////////////////////////////////
  /// create and write a ntuple : //////////////////////////
  //////////////////////////////////////////////////////////
  //tools::wcsv::ntuple ntu(writer,'\t');
  tools::wcsv::ntuple ntu(writer); //default sep is ','

  // create some columns with basic types :
  tools::wcsv::ntuple::column<unsigned int>* col_index = ntu.create_column<unsigned int>("index");
  tools::wcsv::ntuple::column<double>* col_rgauss = ntu.create_column<double>("rgauss");
  tools::wcsv::ntuple::column<float>* col_rbw = ntu.create_column<float>("rbw");
  tools::wcsv::ntuple::column<std::string>* col_str = ntu.create_column<std::string>("strings");

  //ntu.write_hippo_header();

  // fill :
  for(unsigned int count=0;count<entries;count++) {    
    col_index->fill(count);
    col_rgauss->fill(rg.shoot());
    col_rbw->fill(rbw.shoot());
    col_str->fill("str "+tools::to(count));
    ntu.add_row(); // it will write columns data as a row in the file.
  }

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  writer.close();

  return EXIT_SUCCESS;
}
