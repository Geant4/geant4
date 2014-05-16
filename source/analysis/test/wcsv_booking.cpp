// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program produces a out.csv file.
//
//  See rcsv.kumac, rcsv.C for an example of how to manipulate
// the content of this file with CERN-PAW and CERN-ROOT.

#include <tools/wcsv_ntuple>

#include <tools/randd>
#include <tools/randf>

#include <fstream>

#include <iostream>
#include <cstdlib>

int main(int,char**) {

  //////////////////////////////////////////////////////////
  /// create a .csv file : /////////////////////////////////
  //////////////////////////////////////////////////////////
  std::ofstream writer("out_bkg.csv");
  if(writer.fail()) {
    std::cout << "can't open out_bkg.csv." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int entries = 1000;
  tools::rgaussd rg(1,2);
  tools::rbwf rbw(0,1);

  //////////////////////////////////////////////////////////
  /// create a ntuple_booking object : /////////////////////
  //////////////////////////////////////////////////////////
  tools::ntuple_booking nbk;
  nbk.add_column<unsigned int>("index");
  nbk.add_column<double>("rgauss");
  nbk.add_column<float>("rbw");
  //nbk.add_column<bool>("not_handled");

  //////////////////////////////////////////////////////////
  /// create and write a ntuple : //////////////////////////
  //////////////////////////////////////////////////////////
  tools::wcsv::ntuple ntu(writer,std::cout,nbk); //default sep is ','

  if(ntu.columns().size()) {

    tools::wcsv::ntuple::column<unsigned int>* col_index =
      ntu.find_column<unsigned int>("index");
    tools::wcsv::ntuple::column<double>* col_rgauss =
      ntu.find_column<double>("rgauss");
    tools::wcsv::ntuple::column<float>* col_rbw =
      ntu.find_column<float>("rbw");

    // fill :
    for(unsigned int count=0;count<entries;count++) {    
      col_index->fill(count);
      col_rgauss->fill(rg.shoot());
      col_rbw->fill(rbw.shoot());
      ntu.add_row(); // it will write columns data as a row in the file.
    }

  }

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  writer.close();

  return EXIT_SUCCESS;
}
