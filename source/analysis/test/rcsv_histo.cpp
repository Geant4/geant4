// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program read a out_[].csv files produced by wcsv_histo.

#ifdef TOOLS_MEM
#include <tools/mem>
#endif

#include <tools/rcsv_histo>
#include <tools/args>

#include <fstream>

#include <iostream>
#include <cstdlib>

int main(int argc,char** argv) {

#ifdef TOOLS_MEM
  tools::mem::set_check_by_class(true);{
#endif //TOOLS_MEM

  tools::args args(argc,argv);

  std::string file;
  if(!args.file(file)) {
    std::cout << "give a .csv file." << std::endl;
    return EXIT_FAILURE;
  }

  std::ifstream reader(file.c_str());
  if(reader.fail()) {
    std::cout << "can't open " << file << std::endl;
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////
  /// create a csv histo reader : //////////////////////////
  //////////////////////////////////////////////////////////

  bool verbose = args.is_arg("-verbose");

  tools::rcsv::histo rh(reader);
  std::string _class;
  void* obj;
  if(!rh.read(std::cout,_class,obj,verbose)) {
    std::cout << "can't read histo." << std::endl;
    return EXIT_FAILURE;
  }

  if(_class==tools::histo::h1d::s_class()) {
    tools::histo::h1d* h = static_cast<tools::histo::h1d*>(obj);
    h->hprint(std::cout);
    delete h;
  } else if(_class==tools::histo::h2d::s_class()) {
    tools::histo::h2d* h = static_cast<tools::histo::h2d*>(obj);
    h->hprint(std::cout);
    delete h;
  } else if(_class==tools::histo::h3d::s_class()) {
    tools::histo::h3d* h = static_cast<tools::histo::h3d*>(obj);
    h->hprint(std::cout);
    delete h;
  } else if(_class==tools::histo::p1d::s_class()) {
    tools::histo::p1d* h = static_cast<tools::histo::p1d*>(obj);
    h->hprint(std::cout);
    delete h;
  } else if(_class==tools::histo::p2d::s_class()) {
    tools::histo::p2d* h = static_cast<tools::histo::p2d*>(obj);
    h->hprint(std::cout);
    delete h;

  } else {
    std::cout << "unknown class " << _class << std::endl;
  }

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  reader.close();

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif //TOOLS_MEM

  return EXIT_SUCCESS;
}
