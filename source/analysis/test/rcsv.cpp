// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program read a out.csv file.

#include <tools/rcsv_ntuple>
#include <tools/scast>
#include <tools/args>
#include <tools/path>
#include <tools/histo/h1d>

#include <iostream>
#include <cstdlib>

int main(int argc,char** argv) {

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
  /// hippodraw flavour ? //////////////////////////////////
  /// - one header line for the ntuple title. //////////////
  /// - one csv line for column names. /////////////////////
  /// - data at csv format. ////////////////////////////////
  //////////////////////////////////////////////////////////
  bool is_hippo = false;
  if(args.is_arg("-hippo")) {
    is_hippo = true;
  } else {
    if((tools::suffix(file)=="hiptxt")||(tools::suffix(file)=="tnt")) {
      is_hippo = true;
    } else {
      is_hippo = tools::rcsv::ntuple::is_hippo(std::cout,reader);
    }
  }

  if(is_hippo) std::cout << "hippodraw file." << std::endl;

  //////////////////////////////////////////////////////////
  /// create a csv ntuple reader : /////////////////////////
  //////////////////////////////////////////////////////////

  bool verbose = args.is_arg("-verbose");

  char sep = 0; //0 = guessed from first data line.

  unsigned int isep = 0;
  if(args.find("-sep",isep)) sep = isep;

  tools::rcsv::ntuple ntu(reader);
  ntu.set_hippo(is_hippo);
  if(!ntu.initialize(std::cout,sep,"x",verbose)) { //col suffix is x.
    std::cout << "can't initialize ntuple." << std::endl;
    return EXIT_FAILURE;
  }

  if(args.is_arg("-cols")){
    ntu.dump_columns(std::cout);
    return EXIT_SUCCESS;
  }

  std::string scol;
  if(args.find("-col",scol)) {

   {typedef tools::read::icolumn<double> cold_t;
    cold_t* col = ntu.find_column<double>(scol);
    if(col) {
      ntu.start();
      while(ntu.next()){
        double v;
        if(!col->get_entry(v)) {
          std::cout << "get_entry(double) failed." << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << v << std::endl;
      }
      reader.close();
      return EXIT_SUCCESS;
    }}

   {typedef tools::read::icolumn<std::string> cols_t;
    cols_t* col = ntu.find_column<std::string>(scol);
    if(col) {
      ntu.start();
      while(ntu.next()){
	std::string v;
        if(!col->get_entry(v)) {
          std::cout << "get_entry(std::string) failed." << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << v << std::endl;
      }
      reader.close();
      return EXIT_SUCCESS;
    }}

    std::cout << "column " << scol << " not found." << std::endl;
    return EXIT_FAILURE;

  } else if(args.is_arg("-bind")) {

    // read with binding :
    tools::ntuple_binding nbd;
    double v_rgauss;
    nbd.add_column("rgauss",v_rgauss);
    std::string v_string;
    nbd.add_column("strings",v_string);
    if(!ntu.initialize(std::cout,nbd)) {
      std::cout << "can't initialize ntuple with ntuple_binding." << std::endl;
      return EXIT_FAILURE;
    }
    tools::histo::h1d h("rgauss",100,-5,5);
    ntu.start();
    unsigned int count = 0;
    while(ntu.next()){
      if(!ntu.get_row()) {
        std::cout << "get_row() failed." << std::endl;
        return EXIT_FAILURE;
      }
      h.fill(v_rgauss);
      if(count<5) std::cout << "v_string " << v_string << std::endl;
      count++;
    }
    std::cout << " h " << h.mean() << " " << h.rms() << std::endl;

  } else { // read all

    typedef tools::read::icol icol_t;
    typedef tools::read::icolumn<double> cold_t;
    typedef tools::read::icolumn<tools::csv_time> colt_t;
    typedef tools::read::icolumn<std::string> cols_t;
    const std::vector<icol_t*>& cols = ntu.columns();
    std::vector<icol_t*>::const_iterator it;
  
    ntu.start();
    while(ntu.next()){
      for(it=cols.begin();it!=cols.end();++it) {
        if(it!=cols.begin()) std::cout << " ";
        if(cold_t* cold = tools::id_cast<icol_t,cold_t>(*(*it))) {
          double v;
          if(!cold->get_entry(v)) {
            std::cout << "get_entry(double) failed." << std::endl;
            return EXIT_FAILURE;
          }
          std::cout << v;
        } else if(colt_t* colt = tools::id_cast<icol_t,colt_t>(*(*it))) {
          tools::csv_time v;
          if(!colt->get_entry(v)) {
            std::cout << "get_entry(time_t) failed." << std::endl;
            return EXIT_FAILURE;
          }
          std::cout << v.m_l;
        } else if(cols_t* colst = tools::id_cast<icol_t,cols_t>(*(*it))) {
          std::string v;
          if(!colst->get_entry(v)) {
            std::cout << "get_entry(string) failed." << std::endl;
            return EXIT_FAILURE;
          }
          std::cout << v;
        } else {
          std::cout << "column cast failed." << std::endl;
          return EXIT_FAILURE;
        }
      }
      std::cout << std::endl;
    }

  }

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  reader.close();

  return EXIT_SUCCESS;
}
