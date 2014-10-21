// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//exlib_build_use exlib inlib expat

// for example to read file produced by :
//   inlib/example/cpp/waxml.cpp.

#include <tools/raxml>
#include <tools/ntuple_binding>

//#include <tools/gzip>
#include <tools/args>
#include <iostream>

int main(int argc,char** argv) {

#ifdef TOOLS_MEM
  tools::mem::set_check_by_class(true);{
#endif

  tools::args args(argc,argv);

  std::string file;
  if(!args.file(file)) {
    std::cout << " give an AIDA xml file." << std::endl;
    return EXIT_FAILURE;
  }

  bool verbose = false;

  tools::xml::default_factory fac;
  tools::raxml ml(fac,std::cout,verbose);

  //ml.set_compressed_reader(new exlib::gzip_reader());

  std::vector<tools::raxml_out>& objs = ml.objects();
  objs.clear();

  ml.load_file(file,false);

 {std::vector<tools::raxml_out>::const_iterator it;
  for(it=objs.begin();it!=objs.end();++it) {
    std::cout << "---------------------------------" << std::endl;
    std::cout << " obj = " << (*it).object() << std::endl;
    std::cout << " class = " << (*it).cls() << std::endl;
    std::cout << " path = " << (*it).path() << std::endl;
    std::cout << " name = " << (*it).name() << std::endl;
  }}

 {std::vector<tools::raxml_out>::const_iterator it;
  for(it=objs.begin();it!=objs.end();++it) {
    const tools::raxml_out& raxml_out = *it;
    std::cout << "---------------------------------" << std::endl;
    std::cout << "obj name = " << (*it).name() << std::endl;
    if(raxml_out.cls()==tools::histo::h1d::s_class()) {
      tools::histo::h1d* h = (tools::histo::h1d*)raxml_out.object();
      std::cout << "h1d : title " << h->title() << std::endl;
      std::cout << "h1d : entries " << h->entries() << std::endl;
      std::cout << "h1d : mean " << h->mean() << " rms " << h->rms() << std::endl;
    } else if(raxml_out.cls()==tools::histo::h2d::s_class()) {
      tools::histo::h2d* h = (tools::histo::h2d*)raxml_out.object();
      std::cout << "h2d : title " << h->title() << std::endl;
      std::cout << "h2d : entries " << h->entries() << std::endl;
      std::cout << "h2d : mean_x " << h->mean_x() << " rms_x " << h->rms_x() << std::endl;
      std::cout << "h2d : mean_y " << h->mean_y() << " rms_y " << h->rms_y() << std::endl;
    } else if(raxml_out.cls()==tools::aida::ntuple::s_class()) {
      tools::aida::ntuple* nt = (tools::aida::ntuple*)raxml_out.object();
      std::cout << "ntuple : title " << nt->title() << std::endl;
      std::cout << "ntuple : rows " << nt->rows() << std::endl;
     {const std::vector<tools::aida::base_col*>& cols = nt->cols();
      std::vector<tools::aida::base_col*>::const_iterator itc;
      for(itc=cols.begin();itc!=cols.end();++itc) {
	tools::aida::aida_base_col* acol = tools::safe_cast<tools::aida::base_col,tools::aida::aida_base_col>(*(*itc));
        if(acol) {
          std::cout << "ntuple : col name=" << (*itc)->name() << ", aida_type=" << acol->aida_type() << std::endl;
        } else {
          std::cout << "ntuple : col name=" << (*itc)->name() << std::endl;
        }
      }}
      //if from inlib/examples/cpp/waxml out.aida file.
     {tools::aida::aida_col<double>* col = nt->find_column<double>("rgauss");
      if(col) {
        nt->start();
        for(unsigned int row=0;row<5;row++) {
          if(!nt->next()) break;
          double v;
          if(!col->get_entry(v)) {}
          std::cout << " " << v << std::endl;
        }
      }}
     {tools::aida::aida_col<std::string>* col = nt->find_column<std::string>("strings");
      if(col) {
        nt->start();
        for(unsigned int row=0;row<5;row++) {
          if(!nt->next()) break;
	  std::string v;
          if(!col->get_entry(v)) {}
          std::cout << "row = " << row << ", string = " << v << std::endl;
        }
      }}

      ///////////////////////////////////////////////////////
      /// read by using variable column binding : ///////////
      ///////////////////////////////////////////////////////
     {tools::ntuple_binding nbd;
      double rgauss;
      nbd.add_column("rgauss",rgauss);
      std::string sval;
      nbd.add_column("strings",sval);
      if(!nt->set_binding(std::cout,nbd)) {
        std::cout << "set ntuple binding failed." << std::endl;
        return EXIT_FAILURE;
      }
      tools::histo::h1d h("rgauss",100,-2,2);
      nt->start();
      while(nt->next()){
        if(!nt->get_row()) {
          std::cout << "get_row() failed." << std::endl;
          return EXIT_FAILURE;
        }
        //std::cout << energy << std::endl;
        //std::cout << "string : " << sval << std::endl;
        h.fill(rgauss);
      }
      h.hprint(std::cout);}

    }
  }}

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif

  return 0;
}
