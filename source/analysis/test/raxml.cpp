// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//exlib_build_use exlib inlib expat

// for example to read file produced by :
//   inlib/example/cpp/tools_test_waxml.

#include <tools/args>
#include <tools/raxml>

//#include <tools/gzip>

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
    std::cout << " obj = " << (*it).object() << std::endl;
    std::cout << " class = " << (*it).cls() << std::endl;
    std::cout << " path = " << (*it).path() << std::endl;
    std::cout << " name = " << (*it).name() << std::endl;
  }}

 {std::vector<tools::raxml_out>::const_iterator it;
  for(it=objs.begin();it!=objs.end();++it) {
    const tools::raxml_out& raxml_out = *it;
    if(raxml_out.cls()==tools::histo::h1d::s_class()) {
      tools::histo::h1d* h = (tools::histo::h1d*)raxml_out.object();
      std::cout << "---------------------------------" << std::endl;
      std::cout << "h1d : title " << h->title() << std::endl;
      std::cout << "h1d : entries " << h->entries() << std::endl;
      std::cout << "h1d : mean " << h->mean() << " rms " << h->rms() << std::endl;
    } else if(raxml_out.cls()==tools::histo::h2d::s_class()) {
      tools::histo::h2d* h = (tools::histo::h2d*)raxml_out.object();
      std::cout << "---------------------------------" << std::endl;
      std::cout << "h2d : title " << h->title() << std::endl;
      std::cout << "h2d : entries " << h->entries() << std::endl;
      std::cout << "h2d : mean_x " << h->mean_x() << " rms_x " << h->rms_x() << std::endl;
      std::cout << "h2d : mean_y " << h->mean_y() << " rms_y " << h->rms_y() << std::endl;
    } else if(raxml_out.cls()==tools::aida::ntuple::s_class()) {
      tools::aida::ntuple* nt = (tools::aida::ntuple*)raxml_out.object();
      std::cout << "---------------------------------" << std::endl;
      std::cout << "ntuple : title " << nt->title() << std::endl;
      std::cout << "ntuple : rows " << nt->rows() << std::endl;
    }
  }}

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif

  return 0;
}
