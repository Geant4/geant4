// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//tools_build_use tools inlib csz zlib


#include <tools/args>
#include <tools/fileis>

#include <tools/rroot/file>
#include <tools/rroot/rall>

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


  return EXIT_SUCCESS;
}
