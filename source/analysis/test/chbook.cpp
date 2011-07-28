// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

//  This program produces a out.hbook file with
// an histo 1D, an histo 2D and a profile 1D.
//  With CERN-PAW, you can open the file and plot
// the histos with :
//   PAW > cd
//   PAW > h/file 101 out.hbook
//   PAW > cd
//   PAW > h/list T
//   PAW > h/plot 10
//   PAW > h/plot 20
//   PAW > h/plot 30

#include <tools/hbook/CHBOOK>

#ifdef WIN32
extern "C" int __stdcall SETPAWC();
extern "C" void __stdcall OCLOSE(int*);
#else
extern "C" int setpawc_();
extern "C" void oclose_(int*);
#endif

#include <tools/randf>

#include <iostream>

int main(int a_argc,char** a_argv) {

 {

  // initialize HBOOK :
#ifdef WIN32
  tools::hbook::CHLIMIT(SETPAWC());
#else
  tools::hbook::CHLIMIT(setpawc_());
#endif
  tools::hbook::CHCDIR("//PAWC"," ");

  // Create a file for writing. We attach it to the "LUN1" directory.
  int unit = 1;
  int ier = tools::hbook::CHROPEN(unit,"LUN1","out.hbook","NQ",1024);
  if(ier) {
    std::cout << "main :"
              << " error on hropen, code " << ier 
              << std::endl;
    return EXIT_FAILURE;
  }

  // Have a "histo" directory in memory :
  tools::hbook::CHCDIR("//PAWC"," ");
  tools::hbook::CHMDIR("histo"," ");
  tools::hbook::CHCDIR("//PAWC/histo"," ");
  // Then we are under the directory //PAWC/histo.

  // create some histos in memory :
  unsigned int entries = 100000;
  tools::randf::gauss rg(0,1);
  tools::randf::bw rbw(0,1);

  tools::hbook::CHBOOK1(10,"Gauss",100,-5,5);
 {for(unsigned int count=0;count<entries;count++) {
    tools::hbook::CHFILL(10,rg.shoot(),0,1);
  }}
  tools::hbook::CHBOOK2(20,"Gauss_BW",100,-5,5,100,-2,2);
 {for(unsigned int count=0;count<entries;count++) {
    tools::hbook::CHFILL(20,rg.shoot(),rbw.shoot(),1);
  }}
  tools::hbook::CHBPROF(30,"Profile",100,-5,5,-2,2,"");
 {for(unsigned int count=0;count<entries;count++) {
    tools::hbook::CHFILL(30,rg.shoot(),rbw.shoot(),1);
  }}

  // Write directory in file :
  // We go under :
  //   //PAWC/histo in memory
  //   //LUN1 in file
  tools::hbook::CHCDIR("//PAWC/histo"," ");
  tools::hbook::CHCDIR("//LUN1"," ");
  // Do the IO :
  tools::hbook::CHROUT(0,0,"T");
  tools::hbook::CHCDIR("//PAWC/histo"," "); //return to //PAWC/histo

  // close file :
  tools::hbook::CHREND("LUN1");
#ifdef WIN32
  OCLOSE(&unit);
#else
  oclose_(&unit);
#endif
  // remove histo directory in memory :
  tools::hbook::CHCDIR("//PAWC"," ");
  tools::hbook::CHDDIR("histo");
  
  }

  return EXIT_SUCCESS;
}
