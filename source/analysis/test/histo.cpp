// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#include <tools/histo/h1d>
#include <tools/histo/h2d>
#include <tools/histo/p1d>
#include <tools/histo/sliced>

#include <tools/random>

#include <iostream>

int main(int argc,char** argv) {
  bool print = true;
  if(argc==2) {
    std::string s = argv[1];
    if(s=="-noprint") print = false;
  }

  unsigned int entries = 1000000;
 {
   tools::random::gauss rg(1,2);
   tools::histo::h1d h("Gauss",100,-5,5);
   for(unsigned int count=0;count<entries;count++) {
     h.fill(rg.shoot(),1.4);
   }
   if(print) h.hprint(std::cout);
   //std::cout << " mean " << h.mean() << ", rms " << h.rms() << std::endl;

   /*
  {unsigned int nbin = h.axis().bins();
   double mn = h.axis().lower_edge();
   double mx = h.axis().upper_edge();}
   */

 }

 {
   tools::random::bw rbw(0,1);
   tools::histo::h1d h("BW",100,-5,5);
   for(unsigned int count=0;count<entries;count++) {
     h.fill(rbw.shoot(),2.3);
   }
   if(print) h.hprint(std::cout);
 }

 {
   tools::random::gauss rg(1,2);
   tools::random::bw rbw(0,1);
   tools::histo::p1d h("Profile",100,-5,5,-2,2);
   for(unsigned int count=0;count<entries;count++) {
     h.fill(rg.shoot(),rbw.shoot(),1);
   }
   if(print) h.hprint(std::cout);
 }

 {
   tools::random::gauss rg(1,2);
   tools::random::bw rbw(0,1);
   tools::histo::h2d h("Gauss_BW",100,-5,5,100,-2,2);
   for(unsigned int count=0;count<entries;count++) {
     h.fill(rg.shoot(),rbw.shoot(),0.8);
   }
   if(print) h.hprint(std::cout);

   /*
  {unsigned int nbin = h.axis_x().bins();
   double mn = h.axis_x().lower_edge();
   double mx = h.axis_x().upper_edge();}

  {unsigned int nbin = h.axis_y().bins();
   double mn = h.axis_y().lower_edge();
   double mx = h.axis_y().upper_edge();}
   */

  {
    tools::histo::h1d* projection = tools::histo::projection_x(h,"SliceX");
    if(!projection) return -1;
    projection->set_title("Gauss_BW_projectionX");
    if(print) projection->hprint(std::cout);
    delete projection;
  }
   
  {
    tools::histo::h1d* projection = tools::histo::projection_y(h,"SliceY");
    if(!projection) return -1;
    projection->set_title("Gauss_BW_projectionY");
    if(print) projection->hprint(std::cout);
    delete projection;
  }

  {
    tools::histo::h1d* slice = tools::histo::slice_x(h,40,60,"SliceX");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceX");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

  {
    tools::histo::h1d* slice =
      tools::histo::slice_x(h,
              tools::histo::axis<double>::UNDERFLOW_BIN,
              tools::histo::axis<double>::UNDERFLOW_BIN,"SliceX");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceX_UNDER");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

  {
    tools::histo::h1d* slice = 
      tools::histo::slice_x(h,
              tools::histo::axis<double>::OVERFLOW_BIN,
              tools::histo::axis<double>::OVERFLOW_BIN,"SliceX");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceX_OVER");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

  {
    tools::histo::h1d* slice = tools::histo::slice_y(h,30,50,"SliceY");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceY");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

  {
    using namespace tools::histo; //playing with namespaces.
    h1d* slice = slice_y(h,
                         axis<double>::UNDERFLOW_BIN,
                         axis<double>::UNDERFLOW_BIN,"SliceY");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceY_UNDER");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

  {
    namespace tools = tools::histo; //playing with namespaces.
    tools::h1d* slice = slice_y(h,
				tools::axis<double>::OVERFLOW_BIN,
				tools::axis<double>::OVERFLOW_BIN,"SliceY");
    if(!slice) return -1;
    slice->set_title("Gauss_BW_sliceY_OVER");
    if(print) slice->hprint(std::cout);
    delete slice;
  }

 }

  return 0;
}
