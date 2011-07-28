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
//   PAW > cd histo
//   PAW > h/plot 10
//   PAW > h/plot 20
//   PAW > h/plot 30
//   PAW > cd //lun**/ntuple
//   PAW > ntuple/print 40
//   PAW > ntuple/plot 40.rg
//   PAW > ntuple/plot 40.rbw

#include <tools/hbook/wfile>
#include <tools/hbook/h1>
#include <tools/hbook/h2>
#include <tools/hbook/p1>
#include <tools/hbook/wntuple>

#include <tools/randf>
#include <tools/cmp>

#ifdef WIN32
extern "C" int __stdcall SETPAWC();
extern "C" int __stdcall SETNTUC();
#else
extern "C" int setpawc_();
extern "C" int setntuc_();
#endif

#include <iostream>

int main(int a_argc,char** a_argv) {

 {// Initialize HBOOK :
#ifdef WIN32
  tools::hbook::CHLIMIT(SETPAWC());
  SETNTUC(); //for ntuple.
#else
  tools::hbook::CHLIMIT(setpawc_());
  setntuc_(); //for ntuple.
#endif
  tools::hbook::CHCDIR("//PAWC"," ");

  unsigned int unit = 1;
  tools::hbook::wfile hfile(std::cout,"out.hbook",unit);
  if(!hfile.is_valid()) {
    std::cout << "main : open of out.hbook failed." << std::endl;
    return EXIT_FAILURE;
  }
  // At this point, in HBOOK, we should have :
  //   - created a //LUN1 directory attached to the file
  //   - created a //PAWC/LUN1 in memory
  //   - be in the directory //PAWC/LUN1.

  // create an "histo" HBOOK directory both in memory and in the file :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  tools::hbook::CHMDIR("histo"," ");
  tools::hbook::CHCDIR("//LUN1"," ");
  tools::hbook::CHMDIR("histo"," ");
  tools::hbook::CHCDIR("//PAWC/LUN1/histo"," ");
  // the five upper lines could have been done with :
  //hfile.cd_home();
  //hfile.mkcd("histo");

  unsigned int entries = 100000;
  tools::randf::gauss rg(0,1);
  tools::randf::bw rbw(0,1);

  // create some histos under //PAWC/LUN1/histo.
  // At creation, an tools::hbook histo keeps track of the HBOOK
  // directory path in which it had been created (its "booking directory").
  // This permits to have multiple HBOOK histos with same integer
  // id but be in different HBOOK directories.
  tools::hbook::h1 h1(10,"Gauss",100,-5,5);
  tools::hbook::h2 h2(20,"Gauss_BW",100,-5,5,100,-2,2);
  tools::hbook::p1 p(30,"Profile",100,-5,5,-2,2);
  // We should be in //PAWC/LUN1/histo.

  // Return to //PAWC/LUN1 :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");

  // At each fill() an tools::hbook histo saves the current position,
  // goes into its booking directory, do the CHFILL and
  // then returns to the saved position.
 {for(unsigned int count=0;count<entries;count++) {
    h1.fill(rg.shoot());
  }}
 {for(unsigned int count=0;count<entries;count++) {
    h2.fill(rg.shoot(),rbw.shoot());
  }}
 {for(unsigned int count=0;count<entries;count++) {
    p.fill(rg.shoot(),rbw.shoot());
  }}

  // the uppers fill() do a CHCDIR to the booking directory.
  // As this may be costly we have introduced :
  //   fill_beg(), fill_fast(), fill_end()
  // that permits to change directories once.

 {h1.fill_beg();
  for(unsigned int count=0;count<entries;count++) {
    h1.fill_fast(rg.shoot()); //no CHCDIR is done here.
  }h1.fill_end();}

 {h2.fill_beg();
  for(unsigned int count=0;count<entries;count++) {
    h2.fill_fast(rg.shoot(),rbw.shoot());
  }h2.fill_end();}

 {p.fill_beg();
  for(unsigned int count=0;count<entries;count++) {
    p.fill_fast(rg.shoot(),rbw.shoot());
  }p.fill_end();}

  // We should be in //PAWC/LUN1.

  //std::cout << " bins " << h.axis().bins()
  //          << ", entries " << h.all_entries()
  //          << ", bin[50] " << h.bin_height(50)
  //          << std::endl;
  //std::cout << " mean " << h.mean() << ", rms " << h.rms() << std::endl;
  tools::cmp<double>(std::cout,h1.mean(),0,0.01);
  tools::cmp<double>(std::cout,h1.rms(),1,0.01);

  // create an "ntuple" directory both in memory and in the file :
  hfile.cd_home();      //go under //PAWC/LUN1
  //hfile.cd_up();      //could have been ok too.
  hfile.mkcd("ntuple"); //create  //LUN1/ntuple and //PAWC/LUN1/ntuple
  // We should be under //PAWC/LUN1/ntuple

  // create a ntuple under //LUN1/ntuple. Contrary to an tools::hbook
  // histo, an tools::hbook ntuple is attached to a file from creation.
  // At creation, the tools::hbook::wntuple keeps track of two booking
  // directories : the one in memory and the one on file
  // (here //PAWC/LUN1/ntuple and //LUN1/ntuple).
  tools::hbook::wntuple* ntu = new tools::hbook::wntuple(40,"An ntuple");

 {tools::hbook::wntuple::column<double>* col_rg = 
    ntu->create_column<double>("rg");
  tools::hbook::wntuple::column<float>* col_rbw = 
    ntu->create_column<float>("rbw");
  tools::hbook::wntuple::column<int>* col_count = 
    ntu->create_column<int>("count");

  // We should be under //PAWC/LUN1/ntuple.

  // Return to //PAWC/LUN1 :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");

  // At each add_row() the tools::hbook::wntuple saves the current
  // position, goes into its memory and file booking directory
  // (two CHCDIR are done), writes the row with CHFNT on file
  // and then returns to the saved position.

  // first way to loop :
 {for(unsigned int count=0;count<1000;count++) {
    col_rg->fill(rg.shoot());    
    col_rbw->fill(rbw.shoot());    
    col_count->fill(count);    
    ntu->add_row();
  }}

  // the upper add_row() does two CHCDIR to position in the
  // memory and file booking directories. As this may
  // be costly we have introduced :
  //   add_row_beg(), add_row_fast(), add_row_end()
  // that permits to change directories once.
 {ntu->add_row_beg(); //go under //PAWC/LUN1/ntuple and //LUN1/ntuple
  for(unsigned int count=0;count<1000;count++) {
    col_rg->fill(rg.shoot());    
    col_rbw->fill(rbw.shoot());    
    col_count->fill(count);    
    ntu->add_row_fast();  //no CHCDIR is done here.
  }
  ntu->add_row_end();}} //returns to //PAWC/LUN1

  // fill by finding columns (very costly) :
 {ntu->add_row_beg();
  for(unsigned int count=0;count<1000;count++) {
    tools::hbook::wntuple::column<double>* col_rg = 
      ntu->find_column<double>("rg");
    tools::hbook::wntuple::column<float>* col_rbw = 
      ntu->find_column<float>("rbw");
    tools::hbook::wntuple::column<int>* col_count = 
      ntu->find_column<int>("count");

    col_rg->fill(rg.shoot());    
    col_rbw->fill(rbw.shoot());    
    col_count->fill(count);    

    ntu->add_row_fast();
  }
  ntu->add_row_end();}

  // We should be under //PAWC/LUN1

  // The below changes directory under //PAWC/LUN1 and //LUN1
  // and writes memory data (histos) into the file :
  if(!hfile.write()) {
    std::cout << "main : wfile.write() failed." << std::endl;
    return EXIT_FAILURE;
  }

  delete ntu; //WARNING : have to delete the ntuple before closing the file.

  // close the file and delete the //LUN1 directory.
  if(!hfile.close()) {
    std::cout << "main : wfile.close() failed." << std::endl;
    return EXIT_FAILURE;
  }
  
  // The destructor of wfile deletes the directory //PAWC/LUN1.

  }

  return EXIT_SUCCESS;
}
