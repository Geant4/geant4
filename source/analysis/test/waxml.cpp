// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#include <tools/waxml/begend>
#include <tools/waxml/histos>
#include <tools/waxml/ntuple>

#include <tools/randd>

#include <iostream>
#include <cstdlib>

int main(int,char**) {

#ifdef TOOLS_MEM
  tools::mem::set_check_by_class(true);{
#endif

  //////////////////////////////////////////////////////////
  /// create a .aida file : ////////////////////////////////
  //////////////////////////////////////////////////////////
  std::ofstream writer("out.aida");
  if(writer.fail()) {
    std::cout << "can't open out.aida." << std::endl;
    return EXIT_FAILURE;
  }
  tools::waxml::begin(writer);

  //////////////////////////////////////////////////////////
  /// create and write some histos : ///////////////////////
  //////////////////////////////////////////////////////////
  tools::rgaussd rg(1,2);
  tools::rbwd rbw(0,1);
  unsigned int entries = 1000000;

 {tools::histo::h1d h("Gauss",100,-5,5);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),1.4);
  // plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"1.4*entries");
  // write :
  if(!tools::waxml::write(writer,h,"/histo","rg")) {
    std::cout << "can't write h1d." << std::endl;
    return EXIT_FAILURE;
  }}

 {tools::histo::p1d h("Profile",100,-5,5,-2,2);
  for(unsigned int count=0;count<entries;count++) h.fill(rg.shoot(),rbw.shoot(),1);
  if(!tools::waxml::write(writer,h,"/histo","prof")) {
    std::cout << "can't write prof." << std::endl;
    return EXIT_FAILURE;
  }}

 {std::string title = "Gauss_BW";  
  // have XML special characters in the title.
  title += " lower <";
  title += " greater >";
  title += " amp &";
  title += " quote '";
  title += " double quote \"";
  tools::histo::h2d h(title,20,-5,5,20,-2,2);
  for(unsigned int count=0;count<entries;count++) {
    h.fill(rg.shoot(),rbw.shoot(),0.8);
  }
  //plotting hints :
  h.add_annotation(tools::histo::key_axis_x_title(),"rand gauss");
  h.add_annotation(tools::histo::key_axis_y_title(),"rand bw");
  h.add_annotation(tools::histo::key_axis_z_title(),"0.8*entries");
  // write :
  if(!tools::waxml::write(writer,h,"/histo","rgbw")) {
    std::cout << "can't write h2d." << std::endl;
    return EXIT_FAILURE;
  }}

  //////////////////////////////////////////////////////////
  /// create and write a flat ntuple : /////////////////////
  //////////////////////////////////////////////////////////
 {tools::waxml::ntuple ntu(writer);

  // create some columns with basic types :
  tools::waxml::ntuple::column<double>* col_rgauss = ntu.create_column<double>("rgauss");
  tools::waxml::ntuple::column<double>* col_rbw = ntu.create_column<double>("rbw");
  tools::waxml::ntuple::column<std::string>* col_str = ntu.create_column<std::string>("strings");

  ntu.write_header("/tuple","rg_rbw","Randoms");

  // fill :
  for(unsigned int count=0;count<10000;count++) {    
    col_rgauss->fill(rg.shoot());
    col_rbw->fill(rbw.shoot());
    col_str->fill("str "+tools::to(count));
    ntu.add_row(); // it will write columns data as a <row> in the file.
  }
  ntu.write_trailer();}

  ////////////////////////////////////////////////////////////////
  /// create and write a flat ntuple by using ntuple_booking : ///
  ////////////////////////////////////////////////////////////////
 {tools::ntuple_booking nbk;
  nbk.add_column<double>("rgauss");
  nbk.add_column<double>("rbw");
  nbk.add_column<std::string>("strings");
  std::vector<double> user_vec_d;
  nbk.add_column<double>("vec_d",user_vec_d);

  tools::waxml::ntuple ntu(writer,std::cout,nbk);
  if(ntu.columns().size()) {

    tools::waxml::ntuple::column<double>* col_rgauss = ntu.find_column<double>("rgauss");
    tools::waxml::ntuple::column<double>* col_rbw = ntu.find_column<double>("rbw");
    tools::waxml::ntuple::column<std::string>* col_str = ntu.find_column<std::string>("strings");

    ntu.write_header("/tuple","rg_rbw_2","Randoms");

    // fill :
    for(unsigned int count=0;count<100;count++) {    
      col_rgauss->fill(rg.shoot());
      col_rbw->fill(rbw.shoot());
      col_str->fill("str_"+tools::to(count));
     {user_vec_d.clear();
      unsigned int number = count%5;
      for(unsigned int i=0;i<number;i++) {
        user_vec_d.push_back(rbw.shoot());
      }}
      ntu.add_row(); // it will write columns data as a <row> in the file.
    }

    ntu.write_trailer();
  }}

  //////////////////////////////////////////////////////////
  /// close file : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  tools::waxml::end(writer);
  writer.close();

#ifdef TOOLS_MEM
  }tools::mem::balance(std::cout);
#endif

  return EXIT_SUCCESS;
}
