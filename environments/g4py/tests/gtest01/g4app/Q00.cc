// $Id: Q00.cc,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   Q00.cc
//
//                                         2005 Q
// ====================================================================
#include "MyApplication.hh"

// ====================================================================
//     main
// ====================================================================
int main(int argc, char** argv) 
{
  MyApplication* myapp= new MyApplication(argc, argv, "Q00");

  myapp-> Configure();

  if(argc==1) {
    myapp-> StartSession();
  } else {
    myapp-> ExecuteBatch();
  }

  delete myapp;

}

