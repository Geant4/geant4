//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: test13.cc 86624 2014-11-14 08:59:40Z kmura $
// ====================================================================
//   test13.cc
//
//   test for a list argument/return
//
//                                         2005 Q
// ====================================================================
#include <iostream>

int* alloc_int() 
{
  int* intlist= new int [10];
  for (int i=0; i<10; i++) intlist[i]= 0;

  return intlist;
}


void operate_list(int vec[10])
{
  for(int i=0; i<10; i++) {
    std::cout << vec[i] << std::endl;
    vec[i]++;
  }
}

// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;


list f_alloc_int()
{
  int* aaa= alloc_int();
  int n= 10;
  list x;
  for(int i=0; i<n; i++) {
    x.append(aaa[i]);
  }
  return x;
}

void f_operate_list(list& alist)
{
  int* intlist= new int [10];
  for (int i=0; i<10; i++) {
    intlist[i]= extract<int>(alist[i]);
  }

  operate_list(intlist);

  for (int i=0; i<10; i++) {
    alist[i]= intlist[i];
  }

  delete intlist;
}

BOOST_PYTHON_MODULE(test13)
{
  def("alloc_int",    f_alloc_int);
  def("operate_list", f_operate_list);
}

