//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: test13.cc,v 1.3 2006-06-04 21:36:00 kmura Exp $
// $Name: not supported by cvs2svn $
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

