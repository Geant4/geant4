// $Id: test13.cc,v 1.1 2006-02-27 10:05:25 kmura Exp $
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

