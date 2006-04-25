// $Id: test12.cc,v 1.2 2006-04-25 07:54:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test12.cc
//
//   test for indexing a STL vector
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#if BOOST_VERSION < 103300
#include "vector_indexing_suite.hpp"
#else 
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include <vector>
#include <iostream>

class AClass {
private:
  int ival;

public:
  AClass() : ival(0) {
  }

  AClass(const AClass& right) {
    ival= right.ival;
    std::cout << "*** copy constructor is called" << std::endl;
  }

  ~AClass() { }

  void SetIVal(int i) { ival= i; }
  int GetIVal() { return ival; }

  void Print() const {
    std::cout << "*** @" << this << ": i=" 
	      << ival << std::endl;
  }

};

typedef std::vector<AClass*> AVector;


void PrintVector(const AVector& vec) {
  std::cout << "*** size of AVector= " << vec.size() << std::endl;
  for (int i=0; i< vec.size(); i++) {
    vec[i]-> Print();
  }
}


// Boost.Python...

using namespace boost::python;

BOOST_PYTHON_MODULE(test12)
{
  class_<AClass, AClass*>("AClass", "a class")
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .def("Print",         &AClass::Print)
    ;

  class_<AVector> ("AVector", "AClass vector")
    .def(vector_indexing_suite<AVector>())
    ;

  def("PrintVector", PrintVector);
}

