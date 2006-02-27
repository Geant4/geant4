// $Id: test09.cc,v 1.1 2006-02-27 10:05:25 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test09.cc
//
//   test for operators
//
//                                         2005 Q
// ====================================================================
#include <iostream>

class AClass {
private:
  int ival;

public:
  AClass() :ival(999) { }
  AClass(int i) :ival(i) { }
  ~AClass() { }

  void SetIVal(int i) { ival= i; }
  int GetIVal() const { return ival; }

  AClass operator+(const AClass& aclass) {
    AClass atemp;
    atemp.ival= ival + aclass.ival;    
    return atemp;
  }

  AClass& operator+=(const AClass& aclass) {
    ival+= aclass.ival;
    return *this;
  }

  bool operator==(const AClass& aclass) const {
    if(ival == aclass.ival) return true;
    return false;
  }
};

std::ostream& operator<<(std::ostream& ostr, const AClass& aclass)
{
  return ostr << aclass.GetIVal();
}


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test09)
{
  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def(init<int>())
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .def(self + self)
    .def(self += self)
    .def(self == self)
    .def(self_ns::str(self))
    ;
}

