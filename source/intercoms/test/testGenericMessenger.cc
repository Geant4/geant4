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
//
// $Id: testG4VoxelLimits.cc,v 1.4 2006-06-29 18:59:02 gunter Exp $

// testGenericMessanger
//
//

#include "G4UImanager.hh"
#include "G4GenericMessenger.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

#undef NDEBUG
#include <assert.h>

class testGenericMessenger {
public:
  testGenericMessenger();
  ~testGenericMessenger() {}
public:
  bool b;
  int i;
  unsigned int ui;
  long l;
  unsigned long ul;
  float f;
  double d;
  double d_u;
  G4complex c;
  G4String s;
  G4String s2;
  G4ThreeVector v3;
  G4ThreeVector v3_u;
  bool m0_called;
  bool m1_called;
  int  m1_arg;
  bool m2_called;
  G4String m2_arg;
  bool m3_called;
  std::string m3_arg0;
  int m3_arg1;
  void method0() { m0_called = true; }
  void method1(double a) {m1_called = true; m1_arg = a;}
  void method2(const G4String& s) {m2_called = true; m2_arg = s;}
  int  method3(const std::string& s, int n) { m3_called = true; m3_arg0 = s; m3_arg1 = n; return n;}
private:
  G4GenericMessenger messenger;
};

testGenericMessenger::testGenericMessenger() :
b(false), i(0), ui(0), l(0), ul(0), f(0.0), d(0.0), d_u(0.0),
m0_called(false), m1_called(false), m1_arg(0),
m2_called(false), m3_called(false), m3_arg1(0), messenger(this, "/mytest/") {
  messenger.DeclareProperty("bool", b, "bool property");
  messenger.DeclareProperty("int", i, "int property");
  messenger.DeclareProperty("uint", ui, "unsigned int property").SetRange("value>=0");
  messenger.DeclareProperty("long", l, "long property");
  messenger.DeclareProperty("ulong", ul, "unsigned long property");
  messenger.DeclareProperty("float", f, "float property");
  messenger.DeclareProperty("double", d, "double property");
  messenger.DeclareProperty("complex", c, "complex property");
  messenger.DeclareProperty("string", s, "string property");
  messenger.DeclareProperty("string2", s2, "string property with other options")
           .SetParameterName("String", true)
           .SetCandidates("string1 string2 string3")
           .SetDefaultValue("string2");
  messenger.DeclareProperty("3vector", v3, "G4ThreeVector property");
  
  messenger.DeclareProperty("doubleunit", d_u)
           .SetGuidance("double energy property")
           .SetGuidance("additional guidance")
           .SetDefaultUnit("MeV")
           .SetStates(G4State_PreInit, G4State_Idle);
  messenger.DeclareProperty("3vectorunit", v3_u, "3 vector length property")
           .SetUnitCategory("Length");
  
  messenger.DeclareMethod("method0", &testGenericMessenger::method0, "Method with 0 arguments");
  messenger.DeclareMethod("method1", &testGenericMessenger::method1, "Method with 1 arguments");
  messenger.DeclareMethod("method1unit", &testGenericMessenger::method1, "Method with 1 arguments").SetUnit("mm");
  messenger.DeclareMethod("method2", &testGenericMessenger::method2, "Method with 1 const ref argument");
  messenger.DeclareMethod("method3", &testGenericMessenger::method3, "Method with 2 arguments");
}


int main()
{
  //--- create a our test object-------------------------------------
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //---bool----------------------------------------------------------
  assert(tst.b == false);
  
  UImanager->ApplyCommand("/mytest/bool 1");
  assert(tst.b == true);
  UImanager->ApplyCommand("/mytest/bool 0");
  assert(tst.b == false);
  
  UImanager->ApplyCommand("/mytest/bool F");
  assert(tst.b == false);
  UImanager->ApplyCommand("/mytest/bool T");
  assert(tst.b == true);
  
  UImanager->ApplyCommand("/mytest/bool false");
  assert(tst.b == false);
  UImanager->ApplyCommand("/mytest/bool true");
  assert(tst.b == true);
  
  assert(UImanager->GetCurrentIntValue("/mytest/bool") == 1);
  
  //---int-----------------------------------------------------------
  assert(tst.i == 0);
  UImanager->ApplyCommand("/mytest/int 99");
  assert(tst.i == 99);
  UImanager->ApplyCommand("/mytest/int -99");
  assert(tst.i == -99);
  
  assert(UImanager->GetCurrentIntValue("/mytest/int") == -99);
  
  //---unsigned int--------------------------------------------------
  assert(tst.ui == 0);
  UImanager->ApplyCommand("/mytest/uint 88");
  assert(tst.ui == 88);
  assert(UImanager->GetCurrentIntValue("/mytest/uint") == 88);
  
  //---long----------------------------------------------------------
  assert(tst.l == 0);
  UImanager->ApplyCommand("/mytest/long 999");
  assert(tst.l == 999L);
  
  //---unsigned long-------------------------------------------------
  assert(tst.ul == 0);
  UImanager->ApplyCommand("/mytest/ulong 888");
  assert(tst.ul == 888L);

  //---float---------------------------------------------------------
  assert(tst.f == 0.0);
  UImanager->ApplyCommand("/mytest/float 9.9");
  assert(tst.f == 9.9F);
  UImanager->ApplyCommand("/mytest/float -9.9E+10");
  assert(tst.f == -9.9E+10F);

  //---double---------------------------------------------------------
  assert(tst.d == 0.0);
  UImanager->ApplyCommand("/mytest/double 99.99");
  assert(tst.d == 99.99);
  UImanager->ApplyCommand("/mytest/double -99.99E+10");
  assert(tst.d == -99.99E+10);
  
  assert(UImanager->GetCurrentDoubleValue("/mytest/double") == -99.99E+10);
  
  assert(tst.d_u == 0.0);
  UImanager->ApplyCommand("/mytest/doubleunit 77.8 TeV");
  assert(tst.d_u == 77800000.0);

  assert(UImanager->GetCurrentDoubleValue("/mytest/doubleunit") ==  77800000.0);
  
  //---std::complex---------------------------------------------------
  assert(tst.c == G4complex(0.0,0.0));
  UImanager->ApplyCommand("/mytest/complex (8.8,9.9)");
  assert(tst.c == G4complex(8.8,9.9));
  
  //---std::string----------------------------------------------------
  assert(tst.s == G4String());
  UImanager->ApplyCommand("/mytest/string mystring");
  assert(tst.s == "mystring");
  UImanager->ApplyCommand("/mytest/string \"another string\"");
  assert(tst.s == "another string");
  
  //---std::string with other options---------------------------------
  assert(tst.s2 == G4String());
  UImanager->ApplyCommand("/mytest/string2");
  assert(tst.s2 == "string2");
  UImanager->ApplyCommand("/mytest/string2 string3");
  assert(tst.s2 == "string3");

  //---G4ThreeVector--------------------------------------------------
  assert(tst.v3 == G4ThreeVector());
  UImanager->ApplyCommand("/mytest/3vector 1. 2. 3.");
  assert(tst.v3 == G4ThreeVector(1.0, 2.0, 3.0));
  
  assert(tst.v3_u == G4ThreeVector());
  UImanager->ApplyCommand("/mytest/3vectorunit 4.1 5.1 6.1 cm");
  assert(tst.v3_u == G4ThreeVector(41., 51., 61));
  
  //---Method with 0 arguments----------------------------------------
  assert(tst.m0_called == false);
  UImanager->ApplyCommand("/mytest/method0");
  assert(tst.m0_called == true);

  //---Method with 1 arguments----------------------------------------
  assert(tst.m1_called == false);
  UImanager->ApplyCommand("/mytest/method1 999.");
  assert(tst.m1_called == true);
  assert(tst.m1_arg == 999.);
  
  UImanager->ApplyCommand("/mytest/method1unit 9.0 km");
  assert(tst.m1_arg == 9.0E+6);
  
  //---Method with 1 const ref arguments----------------------------------------
  assert(tst.m2_called == false);
  UImanager->ApplyCommand("/mytest/method2 abcdef");
  assert(tst.m2_called == true);
  assert(tst.m2_arg == "abcdef");

  //---Method with 2 arguments----------------------------------------
  assert(tst.m3_called == false);
  UImanager->ApplyCommand("/mytest/method3 abcdef 999");
  assert(tst.m3_called == true);
  assert(tst.m3_arg0 == "abcdef");
  assert(tst.m3_arg1 == 999);
  
  G4cout << "All done" << G4endl;
  return 0;
}

