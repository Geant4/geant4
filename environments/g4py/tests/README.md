Test Modules
============

This module collects a set of tests for python wrapping of C++ elements.

>**Notes:** <br>
"*testXX*" directories contain tests without Geant4 package. <br>
"*gtestXX*" directories contain tests with Geant4 package.


Basic tests
-----------

| Test   | Description |
|:------:|-------------
|test00  | hallo world (K.M)
|test01  | simple class wrapping (K.M)
|test02  | simple inheritance (K.M)
|test03  | singleton w/o public constructor (K.M) <br> can NOT be compiled.
|test04  | call policies (K.M)
|test05  | function overloading & default arguments (K.M)
|test06  | python inheritance from base class (K.M)
|test07  | enums (K.M)
|test08  | static member function (K.M)
|test09  | operators (K.M)
|test10  | call by-reference (K.M)
|test11  | no_init (K.M)
|test12  | indexing a STL vector (K.M)                                
|test13  | test for a list argument/return (K.M)


Function tests
--------------

| Test   | Description |
|:------:|-------------
|gtest01 | user application (K.M)
|gtest02 | test for using site-module packages (K.M)
|gtest03 | test for EZsim package (K.M)
|gtest04 | test for getting command tree and command information (K.M)
|gtest05 | test for constructing CSG geometries in Python (K.M)
|gtest06 | test for constructing/visualizing boolean geoemtries (K.M)
|gtest07 | test for checking overlapped geometries (K.M)

