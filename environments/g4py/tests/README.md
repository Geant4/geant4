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
|test00  | hallo world
|test01  | simple class wrapping
|test02  | simple inheritance
|test03  | singleton w/o public constructor <br> can NOT be compiled.
|test04  | call policies
|test05  | function overloading & default arguments
|test06  | python inheritance from base class
|test07  | enums
|test08  | static member function
|test09  | operators
|test10  | call by-reference
|test11  | no_init
|test12  | indexing a STL vector
|test13  | test for a list argument/return


Function tests
--------------

| Test   | Description |
|:------:|-------------
|gtest01 | user application
|gtest02 | test for using site-module packages
|gtest03 | test for EZsim package
|gtest04 | test for getting command tree and command information
|gtest05 | test for constructing CSG geometries in Python
|gtest06 | test for constructing/visualizing boolean geoemtries
|gtest07 | test for checking overlapped geometries


