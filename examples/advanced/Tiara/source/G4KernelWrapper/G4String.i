# $Id: G4String.i,v 1.1 2003/06/17 15:25:01 mdressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-01 $
# -------------------------------------------------------------------

// coppied from std_string.i

%include exception.i

%{
#include "G4String.hh"
%}


class G4String;

/* Overloading check */

%typemap(typecheck) G4String = char *;
%typemap(typecheck) const G4String & = char *;

%typemap(in) G4String {
  if (PyString_Check($input))
    $1 = G4String(PyString_AsString($input));
  else
    SWIG_exception(SWIG_TypeError, "G4String expected");
}

%typemap(in) const G4String & (G4String temp) {
  if (PyString_Check($input)) {
    temp = G4String(PyString_AsString($input));
    $1 = &temp;
  } else {
    SWIG_exception(SWIG_TypeError, "G4String expected");
  }
}

%typemap(out) G4String {
  $result = PyString_FromStringAndSize($1.data(),$1.size());
}

%typemap(out) const G4String & {
  $result = PyString_FromStringAndSize($1->data(),$1->size());
}

%typemap(inv, parse="s") G4String, const G4String &, G4String & "$1_name.c_str()";

%typemap(inv, parse="s") G4String *, const G4String * "$1_name->c_str()";

%typemap(outv) G4String {
  if (PyString_Check($input))
    $result = G4String(PyString_AsString($input));
  else
    throw SWIG_DIRECTOR_TYPE_MISMATCH("G4String expected");
}

%typemap(outv) const G4String & (G4String temp) {
  if (PyString_Check($input)) {
    temp = G4String(PyString_AsString($input));
    $result = &temp;
  } else {
    throw SWIG_DIRECTOR_TYPE_MISMATCH("G4String expected");
  }
}



