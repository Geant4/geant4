// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
//-------------------------------------------------------------

#include <cctype>
#include <iostream>

namespace {

bool eatwhitespace ( std::istream & is ) {
  // Will discard whitespace until it either encounters EOF or bad input
  // (in which case it will return false) or it hits a non-whitespace.  
  // Will put that non whitespace character back so that after this routine
  // returns true, is.get(c) should always work.
  // If eatwhitespace returns false, is will always be in a fail or bad state.
  char c;
  bool avail = false;	// avail stays false until we know there is a nonwhite
			// character available.
  while ( is.get(c) ) {
    if ( !isspace(c) ) {
      is.putback(c);
      avail = true;
      break;
    }
  }
  return avail;
}

void fouledup() {
  std::cerr << "istream mysteriously lost a putback character!\n";
}


} // end of unnamed namespace


namespace CLHEP  {

void ZMinput3doubles ( std::istream & is, const char * type,
			double & x, double & y, double & z ) {

// Accepted formats are 
// x y z
// x, y, z (each comma is optional, and whitespace ignored if comma present)
// ( x, y, z ) (commas optional)

  char c;
  bool parenthesis = false;

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before trying to input " << type << "\n";
    return;
  }

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == '(' ) {
    parenthesis = true;
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended after ( trying to input " << type << "\n";
      return;
    }
  } else {
    is.putback(c);
  }  

  // At this point, parenthesis or not, the next item read is supposed to
  // be the number x.

  if (!(is >> x)) {
    std::cerr << "Could not read first value in input of " << type << "\n";
    return;
  }

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before second value of " << type << "\n";
    return;
  } 

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == ',' ) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended ater one value and comma in " 
							<< type << "\n";
      return;
    }
  } else {
    is.putback(c);
  }

  // At this point, comma or not, the next item read is supposed to
  // be the number y.

  if (!(is >> y)) {
    std::cerr << "Could not read second value in input of " << type << "\n";
    return;
  }

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before third value of " << type << "\n";
    return;
  } 

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == ',' ) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended ater two values and comma in " 
							<< type << "\n";
      return;
    }
  } else {
    is.putback(c);
  }

  // At this point, comma or not, the next item read is supposed to
  // be the number z.

  if (!(is >> z)) {
    std::cerr << "Could not read third value in input of " << type << "\n";
    return;
  }

  // Finally, check for the closing parenthesis if there was an open paren.

  if (parenthesis) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "No closing parenthesis in input of " << type << "\n";
      return;
    } 
    if ( !is.get(c) ) { fouledup(); return; }
    if ( c != ')' ) {
      std::cerr << "Missing closing parenthesis in input of " 
							<< type << "\n";
      // Now a trick to do (as nearly as we can) what 
      // is.putback(c); is.setstate(std::ios_base::failbit); 
      // would do (because using ios_base will confuse old CLHEP compilers):
      if ( isdigit(c) || (c=='-') || (c=='+') ) {
        is.putback('@');
      } else {
        is.putback('c');
      }
      int m;
      is >> m;  // This fails, leaving the state bad, and the istream
		// otherwise unchanged, except if the next char might
		// have started a valid int, it turns to @
      return;
    }
  }

  return;

}

void ZMinputAxisAngle ( std::istream & is, 
			double & x, double & y, double & z, 
			double & delta ) {
// Accepted formats are 
// parenthesis optional, then
// any acceptable format for a Hep3Vector, then
// optional comma, then
// delta, then
// close parenthesis if opened at start.
//
// But if there is an open parenthesis, it must be for the overall
// object.  That is, if the axis has parentheses, the form must be 
// ( (x,y,z) , delta ) 

  char c;
  bool parenthesis = false;

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before trying to input AxisAngle \n";
    return;
  }

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == '(' ) {
    parenthesis = true;
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended after ( trying to input AxisAngle \n";
      return;
    }
  } else {
    is.putback(c);
  }  

  // At this point, parenthesis or not, the next item read is supposed to
  // be a valid Hep3Vector axis.

  ZMinput3doubles ( is, "axis of AxisAngle", x, y, z );
  if (!is) return;

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before delta of AxisAngle \n";
    return;
  } 

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == ',' ) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended ater axis and comma in AxisAngle \n"; 
      return;
    }
  } else {
    is.putback(c);
  }

  // At this point, comma or not, the next item read is supposed to
  // be the number delta.

  if (!(is >> delta)) {
    std::cerr << "Could not delta value in input of AxisAngle \n";
    return;
  }

  // Finally, check for the closing parenthesis if there was an open paren.

  if (parenthesis) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "No closing parenthesis in input of AxisAngle \n";
      return;
    } 
    if ( !is.get(c) ) { fouledup(); return; }
    if ( c != ')' ) {
      std::cerr << "Missing closing parenthesis in input of AxisAngle \n";
      if ( isdigit(c) || (c=='-') || (c=='+') ) {
        is.putback('@');
      } else {
        is.putback('c');
      }
      int m;
      is >> m;  // This fails, leaving the state bad.
      return;
    }
  }

  return;

}

void ZMinput2doubles ( std::istream & is, const char * type,
			double & x, double & y ) {

// Accepted formats are 
// x y 
// x, y (comma is optional, and whitespace ignored if comma present)
// ( x, y ) (comma optional)

  char c;
  bool parenthesis = false;

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before trying to input " << type << "\n";
    return;
  }

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == '(' ) {
    parenthesis = true;
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended after ( trying to input " << type << "\n";
      return;
    }
  } else {
    is.putback(c);
  }  

  // At this point, parenthesis or not, the next item read is supposed to
  // be the number x.

  if (!(is >> x)) {
    std::cerr << "Could not read first value in input of " << type << "\n";
    return;
  }

  if ( !eatwhitespace(is) ) {
    std::cerr << "istream ended before second value of " << type << "\n";
    return;
  } 

  if ( !is.get(c) ) { fouledup(); return; }
  if ( c == ',' ) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "istream ended ater one value and comma in " 
							<< type << "\n";
      return;
    }
  } else {
    is.putback(c);
  }

  // At this point, comma or not, the next item read is supposed to
  // be the number y.

  if (!(is >> y)) {
    std::cerr << "Could not read second value in input of " << type << "\n";
    return;
  }

  // Finally, check for the closing parenthesis if there was an open paren.

  if (parenthesis) {
    if ( !eatwhitespace(is) ) {
      std::cerr << "No closing parenthesis in input of " << type << "\n";
      return;
    } 
    if ( !is.get(c) ) { fouledup(); return; }
    if ( c != ')' ) {
      std::cerr << "Missing closing parenthesis in input of " 
							<< type << "\n";
      // Now a trick to do (as nearly as we can) what 
      // is.putback(c); is.setstate(std::ios_base::failbit); 
      // would do (because using ios_base will confuse old CLHEP compilers):
      if ( isdigit(c) || (c=='-') || (c=='+') ) {
        is.putback('@');
      } else {
        is.putback('c');
      }
      int m;
      is >> m;  // This fails, leaving the state bad, and the istream
		// otherwise unchanged, except if the next char might
		// have started a valid int, it turns to @
      return;
    }
  }

  return;

}

}  // namespace CLHEP
