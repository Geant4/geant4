#include "TranslateFragment.h"

#include <iostream>

string TranslateFragment(const string& vo)
{
  string translation;
  if      (vo == "neutron")  translation = "n";
  else if (vo == "gamma")    translation = "#gamma";
  else if (vo == "proton")   translation = "p";
  else if (vo == "alpha")    translation = "#alpha";
  else if (vo == "triton")   translation = "t";
  else if (vo == "He3")      translation = "^{3}He";
  else if (vo == "deuteron") translation = "d";
  else if (vo == "pi-")      translation = "#pi^{-}";
  else if (vo == "pi+")      translation = "#pi^{+}";
  else if (vo == "pi0")      translation = "#pi^{0}";
  else if (vo == "pi")       translation = "#pi";
  else if (vo == "eta")      translation = "#eta";
  else 
    {
      string::size_type bra_pos = vo.find("[");
      if (bra_pos != string::npos)
	{
	  string tmp = vo.substr(0,bra_pos);
	  string::size_type idx = tmp.find_first_of("0123456789",0);
	  translation = "^{" + tmp.substr(idx,string::npos) + "}" + tmp.substr(0,idx);
	}
      else 
	{
	  translation="noname";
	  std::cerr << "TranslateFragment: Particle not found " << vo << '\n';
	}
    }
  return translation;
}
