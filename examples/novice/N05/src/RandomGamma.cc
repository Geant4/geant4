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
// --------------------------------------------
// Adapted from RNGAMA function of the KERNLIB
// Shoot at random according to the gamma
// distribution.
// --------------------------------------------
#include "RandomGamma.hh"
#include <CLHEP/Random/Randomize.h>

double RandomGamma(double p)
{
  double h;
  double stor[15];
  
  h = 0.;
  if (p > 15)
    {
      while (h <= 0.)
	{
	  double a;
	  a = RandExponential::shoot();
	  h = p*(1.-1./(9*p)+a/(3*sqrt(p)));
	  h = h*h*h;
	}
    }
  else
    {
      int m = p;
      double f = p - m;
      if (m > 0)
	{
	  double x = 1.;
	  RandFlat::shootArray(m, stor);
	  for (int i = 0; i < m; i++) x *= stor[i];
	  h = -log(x);
	}
      if (f >= 0.00001)
	{
	  double x = RandFlat::shoot();
	  double x1 = -log(x);
	  if (f >= 0.9999) h += x1;
	  else
	    {
	      x = RandFlat::shoot();
	      double wlog;
	      while (1)
		{
		  wlog = log(x)/f;
		  if (wlog > -100.)
		    {
		      double w1 = exp(wlog);
		      x = RandFlat::shoot();
		      wlog = log(x)/(1.-f);
		      if (wlog < -100.)
			{
			  h += x1;
			  break;
			}
		      else
			{
			  double w = w1 + exp(wlog);
			  if (w <= 1.)
			    {
			      h += w1/w;
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
  return h;
}
