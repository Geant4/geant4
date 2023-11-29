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
// G4StatDouble class implementation
//
// Original Author: Giovanni Santin (ESA) - October 2005 in GRAS tool
// Adapted by: John Apostolakis - November 2011
// --------------------------------------------------------------------
#include "G4StatDouble.hh"

G4StatDouble::G4StatDouble() { reset(); }

G4StatDouble::G4StatDouble(G4double x) { fill(x); }

void G4StatDouble::reset()
{
  m_sum_wx  = 0.;
  m_sum_wx2 = 0.;
  m_n       = 0;
  m_sum_w   = 0.;
  m_sum_w2  = 0.;
  m_scale   = 1.;
}

void G4StatDouble::fill(G4double value, G4double weight)
{
  m_sum_wx += value * weight;
  m_sum_wx2 += value * value * weight;
  if(m_n < INT_MAX)
  {
    ++m_n;
  }
  m_sum_w += weight;
  m_sum_w2 += weight * weight;

  if(weight <= 0.)
  {
    G4cout << "[G4StatDouble::fill] WARNING: weight<=0. " << weight << G4endl;
  }
}

void G4StatDouble::scale(G4double value) { m_scale = m_scale * value; }

G4double G4StatDouble::mean() const
{
  G4double mean_val = 0.;
  if(m_sum_w > 0.)
  {
    mean_val = m_sum_wx / m_sum_w;
  }
  return m_scale * mean_val;
}

G4double G4StatDouble::mean(G4double ext_sum_w) const
{
  G4double factor = 0.;
  // factor to rescale the Mean for the requested number
  // of events (or sum of weights) ext_sum_w

  if(ext_sum_w > 0)
  {
    factor = m_sum_w;
    factor /= ext_sum_w;
  }
  return mean() * factor;
}

G4double G4StatDouble::rms(G4double ssum_wx, G4double ssum_wx2, G4double ssum_w,
                           G4int nn)
{
  G4double vrms = 0.0;
  if(nn > 1)
  {
    G4double vmean = ssum_wx / ssum_w;
    G4double xn    = nn;
    G4double tmp =
      // from GNU Scientific Library. This part is equivalent to N/(N-1)
      // when w_i = w
      // ((m_sum_w * m_sum_w) / (m_sum_w * m_sum_w - m_sum_w2))

      // from NIST "DATAPLOT Reference manual", Page 2-66
      // http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
      // rewritten based on: SUM[w(x-m)^2]/SUM[w] = SUM[wx^2]/SUM[w] - m^2
      // and dividing it by sqrt[n] to go from rms of distribution to the
      // rms of the mean value

      (xn / (xn - 1)) * ((ssum_wx2 / ssum_w) - (vmean * vmean));

    tmp  = std::max(tmp, 0.0);  // this avoids observed computation problem
    vrms = std::sqrt(tmp);
    //  G4cout << "[G4StatDoubleElement::rms] m_sum_wx: " << m_sum_wx
    //         << "  m_sum_wx2: " << m_sum_wx2 << "  m_sum_w: " << m_sum_w
    //         << "  m_n: " << m_n << "  tmp: " << tmp<< "  rms: " << rms
    //         << G4endl;
    //  G4cout << "[G4StatDoubleElement::rms] (m_n / (m_n - 1)): " << (xn/(xn -
    //  1))
    // 	   << "  (m_sum_wx2 / m_sum_w): " << (m_sum_wx2 / m_sum_w)
    // 	   << "  (mean * mean): " << (mean * mean)
    // 	   << "  ((m_sum_wx2 / m_sum_w) - (mean * mean)): "
    //         << ((m_sum_wx2 / m_sum_w) - (mean * mean))
    // 	   << G4endl;
  }
  return vrms * m_scale;
}

G4double G4StatDouble::rms()
{
  // this method computes the RMS with "all internal" parameters:
  // all the sums are the internal ones: m_sum_wx, m_sum_wx2, m_sum_w, m_n

  return rms(m_sum_wx, m_sum_wx2, m_sum_w, m_n);
}

G4double G4StatDouble::rms(G4double ext_sum_w, G4int ext_n)
{
  // this method computes the RMS with sum_w and n coming from outside:
  // ext_sum_w and ext_n:
  // this means that the result is normalised to the external events
  // it is useful when, given a number ext_n of events with sum of the weights
  // ext_sum_w, only m_n (with sum of weights m_sum_w) are actually accumulated
  // in the internal summation (e.g. for a dose variable in a volume, because
  // only a few particles reach that volume)

  return rms(m_sum_wx, m_sum_wx2, ext_sum_w, ext_n);
}

void G4StatDouble::add(const G4StatDouble* ptr)
{
  m_n += ptr->n();
  m_sum_w += ptr->sum_w();
  m_sum_w2 += ptr->sum_w2();
  m_sum_wx += ptr->sum_wx();
  m_sum_wx2 += ptr->sum_wx2();
}
