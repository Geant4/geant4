#!/usr/bin/python
# ==================================================================
# An example of ploting by EmCalculator
#
# Plotting photon cross sections and stopping power with ROOT
# ==================================================================
from Geant4 import *
import g4py.NISTmaterials
import g4py.ezgeom

# ==================================================================
# geometry setup
# ==================================================================

# ------------------------------------------------------------------
# setup
# ------------------------------------------------------------------
def Configure():
  g4py.NISTmaterials.Construct()
  g4py.ezgeom.Construct()

# ------------------------------------------------------------------
# constructing geometry
# ------------------------------------------------------------------
def SetMaterial(material_name):
  material= gNistManager.FindOrBuildMaterial(material_name)
  g4py.ezgeom.SetWorldMaterial(material)


# ==================================================================
# plot by ROOT
# ==================================================================
import ROOT
from math import log, log10, sqrt, ceil, floor
from array import array

# ------------------------------------------------------------------
#   caclculate plot range
# ------------------------------------------------------------------
def plot_range(xmin, xmax, xmargin=0.):
  xmaxlog= 10
  xminlog= -10
  
  if(xmax!=0.):
    xmaxlog= log10(xmax)

  if(xmin!=0):
    xminlog= log10(xmin)

  ixmaxlog= xmaxlog+0.5
  ixminlog= xminlog-0.5-xmargin

  return [10**ixminlog, 10**ixmaxlog]

# ------------------------------------------------------------------
#   ROOT init
# ------------------------------------------------------------------
def init_root():
  ROOT.gROOT.Reset()

  # plot style
  ROOT.gStyle.SetTextFont(82)

  ROOT.gStyle.SetTitleFont(82, "X")
  ROOT.gStyle.SetTitleFontSize(0.04)
  ROOT.gStyle.SetLabelFont(82, "X")
  ROOT.gStyle.SetTitleFont(82, "Y")
  ROOT.gStyle.SetLabelFont(82, "Y")

  #ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetErrorX(0)

  canvas= ROOT.TCanvas("g4plot", "g4plot", 620, 30, 600, 600)

  canvas.SetLogy()
  canvas.SetLogx()
  canvas.SetGrid()

  return canvas

# ------------------------------------------------------------------ 
#   do a plot
# ------------------------------------------------------------------ 
def make_plot(xlist, user_title, axis_titile, q_super_impose=0):

  ekin_array, y_array = array('d'), array('d')

  for x in xlist:
    ekin_array.append(x[0])
    y_array.append(x[1])
    
  # plot range
  xmin= min(ekin_array)
  xmax= max(ekin_array)
  xrange= plot_range(xmin, xmax)
  
  ymin= min(y_array)
  ymax= max(y_array)
  yrange= plot_range(ymin, ymax, 2)

  if(q_super_impose==0):
    htit= user_title
    global frame
    frame= ROOT.TH1F("dumy", htit, 1, xrange[0], xrange[1]);
    frame.SetMinimum(yrange[0]);
    frame.SetMaximum(yrange[1]);
    frame.SetXTitle("Kinetic Energy (MeV)")
    frame.GetXaxis().SetLabelSize(0.025)
    frame.GetXaxis().SetTitleSize(0.03)
    frame.SetYTitle(axis_titile)
    frame.GetYaxis().SetLabelSize(0.025)
    frame.GetYaxis().SetTitleSize(0.03)
    frame.SetStats(0)
    frame.Draw()

  plot= ROOT.TGraph(len(ekin_array), ekin_array, y_array)
  plot.Draw("L")
  plot.SetLineColor(q_super_impose+1)

  return plot

