#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ==================================================================
# EM Calculator (GTK)
#
# Koichi Murakami (KEK/CRC)
# ==================================================================
from Geant4 import *
import g4py.EMSTDpl
import g4py.emcalculator
import EmPlot
import ROOT
import sys
from cStringIO import StringIO

# ==================================================================
# Geant4
# ==================================================================
def g4_configure() :
  EmPlot.Configure()
  g4py.EMSTDpl.Construct()

# -------------------------------------------------------------------
# plot for chaged particles
# -------------------------------------------------------------------
def plot_charged(material, pname) :
  EmPlot.SetMaterial(material)

  # initialize G4 kernel
  gRunManager.Initialize()
  gRunManagerKernel.RunInitialization()

  # energy
  elist= []
  for n in range(-3, 3):
    for i in range(10,99):
      elist.append(i/10.*10.**n *MeV)

  # calculate stopping power
  global mycout
  mycout.close()
  sys.stdout = mycout = StringIO()
  dedx_list= g4py.emcalculator.CalculateDEDX(pname, material, elist, 1)
  xlist_tot=[]
  xlist_ioni=[]
  xlist_brems=[]

  for x in dedx_list:
    xlist_tot.append((x[0], x[1]["tot"]/(MeV*cm2/g)))
    xlist_ioni.append((x[0], x[1]["ioni"]/(MeV*cm2/g)))
    xlist_brems.append((x[0], x[1]["brems"]/(MeV*cm2/g)))

  # make plot
  global myCanvas, aplot, bplot, cplot
  myCanvas = EmPlot.init_root()
  aplot = EmPlot.make_plot(xlist_tot, pname+" Stopping Power ("+material+")",
                           "dE/dX (MeV cm^{2}/g)")
  bplot = EmPlot.make_plot(xlist_ioni, "Stopping Power ("+material+")",
                           "dE/dX (MeV cm^{2}/g)", 1)
  cplot = EmPlot.make_plot(xlist_brems, "Stopping Power ("+material+")",
                           "dE/dX (MeV cm^{2}/g)", 3)
  myCanvas.SaveAs("/tmp/sp.png")


# -------------------------------------------------------------------
# plot for gamma
# -------------------------------------------------------------------
def plot_gamma(material) :
  EmPlot.SetMaterial(material)

  # initialize G4 kernel
  gRunManager.Initialize()
  gRunManagerKernel.RunInitialization()

  # energy
  elist= []
  for n in range(-3, 4):
    for i in range(10,99):
      elist.append(i/10.*10.**n *MeV)

  # calculate cross sections
  global mycout
  mycout.close()
  sys.stdout = mycout = StringIO()
  xsection_list= g4py.emcalculator.CalculatePhotonCrossSection(material,
                                                               elist, 1)
  xlist_tot=[]
  xlist_comp=[]
  xlist_pe=[]
  xlist_conv=[]
  for x in xsection_list:
    xlist_tot.append((x[0]/MeV, x[1]["tot"]/(cm2/g)))
    xlist_comp.append((x[0]/MeV, x[1]["compt"]/(cm2/g)))
    xlist_pe.append((x[0]/MeV, x[1]["phot"]/(cm2/g)))
    xlist_conv.append((x[0]/MeV, x[1]["conv"]/(cm2/g)))

  # make plots
  global myCanvas, aplot, bplot, cplot, dplot
  myCanvas = EmPlot.init_root()
  aplot = EmPlot.make_plot(xlist_tot,  "Photon Cross Section ("+material+")",
                             "Cross Section (cm^{2}/g)")
  bplo = EmPlot.make_plot(xlist_comp, "Photon Cross Section ("+material+")",
                             "Cross Section (cm^{2}/g)", 1)
  cplot = EmPlot.make_plot(xlist_pe,   "Photon Cross Section ("+material+")",
                             "Cross Section (cm^{2}/g)", 7)
  dplot = EmPlot.make_plot(xlist_conv, "Photon Cross Section ("+material+")",
                             "Cross Section (cm^{2}/g)", 3)
  myCanvas.SaveAs("/tmp/cs.png")


# ==================================================================
# GUI
# ==================================================================
import pygtk
pygtk.require20()
import gtk

# -------------------------------------------------------------------
# main window
# -------------------------------------------------------------------
class MainWindow :
  def __init__(self) :
    self.__margin = 8

    # main window
    self.mainwindow = gtk.Window(gtk.WINDOW_TOPLEVEL)
    self.mainwindow.set_title('Em Calculator')
    self.mainwindow.set_position(gtk.WIN_POS_MOUSE)
    self.mainwindow.set_default_size(500, 300)
    self.mainwindow.connect('delete-event', lambda w,d: gtk.main_quit())
    self.mainwindow.connect('destroy-event',lambda w,d: gtk.main_quit())

    # components
    header = self.create_header()

    particle_frame = self.create_particle_frame()
    particle_frame.set_border_width(self.__margin)

    material_frame = self.create_material_frame()
    material_frame.set_border_width(self.__margin)

    separator = gtk.HSeparator()

    action_box = self.create_action_box()
    action_box.set_border_width(self.__margin)

    # layout
    vbox = gtk.VBox()
    vbox.pack_start(header)
    vbox.pack_start(particle_frame)
    vbox.pack_start(material_frame)
    vbox.pack_start(separator)
    vbox.pack_start(action_box)

    self.mainwindow.add(vbox)
    self.mainwindow.show_all()

    # text view
    self.textview = TextView()

    # error dialog
    self.error_dialog = gtk.MessageDialog(parent=self.mainwindow,
      buttons=gtk.BUTTONS_CLOSE, type=gtk.MESSAGE_ERROR,
      message_format="Material is not defined in G4Nist materials")
    self.error_dialog.connect("response", self.cb_close_dialog)

  def create_header(self) :
    hbox = gtk.HBox()

    label = gtk.Label()
    label.set_markup("<big><b>EM Calculator</b></big>")
    hbox.pack_start(label)

    label = gtk.Label()
    text = """

Shows
 - stopping power for e/mu/proton
 - cross sections for gamma
"""
    label.set_markup(text)
    hbox.pack_start(label)

    return hbox

  def create_particle_frame(self) :
    frame = gtk.Frame("Particle")

    hbox = gtk.HBox()
    frame.add(hbox)
    hbox.set_border_width(self.__margin)

    button = gtk.RadioButton(None, "electron")
    button.connect("toggled", self.cb_select_particle, "e-")
    hbox.pack_start(button, True, True, 0)
    button.show()
    self.particle = "e-"

    button = gtk.RadioButton(button, "positron")
    button.connect("toggled", self.cb_select_particle, "e+")
    hbox.pack_start(button, True, True, 0)
    button.show()

    button = gtk.RadioButton(button, "mu-")
    button.connect("toggled", self.cb_select_particle, "mu-")
    hbox.pack_start(button, True, True, 0)
    button.show()

    button = gtk.RadioButton(button, "mu+")
    button.connect("toggled", self.cb_select_particle, "mu+")
    hbox.pack_start(button, True, True, 0)
    button.show()

    button = gtk.RadioButton(button, "proton")
    button.connect("toggled", self.cb_select_particle, "proton")
    hbox.pack_start(button, True, True, 0)
    button.show()

    button = gtk.RadioButton(button, "gamma")
    button.connect("toggled", self.cb_select_particle, "gamma")
    hbox.pack_start(button, True, True, 0)
    button.show()

    return frame

  def create_material_frame(self) :
    frame = gtk.Frame("Material (G4Nist)")

    hbox = gtk.HBox()
    frame.add(hbox)
    hbox.set_border_width(self.__margin)

    self.material_list = [ "G4_Al", "G4_Si", "G4_Ar", "G4_Cu", "G4_Fe",
                           "G4_Ge", "G4_Ag", "G4_W", "G4_Au", "G4_Pb",
                           "G4_AIR", "G4_Galactic", "G4_WATER", "G4_CESIUM_IODIDE",
                           "G4_SODIUM_IODIDE", "G4_PLASTIC_SC_VINYLTOLUENE",
                          "G4_MYLAR" ]

    self.material_combo = gtk.combo_box_entry_new_text()
    hbox.pack_start(self.material_combo)
    for name in self.material_list :
      self.material_combo.append_text(name)
    self.material_combo.set_active(0)
    self.material = self.material_list[0]
    self.material_combo.connect("changed", self.cb_select_material)

    return frame

  def create_action_box(self) :
    box = gtk.HButtonBox()
    box.set_layout(gtk.BUTTONBOX_END)
    box.set_spacing(self.__margin)

    exec_button = gtk.Button(stock = gtk.STOCK_EXECUTE)
    text_button = gtk.Button("Text View")
    quit_button = gtk.Button(stock = gtk.STOCK_QUIT)
    box.add(exec_button)
    box.add(text_button)
    box.add(quit_button)

    exec_button.connect("clicked", self.cb_show_plot)
    text_button.connect("clicked", self.cb_show_textview)
    quit_button.connect("clicked", lambda w: gtk.main_quit())

    return box

  # callbacks
  def cb_show_textview(self, widget, data=None) :
    window = self.textview.get_window()
    window.show_all()

  def cb_select_particle(self, widget, data=None) :
    self.particle = data

  def cb_select_material(self, widget, data=None) :
    entry = widget.get_child()
    self.material = entry.get_text()

  def cb_show_plot(self, widget, data=None) :
    g4mate = gNistManager.FindOrBuildMaterial(self.material)
    if (g4mate == None) :
      self.error_dialog.show_all()
      return

    if (self.material_list.count(self.material) == 0) :
      self.material_combo.append_text(self.material)

    if (self.particle == "gamma" ) :
      plot_gamma(self.material)
    else :
      plot_charged(self.material, self.particle)
    self.textview.textbuffer.set_text(mycout.getvalue())

  def cb_close_dialog(self, widget, data=None) :
    widget.hide_all()


# -------------------------------------------------------------------
# text view
# -------------------------------------------------------------------
class TextView :
  def __init__(self) :
    self.__margin = 8
    self.text_window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    self.text_window.set_title('Value with Text')
    self.text_window.set_position(gtk.WIN_POS_MOUSE)
    self.text_window.set_default_size(500, 300)

    vbox = gtk.VBox()
    self.text_window.add(vbox)

    sw = gtk.ScrolledWindow()
    sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
    sw.set_border_width(self.__margin)
    vbox.pack_start(sw)

    textview = gtk.TextView()
    self.textbuffer = textview.get_buffer()
    sw.add(textview)

    hbox = gtk.HButtonBox()
    hbox.set_layout(gtk.BUTTONBOX_END)
    hbox.set_border_width(self.__margin)
    vbox.pack_start(hbox, expand=False)

    close_button = gtk.Button(stock = gtk.STOCK_CLOSE)
    close_button.connect("clicked", self.cb_hide_window)
    hbox.add(close_button)

  def get_window(self) :
    return self.text_window

  def cb_hide_window(self, widget, data=None) :
    self.text_window.hide_all()
    return False


# ==================================================================
# main
# ==================================================================
def main() :
  SetG4PyCoutDestination()

  default_stdout = sys.stdout
  global mycout
  sys.stdout = mycout = StringIO()

  # G4 setup
  g4_configure()

  # start GUI
  application = MainWindow()
  gtk.main()

  gTerminate()
  sys.stdout = default_stdout


if __name__ == "__main__":
  try :
    main()
  except KeyboardInterrupt :
    pass
