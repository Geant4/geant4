# tell HistoManager where to find the histograms:
# change the name to access hbook files different from the first run
hm.selectStore("gammaraytel0.hbook")

# ... and load them into memory:
hE    = hm.load1D(1)
hPl   = hm.load1D(2)
hXZ   = hm.load2D(3)
hYZ   = hm.load2D(4)

# set plotter to 3*2 zones
pl.zone(2,2)
# ... and plot the histograms
hplot(hE)
hplot(hPl)
hplot(hXZ)
hplot(hYZ)

# get the primary ntuple from the NtupleManager
nt1 = ntm.findNtuple("gammaraytel0.hbook::1" )

# see which attributes there are:
nt1.listAttributes()









