global tree, hm

# Here use the name of the hbook you want to analyze
tree = tf.create("gammaraytel.hbook",1,1,"hbook")
tree.thisown=1
hm = af.createHistogramFactory(tree)

# Retrieve histograms and load them into memory:
hE    = tree.findH1D("10")
hPl   = tree.findH1D("20")
hXZ   = tree.findH2D("30")
hYZ   = tree.findH2D("40")

# set plotter to 2*2 zones
pl.createRegions(2,2)

# ... and plot the histograms
pl.plot(hE      ) ; pl.show() ; pl.next()
pl.plot(hPl      ) ; pl.show() ; pl.next()
pl.plot(hXZ   ) ; pl.show() ; pl.next()
pl.plot(hYZ  ) ; pl.show() ; pl.next()

# get the primary ntuple from the NtupleManager (just a check)
nt1 = tree.findTuple("1")










