global tree, hm

tree = tf.create("gammaraytel0.hbook",1,1,"hbook")
tree.thisown=1
hm = af.createHistogramFactory(tree)

# ... and load them into memory:
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

# get the primary ntuple from the NtupleManager
nt1 = tree.findTuple("1")

# see which attributes there are:
nt1.listAttributes()









