global tree, hf

tree = tf.create("comptonhisto.hbook",1,1,"hbook")
tree.thisown=1
hf = af.createHistogramFactory(tree)

# ... and load them into memory:
hEKin    = tree.findH1D("10")
hP       = tree.findH1D("20")
hNSec    = tree.findH1D("30")
hDeposit = tree.findH1D("40")
hTheta   = tree.findH1D("50")
hPhi     = tree.findH1D("60")

# set plotter to 3*2 zones
pl.createRegions(3,2)

# ... and plot the histograms
pl.plot(hEKin   ) ; pl.show() ; pl.next()
pl.plot(hP      ) ; pl.show() ; pl.next()
pl.plot(hNSec   ) ; pl.show() ; pl.next()
pl.plot(hDeposit) ; pl.show() ; pl.next()
pl.plot(hTheta  ) ; pl.show() ; pl.next()
pl.plot(hPhi    ) ; pl.show() ; pl.next()

wait()

# reset number of zones
pl.createRegions(2,1)

# helper function
def nplot(tup, col, cut, xmin, xmax, nbins=100):
  global hf
  h1 = hf.create1D("1000000","temp hist", nbins, xmin, xmax)
  colId = tup.findColumn(col)
  bNextRow=tup.start()                                         # Looping over the tuple entries
  while bNextRow:
    h1.fill(tup.getFloat(colId))
    bNextRow=tup.next()                                        # Retrieving the subsequent row
  pl.plot(h1) ; pl.show() ; pl.next()
  tree.rm("1000000")
  return

# get the primary ntuple from the NtupleManager
nt1 = tree.findTuple("1" )

# plot a few quantities using the shortcut
nplot(nt1,"initen" ,"",0,20)
nplot(nt1,"dedx"   ,"",0.,100)

# prompt user for <return>
wait()

# get the secondaries ntuple and plot a few attributes from it
nt2 = tree.findTuple("2" )
nplot(nt2, "parttyp", "", 0., 42)
nplot(nt2, "e"      , "", 0., 100.)

pl.write("secondaries.ps", "ps")

tree.commit()
tree.close()
del tree




