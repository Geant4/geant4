global tree, hf

tree = tf.create("gamma_2000.his",1,1,"hbook")
tree.thisown=1
hf = af.createHistogramFactory(tree)

# ... and load them into memory:
hEsource       = tree.findH1D("10")
hDeposit       = tree.findH1D("20")
hNuclearDep    = tree.findH1D("30")
hNumberPhLow   = tree.findH1D("40")
hNumberPhHigh  = tree.findH1D("50")
hPhArrivalTime = tree.findH1D("60")

# set plotter to 3*2 zones
pl.createRegions(3,2)

# ... and plot the histograms
pl.plot(hEsource      ) ; pl.show() ; pl.next()
pl.plot(hDeposit      ) ; pl.show() ; pl.next()
pl.plot(hNuclearDep   ) ; pl.show() ; pl.next()
pl.plot(hNumberPhLow  ) ; pl.show() ; pl.next()
pl.plot(hNumberPhHigh ) ; pl.show() ; pl.next()
pl.plot(hPhArrivalTime) ; pl.show() ; pl.next()

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
nt1 = tree.findTuple("3" )

# plot a few quantities using the shortcut
nplot(nt1,"tot_e"  ,"",0.,20.)
nplot(nt1,"xe_time","",0.,100.)

# prompt user for <return>
wait()

# get the secondaries ntuple and plot a few attributes from it
nt2 = tree.findTuple("2" )
nplot(nt2, "xpos", "", 0., 50)
nplot(nt2, "zpos", "", 0., 100.)

pl.write("pmtpos.ps", "ps")

tree.commit()
tree.close()
del tree




