
# Creating a memory-mapped tree

tree=tf.create()                                                 

# Creating a histogram factory mapped to the tree

hf=af.createHistogramFactory( tree )                             

# Open an existing HBook file

treeHBook=tf.create("gammaraytel.hbook", "hbook", 1, 0 )              

# Mounting the hbook tree under the master memory tree.

tree.mkdir( "hbook" )
tree.mount( "/hbook", treeHBook, "/" )

# Retrieve histograms and load them into memory:

# Fetching the histograms from hbook

hE=tree.findH1D( "/hbook/10" )
hPl=tree.findH1D( "/hbook/20" )
hXZ=tree.findH2D("/hbook/30")
hYZ=tree.findH2D("/hbook/40")

# set plotter to 2*2 zones
#pl.createRegions(2,2)
#>>> NOTE! ONLY SINGLE REGION IN THIS VERSION!

# plot the histograms

pr = pl.currentRegion()
pr.plot(hE,"")
pl.refresh()

wait()

pr.clear()
pr.plot(hPl,"")
pl.refresh()

wait()

del hf
del tree
del treeHBook
exit()








