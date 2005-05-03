from AidaProxy import *
treeSim = Proxy_Store("xrayfluo.xml","xml",3)
histoSim = treeSim.retrieveH1D("1")
from rootPlotter2 import *  
pl = RootPlotter()          
pl.plot(histoSim)
