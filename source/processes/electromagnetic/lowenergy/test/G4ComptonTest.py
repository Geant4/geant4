
# tell HistoManger where to find the histograms:
hm.selectStore("comptonhisto.hbook")
# ... and load them into memory:
hEKin    = hm.load1D(10)
hP       = hm.load1D(20)
hNSec    = hm.load1D(30)
hDeposit = hm.load1D(40)
hTheta   = hm.load1D(50)
hPhi     = hm.load1D(60)

# set plotter to 3*2 zones
pl.zone(3,2)
# ... and plot the histograms
hplot(10)
hplot(20)
hplot(30)
hplot(40)
hplot(50)
hplot(60)

# --------------------------------------------------------------------------------
# define a helper function for simple ntuple plotting
def nplot(_nt,_var,_cut="", _min, _max, _first=0, _num=1000000000) :
    _histNtPlot = hm.create1D(1000000,_var,100, _min, _max)
    _nt.project1D(_histNtPlot, _var, _cut, _first, _num)
    vTemp=vm.from1D(_histNtPlot)
    pl.plot(vTemp)
    del vTemp
    return 
# --------------------------------------------------------------------------------

# set plotter to 2*1 zone (x*y)
pl.zone(2,1)

# get the primary ntuple from the NtupleManager
nt1 = ntm.findNtuple("comptonhisto1.hbook::1" )

# see which attributes there are:
nt1.listAttributes()

# plot a few quantities using the shortcut
nplot(nt1,"energyf","",0,20)
nplot(nt1,"dedx","",0.,100)

# get the secondaries ntuple and plot a few attributes from it
nt2 = ntm.findNtuple("comptonhisto2.hbook::2" )
nplot(nt2, "type", "", 0., 42)
nplot(nt2, "e", "", 0., 100.)





