// Plotting with Asymptote (http://asymptote.sourceforge.net)

import graph;

size(200, 150, IgnoreAspect);


file in=line(input("plot1.dat"));
real[][] a=transpose(dimension(in, 0, 0));

real[] x = a[0];
real[] y = a[1];
real[] z = a[2];
real[] d = a[3];

//    std::setw(15) << averageMultiplicity / eventNumber << 
//    std::setw(15) << averageProtonNumber / eventNumber <<
//    std::setw(15) << averageNeutronNumber / eventNumber << " " <<
//    std::setw(15) << averageNucleonKinEnergy / (averageProtonNumber + averageNeutronNumber) << " " <<
//    std::setw(15) << averageProtonKinEnergy / (averageProtonNumber + 1.0e-10) << " " <<
//    std::setw(15) << averageNeutronKinEnergy / (averageNeutronNumber + 1.0e-10) << " " <<
//    std::setw(15) << averagePionNumber / eventNumber << " " <<
//    std::setw(15) << averagePionKinEnergy / (averagePionNumber + 1.0e-10) << G4endl;

draw(graph(x, y), red);
draw(graph(x, z), blue);
draw(graph(x, d), black);

xaxis("$E$", BottomTop, LeftTicks);
yaxis("$avg. num. of particles$", LeftRight, RightTicks);

