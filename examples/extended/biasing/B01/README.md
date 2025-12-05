\page ExampleB01 Example B01 

## Importance sampling

The example uses importance sampling or the weight window technique 
according to an input parameter. It uses scoring in both cases. 
Importance values or weight windows are defined according to the mass 
geometry. In this example the weight window technique is configured such 
that it behaves equivalent to importance sampling: The window is actually 
not a window but simply the inverse of the importance value and only
one energy region is used that covers all energies in the problem.
The user may change the weight window configuration by changing the
initialization of the weight window algorithm in example,cc. 
Different energy bounds for the weight window technique may be specified 
in B01DetectorConstruction.

The executable takes one optional argument: 0 or 1. Without argument or
with argument: 0, the importance sampling is applied with argument: 1,
the weight window technique is applied.

See also [Category "biasing"](../../html/Examples_biasing.html) documentation.

