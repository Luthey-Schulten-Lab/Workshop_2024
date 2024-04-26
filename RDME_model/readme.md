# Instruction Guide

# Example of usage


```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
from jLM.RDME import Sim as RDMESim                   # Main simulation class
from jLM.RegionBuilder import RegionBuilder           # Helper class to design the site lattice
from jLM.Solvers import ConstBoundaryConc, makeSolver # Run the simulation
from lm import MpdRdmeSolver                          # lm::rdme::MpdRdmeSolver
import jLM                                            # Set up the Jupyter environment
```

Our example system will be a $128\times 128\times 128$ lattice with a $32\ nm$ lattice spacing. We are 
going to name it "Example", and it will write the Lattice Microbes data file to "tmp.lm". The simulation 
time step is $25\ \mu s$, and the default region is named "external".


```python
lmFile =  "tmp.lm"

sim = RDMESim("Example",
              lmFile,
              [128,128,128],
              32e-9,
              "external",
              dt=25e-6)
```

We define a second region, "free", where the chemical reactions will occur. Here `sim.region()` either looks 
up an existing region type, or creates a new type if the name isn't already taken. We get handles to the `free` 
and `external` regions. These variables will be used as nouns in our computational description of the system.


```python
free = sim.region("free")
external = sim.region("external")
```

In order to design the site lattice geometry, we create a `RegionBuilder` object, which reads the lattice dimensions 
from our simulation object `sim`.


```python
B = RegionBuilder(sim)
```

To describe the geometry of a compartment in a lattice based representation, it's easiest to think in terms of set 
operations. Our simulation geometry here is the space between two concentric spheres. We define 3-D boolean arrays 
of the same dimensions as the site lattice (called binary masks) such that `True` elements indicate that the 
corresponding subvolume is included in that compartment. The `external` region is defined by  as
$$\mathtt{external} = \mathtt{bigSphere}^\complement \cup \mathtt{littleSphere},  $$
and naturally
$$\mathtt{free} = \mathtt{external}^\complement.$$
Finally, to create the site lattice in our `sim` object, we call `B.compose`. This function takes a series of tuples 
where the first element is the region object that we are defining the geometry for, and the second element is the 
binary mask. To check our work, we call `sim.showRegionStack()` to get an interactive visualization of the site 
lattice.


```python
outerSph = B.ellipsoid(radius=60,center=(64,64,64))
innerSph = B.ellipsoid(radius=10,center=(64,64,64))
extMask = ~outerSph | innerSph
freeMask = ~extMask
B.compose((free, freeMask), (external, extMask))
sim.showRegionStack()
```

Now we add some simple chemical reactions. Species, rate constant, and reaction types are created using the same 
semantics as for regions. When rate constants are defined, their reaction order must be specified after their 
reaction rate. Reaction rates are all given in terms of deterministic rate constants.

Chemical reactions within a region are specified by two lists, reactants and products, composed of species objects, 
and the reaction rate object. An empty list denotes nothing on that side of the reaction.

We place a molecule at a specific location using `Species.placeParticle`, and distribute a number of particles in a 
region uniformly using `Species.placeNumberInto`.

After constructing the reactions, we can check our work using `sim.showReactions`


```python
outerSpecies = sim.species("Outer")
innerSpecies = sim.species("Inner")
productSpecies = sim.species("Product")
dmrRate = sim.rateConst("dmr", 1e6, 2)
free.addReaction([innerSpecies, outerSpecies], [productSpecies], dmrRate)

sourceSpecies = sim.species("Source")
sourceSpecies.placeParticle(64,64,32,1)

degSpecies = sim.species("Deg")
degSpecies.placeNumberInto(free, 10)

birthRate = sim.rateConst("degBirth", 100., 1)
free.addReaction([sourceSpecies], [degSpecies, sourceSpecies], birthRate)

dcyRate = sim.rateConst("degDeath", 1.0, 1)
free.addReaction([degSpecies], [], dcyRate)

anhRate = sim.rateConst("annhil", 1e7, 2)
free.addReaction([degSpecies, productSpecies], [], anhRate)

sim.showReactions()
```

To define the diffusive properties of the chemical species, we use two basic techniques. To set multiple entries in 
the diffusion matrix, we use `sim.transitionRate(species, fromRegion, toRegion, diffusionRate`. If the species or 
region terms are `None`, that acts as a wild card which sets all entries along that axis of the diffusion matrix. 
This allows an easy way to set up the default diffusion parameters first, then specialize the diffusion rates for
specific species.

To define the diffusion rate for specific species, the method `Species.diffusionRate` is called.

To check our work so far, we use `sim.showAllSpecies`. This provides a dropdown menu that selects the species type. 
If the diffusion matrix is underspecified, the missing elements will be highlighted here.


```python
sim.transitionRate(None, external, external, sim.diffusionFast)
sim.transitionRate(None, free, external, sim.diffusionZero)
sim.transitionRate(None, external, free, sim.diffusionFast)

diff = sim.diffusionConst("normal", 1e-13)
diffSlow = sim.diffusionConst("slow", 1e-15)

outerSpecies.diffusionRate(free, diff)
innerSpecies.diffusionRate(free, diff)
productSpecies.diffusionRate(free, diffSlow)
sourceSpecies.diffusionRate(free, sim.diffusionZero)
degSpecies.diffusionRate(free, diffSlow)

sim.showAllSpecies()
```

To be sure that we placed the $\mathrm{Source}$ molecule in the right place with respect to the lattice geometry, we 
can use the 3-D lattice viewer, `sim.displayGeometry`. The viewer has a key which shows which colors correspond to 
which region. Included in the key is a button which toggles the visibility of the region or species. After making the 
`external` region transparent and hiding the `free` region, we can see that the $\mathrm{Source}$ species is placed
reasonably.


```python
sim.displayGeometry()
```

Now that we've defined the species, reactions, geometry, and diffusion constants, we can finalize the simulation and 
write a LM data file. First we need to specify the rate that the full particle lattice is written to disk. Here we 
set `sim.latticeWriteInterval` to 0.01 seconds, meaning that every 10 ms of simulated time, the state of the particle 
lattice will be written. Similarly, we set `sim.speciesWriteInterval` which governs the rate in which only the total 
particle counts are written. We set the total simulated time to 2 seconds, and write the simulation file to disk.


```python
sim.latticeWriteInterval = 0.01
sim.speciesWriteInterval = sim.latticeWriteInterval
sim.simulationTime = 2
sim.finalize()
sim
```

At this point, we could run the simulation, however we did not specify any concentration for the $\mathrm{Inner}$ 
and $\mathrm{Outer}$ species. We are going to introduce these species into the simulation through constant 
concentration boundary conditions. The $\mathrm{Outer}$ species will be created at the outer surface, and the 
$\mathrm{Inner}$ species at the inner surface. 

To describe the boundary, we use a binary mask which indicate whether or not a subvolume is included in the boundary. 
To create a shell that overlaps with the `free` region and touches the `external` region, we use binary morphological 
operations on our previously created sphere masks. The outer boundary is simply the outer sphere minus the 
complement of the dilated outer sphere. Here we use a structuring element which includes only the six nearest
neighbors of each subvolume. Similarly for the inner boundary, it is simply the dilated inner sphere minus the 
original inner sphere.

To make sure we described this correctly, we will use `RegionBuilder.showBinaryLattices` to view our construction. To
allow us to see into the spheres, we will only display half of them. In order to only show the half spheres, the 
`filterFunction` option is used. This takes a dictionary where the keys are the names of regions specified in the 
function call, and the values are functions such that $f(x,y,z)\to\mathrm{Bool}$ is a function that returns true when 
the site index $(x,y,z)$ should be visible.


```python
outerBoundary = outerSph & ~ B.erode(outerSph, se=B.se6)
innerBoundary = B.dilate(innerSph, se=B.se6) & ~ innerSph

B.showBinaryLattices(dict(external=extMask, outerBoundary=outerBoundary, innerBoundary=innerBoundary),
                    filterFunctions=dict(external=lambda x,y,z: z>64,
                                         outerBoundary=lambda x,y,z: z>64,
                                         innerBoundary=lambda x,y,z:z > 64))
```

To actually use the boundary conditions, we must use a derived solver class. `ConstBoundaryConc` is a solver mixin 
defined in `jLM` which overrides `hookSimulation` in order to reset the concentration of species at the boundaries. 
`hookSimulation` is only called when the particle lattice is written to disk, so the lattice should be written 
sufficiently frequently to ensure that the boundary concentrations are indeed constant.

To create the solver, we use the single GPU MPD-RDME solver, `lm.MpdRdmeSolver`, however other solvers could be used 
as well (e.g. `lm.MGPUMpdRdmeSolver`.) The convenience function `makeSolver` composes the custom solver and the LM 
solver.

After instantiating the solver, we specify the boundary conditions in terms of the list of fixed concentration 
species, the list of requested concentrations, and the binary mask describing the boundary.


```python
Solver = makeSolver(MpdRdmeSolver, ConstBoundaryConc)
solver = Solver(sim)
solver.setBoundary([innerSpecies], [1e-5], innerBoundary)
solver.setBoundary([outerSpecies], [1e-5], outerBoundary)
```

Finally, we run the simulation in the notebook. This returns a new simulation object with the results of the 
simulation.


```python
traj = sim.run(solver=solver, cudaDevices=[0])
```

From the new simulation object, we extract the total particle numbers as a function of time, the spatial 
distribution of particles as a function of time, and the spatial distribution of particles. It should be noted here 
that the `Species`, `Region`, etc. data types from before are no longer valid here since they were associated with 
the `sim` object, not the new `traj` object.


```python
fig,((ax0,ax1),(ax2,ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(12,10))

def plotNum(sp):
    ts, ns = sp.getNumberTrajectory()
    ax0.semilogy(ts, ns, label=sp.name)
    
plotNum(traj.sp.Outer)
plotNum(traj.sp.Inner)
plotNum(traj.sp.Product)
plotNum(traj.sp.Source)
plotNum(traj.sp.Deg)
ax0.legend()
ax0.set(title="Absolute number", xlabel="t/s", ylabel="Particle count")
    
ts, degN =  traj.sp.Deg.getLatticeTrajectory(integrated='xy')
ax1.imshow(degN.T,aspect='auto',extent=(0,ts[-1],0,degN.shape[1]))
ax1.set(title="Deg", xlabel="t/s", ylabel="z-index")

ts, outerN =  traj.sp.Outer.getLatticeTrajectory(integrated='xy')
ax2.imshow(outerN.T,aspect='auto',extent=(0,ts[-1],0,outerN.shape[1]))
ax2.set(title="Outer", xlabel="t/s", ylabel="z-index")

prodHist = traj.sp.Product.getLatticeHistogram(integrated='z',timeStart=1.8, timeEnd=2.0)
ax3.imshow(prodHist,aspect='auto')
ax3.set(title="Product", xlabel="x-index", ylabel="y-index")

fig.tight_layout()
```
