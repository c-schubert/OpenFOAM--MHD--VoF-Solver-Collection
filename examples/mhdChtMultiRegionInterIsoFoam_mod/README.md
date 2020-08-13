# Example of mhdChtMultiRegionInterIsoFoam_mod solver

## Steps

- Generate Elmer Mesh and copy to `./meshElmer` (mesh.boundary, mesh.elements, mesh.header, mesh.nodes) 
  - Mesh can be generates from `./meshElmer/mesh.grd` file (using ElmerGrind oder ElmerGUI)
- Run `Allrun.pre`
- Adapt `./system/controlDict`
- Adapt number of cores in `system/**` decomposeParDicts
- Adapt and run `Allrun.eof`