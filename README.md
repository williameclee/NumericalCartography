# Numerical Cartography

This repository contains the code for the numerical cartography project. It contains some basic functionalities to create lighting effects on a 3D terrain model. The code is written in MATLAB and with a CPU-based approach.

## Demos

Try to run the demo of the function `shadow.m`:

```matlab
shadow('demo')
```

It will generate a simple 3D terrain model with a light source and a shadow.

## Functions

### Shading

- `grainynoise`: Adds a grainy noise to a normal map.
- `hillshade`: Computes the hill shade of a 3D terrain model with various Lambert models.
- `normal`: Computes the normal vector of a 3D terrain model.
- `sampledem` and `sampledem2`: Generate sample 3D terrain models.
- `shadow`: Computse the shadow of a 3D terrain model with a light source.

### Projection

- `equalearth` & `equalearthd`: Equal Earth projection.
- `gridearth` & `gridearthd`: Gridded Earth projection, which is an equal-area projection that has its bounding box fit a given aspect ratio (and more specifications)!

## Objects

- `LightSource`: A class to define a light source.

## Contributors

Note that some of the functions are migrated from the [slepian_ulmo](https://github.com/williameclee/slepian_ulmo) project.

- [En-Chi Lee](mailto:williameclee@gmail.com)

Last modified by:
[En-Chi Lee](mailto:williameclee@gmail.com) 2024-10-18
