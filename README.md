# Simplex 2D / 3D in Rust

Provides basic 2D and 3D simplex noise functions.

* https://en.wikipedia.org/wiki/Simplex_noise

This Rust version is based on this public domain Java implementation:

```
Simplex noise demystified
Stefan Gustavson, Linköping University, Sweden (stegu@itn.liu.se), 2005-03-22
```

* https://github.com/stegu/perlin-noise/blob/master/simplexnoise.pdf

* https://web.archive.org/web/20210506195658/http://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java

## Usage

```rust
use simplex_23d::SimplexNoise;
```

## Visualizations

![simplex noise 2d](noise2d.png)

![simplex noise 3d](noise3d.gif)