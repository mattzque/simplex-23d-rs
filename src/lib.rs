#[cfg(test)]
mod tests;

// Rust implementation of:
//
// Simplex noise demystified
// Stefan Gustavson, Linköping University, Sweden (stegu@itn.liu.se), 2005-03-22
//
// https://github.com/stegu/perlin-noise/blob/master/simplexnoise.pdf
// 
// https://web.archive.org/web/20210506195658/http://weber.itn.liu.se/~stegu/simplexnoise/SimplexNoise.java
//

const SQRT3: f32 = 1.732_050_8; //  f32::sqrt(3.0);
const F2: f32 = 0.5 * (SQRT3 - 1.0);
const G2: f32 = (3.0 - SQRT3) / 6.0;
const F3: f32 = 1.0 / 3.0;
const G3: f32 = 1.0 / 6.0; // Very nice and simple unskew factor, too

const GRAD3: [[f32; 3]; 12] = [
    [1.0, 1.0, 0.0],
    [-1.0, 1.0, 0.0],
    [1.0, -1.0, 0.0],
    [-1.0, -1.0, 0.0],
    [1.0, 0.0, 1.0],
    [-1.0, 0.0, 1.0],
    [1.0, 0.0, -1.0],
    [-1.0, 0.0, -1.0],
    [0.0, 1.0, 1.0],
    [0.0, -1.0, 1.0],
    [0.0, 1.0, -1.0],
    [0.0, -1.0, -1.0],
];

fn generate_perm_mod12(seed: u64) -> [usize; 512] {
    use rand::{seq::SliceRandom, SeedableRng};
    // random number generator with seed:
    let mut rng: rand::rngs::StdRng = SeedableRng::seed_from_u64(seed);
    let mut seq: Vec<_> = (0usize..256).collect();
    seq.shuffle(&mut rng);
    let mut perm = [0; 512];
    // To remove the need for index wrapping, double the permutation table length
    // mod 12 to lookup the gradients of the simplex corners
    for i in 0..512 {
        perm[i] = seq[i & 255] % 12;
    }
    perm
}

#[inline(always)]
fn fastfloor(x: f32) -> i32 {
    if x > 0.0 {
        x as i32
    } else {
        (x as i32) - 1
    }
}

#[inline(always)]
fn dot2(x0: f32, y0: f32, x1: f32, y1: f32) -> f32 {
    x0 * x1 + y0 * y1
}

#[inline(always)]
fn dot3(x0: f32, y0: f32, z0: f32, x1: f32, y1: f32, z1: f32) -> f32 {
    x0 * x1 + y0 * y1 + z0 * z1
}

/// Simplex noise generator instance that keeps a permutation table internally.
pub struct Simplex {
    seed: u64,
    perm_mod12: [usize; 512],
}

impl Simplex {
    /// Create a new simplex noise generator for the given seed value.
    /// 
    /// Creates a pre-computed permutation table using the seed value.
    pub fn new(seed: u64) -> Self {
        Self {
            seed,
            perm_mod12: generate_perm_mod12(seed),
        }
    }

    /// Sample the noise function in 2D at the given coordinates.
    /// 
    /// For frequency, multiply the coordinates by the desired frequency.
    pub fn sample2d(&self, x: f32, y: f32) -> f32 {
        // Skew the input space to determine which simplex cell we're in

        let s = (x + y) * F2; // Hairy factor for 2D
        let i = fastfloor(x + s);
        let j = fastfloor(y + s);
        let t = ((i + j) as f32) * G2;
        let x0 = (i as f32) - t; // Unskew the cell origin back to (x,y) space
        let y0 = (j as f32) - t;
        let x0 = x - x0; // The x,y distances from the cell origin
        let y0 = y - y0;

        // For the 2D case, the simplex shape is an equilateral triangle.
        // Determine which simplex we are in.
        // Offsets for second (middle) corner of simplex in (i,j) coords
        let (i1, j1) = if x0 > y0 {
            (1, 0) // lower triangle, XY order: (0,0)->(1,0)->(1,1)
        } else {
            (0, 1) // upper triangle, YX order: (0,0)->(0,1)->(1,1)
        };

        // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
        // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
        // c = (3-sqrt(3))/6
        let x1 = x0 - (i1 as f32) + G2; // Offsets for middle corner in (x,y) unskewed coords
        let y1 = y0 - (j1 as f32) + G2;
        let x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
        let y2 = y0 - 1.0 + 2.0 * G2;

        // Work out the hashed gradient indices of the three simplex corners
        let ii: usize = (i & 255) as usize;
        let jj: usize = (j & 255) as usize;
        let gi0 = self.perm_mod12[ii + self.perm_mod12[jj]];
        let gi1 = self.perm_mod12[ii + i1 + self.perm_mod12[jj + j1]];
        let gi2 = self.perm_mod12[ii + 1 + self.perm_mod12[jj + 1]];

        // Calculate the contribution from the three corners
        // double n0, n1, n2; // Noise contributions from the three corners
        let mut t0 = 0.5 - x0 * x0 - y0 * y0;
        let n0 = if t0 < 0.0 {
            0.0
        } else {
            t0 *= t0;
            t0 * t0 * dot2(GRAD3[gi0][0], GRAD3[gi0][1], x0, y0) // (x,y) of grad3 used for 2D gradient
        };

        let mut t1 = 0.5 - x1 * x1 - y1 * y1;
        let n1 = if t1 < 0.0 {
            0.0
        } else {
            t1 *= t1;
            t1 * t1 * dot2(GRAD3[gi1][0], GRAD3[gi1][1], x1, y1)
        };

        let mut t2 = 0.5 - x2 * x2 - y2 * y2;
        let n2 = if t2 < 0.0 {
            0.0
        } else {
            t2 *= t2;
            t2 * t2 * dot2(GRAD3[gi2][0], GRAD3[gi2][1], x2, y2)
        };

        // Add contributions from each corner to get the final noise value.
        // The result is scaled to return values in the interval [-1,1].
        n0 + n1 + n2
    }

    /// Sample the noise function in 3D at the given coordinates.
    /// 
    /// For frequency, multiply the coordinates by the desired frequency.
    pub fn sample3d(&self, x: f32, y: f32, z: f32) -> f32 {
        // Skew the input space to determine which simplex cell we're in
        let s = (x + y + z) * F3; // Very nice and simple skew factor for 3D
        let i = fastfloor(x + s);
        let j = fastfloor(y + s);
        let k = fastfloor(z + s);
        let t = (i + j + k) as f32 * G3;
        let x0 = i as f32 - t; // Unskew the cell origin back to (x,y,z) space
        let y0 = j as f32 - t;
        let z0 = k as f32 - t;
        let x0 = x - x0; // The x,y,z distances from the cell origin
        let y0 = y - y0;
        let z0 = z - z0;
        // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
        // Determine which simplex we are in.

        // Offsets for second corner of simplex in (i,j,k) coords
        let i1;
        let j1;
        let k1;

        // Offsets for third corner of simplex in (i,j,k) coords
        let i2;
        let j2;
        let k2;

        if x0 >= y0 {
            if y0 >= z0 {
                i1 = 1;
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            }
            // X Y Z order
            else if x0 >= z0 {
                i1 = 1;
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            }
            // X Z Y order
            else {
                i1 = 0;
                j1 = 0;
                k1 = 1;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            } // Z X Y order
        } else {
            // x0<y0
            if y0 < z0 {
                i1 = 0;
                j1 = 0;
                k1 = 1;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            // Z Y X order
            else if x0 < z0 {
                i1 = 0;
                j1 = 1;
                k1 = 0;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            // Y Z X order
            else {
                i1 = 0;
                j1 = 1;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            } // Y X Z order
        }
        // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
        // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
        // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
        // c = 1/6.
        let x1 = x0 - i1 as f32 + G3; // Offsets for second corner in (x,y,z) coords
        let y1 = y0 - j1 as f32 + G3;
        let z1 = z0 - k1 as f32 + G3;
        let x2 = x0 - i2 as f32 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
        let y2 = y0 - j2 as f32 + 2.0 * G3;
        let z2 = z0 - k2 as f32 + 2.0 * G3;
        let x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
        let y3 = y0 - 1.0 + 3.0 * G3;
        let z3 = z0 - 1.0 + 3.0 * G3;
        // Work out the hashed gradient indices of the four simplex corners
        let ii = (i & 255) as usize;
        let jj = (j & 255) as usize;
        let kk = (k & 255) as usize;
        let gi0 = self.perm_mod12[ii + self.perm_mod12[jj + self.perm_mod12[kk]]];
        let gi1 = self.perm_mod12[ii + i1 + self.perm_mod12[jj + j1 + self.perm_mod12[kk + k1]]];
        let gi2 = self.perm_mod12[ii + i2 + self.perm_mod12[jj + j2 + self.perm_mod12[kk + k2]]];
        let gi3 = self.perm_mod12[ii + 1 + self.perm_mod12[jj + 1 + self.perm_mod12[kk + 1]]];
        // Calculate the contribution from the four corners
        // let n0, n1, n2, n3;
        // Noise contributions from the four corners
        let mut t0 = 0.5 - x0 * x0 - y0 * y0 - z0 * z0;
        let n0 = if t0 < 0.0 {
            0.0
        } else {
            t0 *= t0;
            t0 * t0 * dot3(GRAD3[gi0][0], GRAD3[gi0][1], GRAD3[gi0][2], x0, y0, z0)
        };

        let mut t1 = 0.5 - x1 * x1 - y1 * y1 - z1 * z1;
        let n1 = if t1 < 0.0 {
            0.0
        } else {
            t1 *= t1;
            t1 * t1 * dot3(GRAD3[gi1][0], GRAD3[gi1][1], GRAD3[gi1][2], x1, y1, z1)
        };

        let mut t2 = 0.5 - x2 * x2 - y2 * y2 - z2 * z2;
        let n2 = if t2 < 0.0 {
            0.0
        } else {
            t2 *= t2;
            t2 * t2 * dot3(GRAD3[gi2][0], GRAD3[gi2][1], GRAD3[gi2][2], x2, y2, z2)
        };

        let mut t3 = 0.5 - x3 * x3 - y3 * y3 - z3 * z3;
        let n3 = if t3 < 0.0 {
            0.0
        } else {
            t3 *= t3;
            t3 * t3 * dot3(GRAD3[gi3][0], GRAD3[gi3][1], GRAD3[gi3][2], x3, y3, z3)
        };

        // Add contributions from each corner to get the final noise value.
        // The result is scaled to stay just inside [-1,1]
        n0 + n1 + n2 + n3
    }

    /// Returns the seed value that was used to generate the permutation table.
    pub fn seed(&self) -> u64 {
        self.seed
    }
}