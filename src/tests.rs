use std::{collections::HashMap, time::Instant};

use crate::Simplex;
#[derive(Hash, Eq, PartialEq)]

struct Value {
    integral: u64,
    fractional: u64,
}

type VMap2d = HashMap<(usize, usize), f32>;
type VMap3d = HashMap<(usize, usize, usize), f32>;

fn plot_2d_values(filename: &'static str, values: &VMap2d, n: usize) {
    // create image from values:
    let mut imgbuf = image::ImageBuffer::new(n as u32, n as u32);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let v = values[&(x as usize, y as usize)];
        let v = (v * 255.0) as u8;
        *pixel = image::Rgb([v, v, v]);
    }
    imgbuf.save(filename).unwrap();
}

fn plot_3d_values(filename: &'static str, values: &VMap3d, n: usize) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let root = BitMapBackend::gif(filename, (512, 512), 100)?.into_drawing_area();

    for z in 0..n {
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("3D Simplex Noise", ("sans-serif", 20))
            .build_cartesian_3d(0.0..n as f64, 0.0..n as f64, 0.0..n as f64)?;

        chart.with_projection(|mut p| {
            p.pitch = 0.0_f64.to_radians();
            p.yaw = 0.0_f64.to_radians();
            p.scale = 1.0;
            p.into_matrix() // build the projection matrix
        });

        chart.configure_axes().light_grid_style(BLACK.mix(0.15)).max_light_lines(3).draw()?;

        chart.draw_series(
            SurfaceSeries::xoy(
                (0..n).map(|x| x as f64),
                (0..n).map(|y| y as f64),
                |x, y| *values.get(&(x as usize, y as usize, z)).unwrap_or(&0.0) as f64
            )
            .style_func(&|&v| (VulcanoHSL::get_color(v)).into()),
        )?;

        root.present()?;
    }

    Ok(())
}

#[test]
fn test_noise2d_plot() {
    let mut values: VMap2d = HashMap::new();
    let start = Instant::now();
    let noise = Simplex::new(42);
    let freq = 0.007234;
    let n: usize = 512;
    for x in 0..n {
        for y in 0..n {
            let key = (x, y);
            let x = x as f32 * freq;
            let y = y as f32 * freq;
            let v = noise.sample2d(x, y);
            values.insert(key, v);
        }
    }
    let min = values.values().fold(f32::INFINITY, |a, &b| a.min(b));
    let max = values.values().fold(-f32::INFINITY, |a, &b| a.max(b));
    println!("min: {:?}", min);
    println!("max: {:?}", max);
    println!("avg: {:?}", values.values().sum::<f32>() / values.len() as f32);
    println!("duration: {:?} for {} values", start.elapsed(), values.len());

    // normalize value between 0.0 and 1.0 using the min and max of values:
    let values: VMap2d = values.iter().map(|(k, v)| (*k, (v - min) / (max - min))).collect();

    plot_2d_values("noise2d.png", &values, n);
}

#[test]
fn test_noise3d_plot() {
    let mut values: VMap3d = HashMap::new();
    let start = Instant::now();
    let noise = Simplex::new(42);
    let freq = 0.02234;
    let n: usize = 100;
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let key = (x, y, z);
                let x = x as f32 * freq;
                let y = y as f32 * freq;
                let z = z as f32 * freq;
                let v = noise.sample3d(x, y, z);
                values.insert(key, v);
            }
        }
    }
    let min = values.values().fold(f32::INFINITY, |a, &b| a.min(b));
    let max = values.values().fold(-f32::INFINITY, |a, &b| a.max(b));
    println!("min: {:?}", min);
    println!("max: {:?}", max);
    println!("avg: {:?}", values.values().sum::<f32>() / values.len() as f32);
    println!("duration: {:?} for {} values", start.elapsed(), values.len());

    // normalize value between 0.0 and 1.0 using the min and max of values:
    let values: VMap3d = values.iter().map(|(k, v)| (*k, (v - min) / (max - min))).collect();

    plot_3d_values("noise3d.gif", &values, n).unwrap();
}
