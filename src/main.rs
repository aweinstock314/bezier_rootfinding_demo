use macroquad::prelude::{
    clear_background, draw_circle, draw_line, draw_text, next_frame, screen_height, screen_width,
    BLUE, DARKGRAY, GREEN, LIGHTGRAY, RED,
};
use std::ops::{Add, Mul, Neg, Sub};
use vek::*;

/*

f(t) = (1-t)^3 A + 3(1-t)^2 t B + 3(1-t)t^2 C + t^3 D
f'(t) = 3(1-t)^2 (B - A) + 6(1-t)t(C - B) + 3t^2 (D - C)
f''(t) = 6(1-t)(C - 2B + A) + 6t(D - 2C + B)
f'''(t) = -6(a - b - c + d)
---
f'(t) = 3(1-t)^2 (B - A) + 6(1-t)t(C - B) + 3t^2 (D - C) = 3(1 - t)
>>> f1.expand()
-3*a*t**2 + 6*a*t - 3*a + 3*b*t**2 - 4*b*t + b - 3*c*t**2 + 2*c*t + 3*d*t**2
---
-3*a*t**2 + 3*b*t**2 - 3*c*t**2 + 3*d*t**2 + 6*a*t - 4*b*t + 2*c*t - 3*a + b
(-3*a + 3*b - 3*c + 3*d)*t**2 + (6*a - 4*b + 2*c)*t - 3*a + b
*/
fn quadform(a: f32, b: f32, c: f32) -> [f32; 2] {
    //println!("{:?} {:?} {:?}", a, b, c);
    let disc = (b * b - 4.0 * a * c).sqrt();
    [(-b + disc) / (2.0 * a), (-b - disc) / (2.0 * a)]
}

/*
>>> f = (1-t)**3 * a + (1-t)**2*t * b + (1-t)*t**2 * c + t**3 * d
>>> f
a*(1 - t)**3 + b*t*(1 - t)**2 + c*t**2*(1 - t) + d*t**3
>>> f.as_poly(t)
Poly((-a + b - c + d)*t**3 + (3*a - 2*b + c)*t**2 + (-3*a + b)*t + a, t, domain='ZZ[a,b,c,d]')
*/

#[derive(Clone, Debug)]
enum SFormula {
    T,
    F(f32),
    Pow(Box<SFormula>, f32),
    SoP(Vec<Vec<SFormula>>),
    Dot(Box<(VFormula, VFormula)>),
}

#[derive(Clone, Debug)]
enum VFormula {
    P(usize),
    SMul(SFormula, Box<VFormula>),
    Sum(Vec<VFormula>),
}

impl SFormula {
    fn pow(&self, x: f32) -> SFormula {
        SFormula::Pow(Box::new(self.clone()), x)
    }
    fn simplify(self) -> SFormula {
        self
    }
}

impl VFormula {
    fn dot(&self, other: &Self) -> SFormula {
        SFormula::Dot(Box::new((self.clone(), other.clone())))
    }
    fn flatten(self) -> VFormula {
        match self {
            VFormula::P(i) => VFormula::P(i),
            VFormula::SMul(s, f) => {
                let g = f.flatten();
                if let VFormula::Sum(gs) = g {
                    VFormula::Sum(
                        gs.into_iter()
                            .map(|g| VFormula::SMul(s.clone(), Box::new(g)).flatten())
                            .collect(),
                    )
                } else {
                    VFormula::SMul(s, Box::new(g))
                }
            }
            VFormula::Sum(fs) => {
                let mut ret = Vec::new();
                for f in fs {
                    let g = f.flatten();
                    if let VFormula::Sum(gs) = g {
                        ret.extend(gs)
                    } else {
                        ret.push(g)
                    }
                }
                VFormula::Sum(ret)
            }
        }
    }
    fn simplify(self) -> VFormula {
        match self {
            VFormula::P(i) => VFormula::P(i),
            VFormula::SMul(s, f) => {
                let t = s.simplify();
                if let SFormula::F(1.0) = t {
                    f.simplify()
                } else {
                    VFormula::SMul(t, Box::new(f.simplify()))
                }
            }
            VFormula::Sum(fs) => {
                let mut ret = Vec::new();
                for f in fs {
                    let g = f.simplify();
                    match g {
                        VFormula::P(i) => ret.push(VFormula::P(i)),
                        VFormula::SMul(s, f) => {
                            if let SFormula::F(0.0) = s {
                            } else {
                                ret.push(VFormula::SMul(s, Box::new(f.simplify())));
                            }
                        }
                        VFormula::Sum(gs) => {
                            for g in gs {
                                let h = g.simplify();
                                ret.push(h);
                            }
                        }
                    }
                }
                VFormula::Sum(ret)
            }
        }
    }
}

impl Neg for VFormula {
    type Output = VFormula;
    fn neg(self) -> VFormula {
        VFormula::SMul(SFormula::F(-1.0), Box::new(self))
    }
}

impl Add<VFormula> for VFormula {
    type Output = VFormula;
    fn add(self, other: VFormula) -> VFormula {
        VFormula::Sum(vec![self, other]).flatten()
    }
}
impl Mul<SFormula> for VFormula {
    type Output = VFormula;
    fn mul(self, other: SFormula) -> VFormula {
        VFormula::SMul(other, Box::new(self))
    }
}
impl Mul<VFormula> for f32 {
    type Output = VFormula;
    fn mul(self, other: VFormula) -> VFormula {
        VFormula::SMul(SFormula::F(self), Box::new(other))
    }
}
impl Sub<VFormula> for VFormula {
    type Output = VFormula;
    fn sub(self, other: VFormula) -> VFormula {
        self + (-other)
    }
}

#[test]
fn test_derivative() {
    use SFormula::T;
    use VFormula::P;
    let f = (-P(0) + P(1) - P(2) + P(3)) * T.pow(3.0)
        + (3.0 * P(0) - 2.0 * P(1) + P(2)) * T.pow(2.0)
        + (-3.0 * P(0) + P(1)) * T;
    println!("{:?}", f);
    println!("{:?}", f.simplify());
}

const EPSILON: f32 = 0.0001;

fn secant<F: Fn(f32) -> f32>(f: &F, mut x0: f32, mut x1: f32) -> f32 {
    let mut y0 = f(x0);
    let mut y1 = f(x1);
    let mut x2 = x1 - y1 * (x1 - x0) / (y1 - y0);
    while y1.abs() > EPSILON {
        x0 = x1;
        x1 = x2;
        y0 = f(x0);
        y1 = f(x1);
        x2 = x1 - y1 * (x1 - x0) / (y1 - y0);
    }
    x2
}

fn newton_raphson<F: Fn(f32) -> f32, F1: Fn(f32) -> f32>(f: &F, f1: &F1, mut x0: f32) -> f32 {
    let mut y0 = f(x0);
    let mut x1 = x0 - y0 / f1(x0);
    while y0.abs() > EPSILON {
        x0 = x1;
        y0 = f(x0);
        x1 = x0 - y0 / f1(x0);
    }
    x1
}

#[test]
fn test_rootfinding() {
    let f = |x: f32| -x * x + 3.0 * x + 5.0;
    let f1 = |x: f32| -2.0 * x + 3.0;
    let f2 = |x: f32| -2.0;
    let x = secant(&f, 0.0, 1.0);
    println!("secant f: {:?} {:?} {:?}", x, f(x), f1(x));
    let x1 = secant(&f1, 0.0, 1.0);
    println!("secant f1: {:?} {:?} {:?}", x1, f(x1), f1(x1));
    let x2 = newton_raphson(&f, &f1, 1.0);
    println!("nr f f1: {:?} {:?} {:?}", x2, f(x2), f1(x2));
    let x3 = newton_raphson(&f1, &f2, 1.0);
    println!("nr f1 f2: {:?} {:?} {:?}", x3, f(x3), f1(x3));
}

fn gradient_descent<F: FnMut(f32) -> f32>(f: &mut F, alpha: f32, mut x0: f32) -> f32 {
    let mut y0 = f(x0);
    loop {
        let grad = (f(x0 + EPSILON) - y0) / EPSILON;
        let x1 = x0 - alpha * grad;
        let y1 = f(x1);
        if grad == 0.0 || y0 < y1 {
            return x0;
        }
        x0 = x1;
        y0 = y1;
    }
}

#[test]
fn test_gradientdescent() {
    let f = |x: f32| -x * x + 3.0 * x + 5.0;
    // Objective is negated because gradient_descent finds minima and we're trying to find a maximum
    let x = gradient_descent(&mut |x| -f(x), 0.1, 1.0);
    println!("gradient_descent: {} {}", x, f(x));
}

#[test]
fn test_gradientdescent_bezier() {
    let w = 1.0;
    let h = 1.0;
    let a = Vec2::new(0.05 * w, h / 2.0);
    let d = Vec2::new(0.9 * w, h / 2.0);
    let b = Lerp::lerp(a, d, 0.5) + Vec2::new(0.0, h / 2.0);
    let c = Lerp::lerp(a, d, 0.9) + Vec2::new(0.0, -h / 3.0);

    let f = |t: f32| {
        let f1 = 3.0 * (1.0 - t) * (1.0 - t) * (b - a)
            + 6.0 * (1.0 - t) * t * (c - b)
            + 3.0 * t * t * (d - c);
        let f2 = 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
        let k = (f1.magnitude_squared() * f2.magnitude_squared() - f1.dot(f2).powf(2.0)).sqrt()
            / f1.magnitude_squared().powf(3.0);
        k
    };
    let t0 = gradient_descent(&mut |x| -f(x), 0.001, 0.0);
    println!("gradient_descent_bezier: {} {}", t0, f(t0));
    let t1 = gradient_descent(&mut |x| -f(x), 0.001, 1.0);
    println!("gradient_descent_bezier: {} {}", t1, f(t1));
}

#[test]
fn test_bezier() {
    let a = Vec2::new(0.0, 0.0);
    let b = Vec2::new(1.0, 0.0);
    let c = Vec2::new(1.0, 1.0);
    let d = Vec2::new(0.0, 1.0);
    let f = |t: f32| {
        (1.0 - t).powf(3.0) * a
            + 3.0 * (1.0 - t).powf(2.0) * t * b
            + 3.0 * (1.0 - t) * t.powf(2.0) * c
            + t.powf(3.0) * d
    };
    let f1 = |t: f32| {
        3.0 * (1.0 - t) * (1.0 - t) * (b - a)
            + 6.0 * (1.0 - t) * t * (c - b)
            + 3.0 * t * t * (d - c)
    };
    let f2 = |t: f32| 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
    println!("f(0) = {}", f(0.0));
    println!("f'(0) = {}", f1(0.0));
    println!("f''(0) = {}", f2(0.0));
    println!("f(1) = {}", f(1.0));
    println!("f'(1) = {}", f1(1.0));
    println!("f''(1) = {}", f2(1.0));
}

#[macroquad::main("BasicShapes")]
async fn main() {
    loop {
        clear_background(LIGHTGRAY);

        /*draw_line(40.0, 40.0, 100.0, 200.0, 15.0, BLUE);
        draw_rectangle(screen_width() / 2.0 - 60.0, 100.0, 120.0, 60.0, GREEN);
        draw_circle(screen_width() - 30.0, screen_height() - 30.0, 15.0, YELLOW);

        draw_text("HELLO", 20.0, 20.0, 30.0, DARKGRAY);*/

        let (w, h) = (screen_width(), screen_height());
        let a = Vec2::new(0.05 * w, h / 2.0);
        let d = Vec2::new(0.9 * w, h / 2.0);
        let b = Lerp::lerp(a, d, 0.5) + Vec2::new(0.0, h / 2.0);
        let c = Lerp::lerp(a, d, 0.9) + Vec2::new(0.0, -h / 3.0);
        draw_circle(a.x, a.y, 5.0, BLUE);
        draw_circle(b.x, b.y, 5.0, BLUE);
        draw_circle(c.x, c.y, 5.0, BLUE);
        draw_circle(d.x, d.y, 5.0, BLUE);
        let bz = CubicBezier2 {
            start: a,
            ctrl0: b,
            ctrl1: c,
            end: d,
        };

        let (qa, qb, qc) = (
            (-3.0 * a + 3.0 * b - 3.0 * c + 3.0 * d),
            (6.0 * a - 4.0 * b + 2.0 * c),
            (-3.0 * a + b),
        );
        let mut prev = a;
        let mut ks = Vec::new();
        for i in 1..=64 {
            let t = i as f32 / 64.0;
            let next = bz.evaluate(t);
            draw_circle(next.x, next.y, 2.5, GREEN);
            //let deriv = 40.0 * bz.normalized_tangent(i as f32 / 16.0);
            draw_line(prev.x, prev.y, next.x, next.y, 2.0, GREEN);
            let mid = Lerp::lerp(prev, next, 0.5);
            //draw_text(&format!("{:?}", bz.evaluate_derivative(i as f32/16.0)), mid.x, mid.y - 20.0, 10.0, DARKGRAY);
            //let f1 = t * t * qa + t * qb + qc;
            let sz = Vec2::new(w, h);
            let (a, b, c, d) = (a / sz, b / sz, c / sz, d / sz);
            let f1 = 3.0 * (1.0 - t) * (1.0 - t) * (b - a)
                + 6.0 * (1.0 - t) * t * (c - b)
                + 3.0 * t * t * (d - c);
            let f2 = 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
            let k = (f1.magnitude_squared() * f2.magnitude_squared() - f1.dot(f2).powf(2.0)).sqrt()
                / f1.magnitude_squared().powf(3.0);
            ks.push((t, k));
            //draw_text(&format!("f': {:0.2} {:0.2}", f1.x, f1.y), next.x, next.y - 40.0, 10.0, DARKGRAY);
            //draw_text(&format!("f'': {:0.2} {:0.2}", f2.x, f2.y), next.x, next.y - 30.0, 10.0, DARKGRAY);
            //draw_text(&format!("k: {:0.4}", k), next.x, next.y - 20.0, 10.0, DARKGRAY);
            prev = next;
        }
        for w in ks.windows(3) {
            if w[1].1 > w[0].1 && w[1].1 > w[2].1 {
                let t = w[1].0;
                let f = bz.evaluate(t);
                draw_circle(f.x, f.y, 2.5, RED);
            }
        }
        let f = |t: f32| {
            let sz = Vec2::new(w, h);
            let (a, b, c, d) = (a / sz, b / sz, c / sz, d / sz);
            let f1 = 3.0 * (1.0 - t) * (1.0 - t) * (b - a)
                + 6.0 * (1.0 - t) * t * (c - b)
                + 3.0 * t * t * (d - c);
            let f2 = 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
            let k = (f1.magnitude_squared() * f2.magnitude_squared() - f1.dot(f2).powf(2.0)).sqrt()
                / f1.magnitude_squared().powf(3.0);
            k
        };
        let mut steps0 = 0;
        let t0 = gradient_descent(
            &mut |x| {
                steps0 += 1;
                -f(x)
            },
            0.00001,
            0.0,
        );
        let mut steps1 = 0;
        let t1 = gradient_descent(
            &mut |x| {
                steps1 += 1;
                -f(x)
            },
            0.00001,
            1.0,
        );
        let p0 = bz.evaluate(t0);
        let p1 = bz.evaluate(t1);
        draw_circle(p0.x, p0.y, 2.5, DARKGRAY);
        draw_circle(p1.x, p1.y, 2.5, DARKGRAY);
        draw_text(
            &format!("gradient descent steps: {} {}", steps0, steps1),
            20.0,
            20.0,
            30.0,
            DARKGRAY,
        );
        let g = |t: f32| {
            let sz = Vec2::new(w, h);
            let (a, b, c, d) = (a / sz, b / sz, c / sz, d / sz);
            let f1 = 3.0 * (1.0 - t) * (1.0 - t) * (b - a)
                + 6.0 * (1.0 - t) * t * (c - b)
                + 3.0 * t * t * (d - c);
            let f2 = 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
            let f3 = -6.0 * (a - b - c + d);
            let k = (f2.magnitude_squared() * f3.magnitude_squared() - f2.dot(f3).powf(2.0)).sqrt()
                / f2.magnitude_squared().powf(3.0);
            k
        };
        /*let t2 = secant(&|x: f32| { let e = 0.01; (f(x + e/2.0) - f(x - e/2.0)) / e }, 0.5, 0.6);
        let t2 = secant(&g, 0.9, 1.0);
        let p2 = bz.evaluate(t2);
        draw_text(&format!("secant t: {} {:?}", t2, p2), 20.0, 20.0, 30.0, DARKGRAY);
        draw_circle(p2.x, p2.y, 2.5, BLUE);*/

        /*let mut prev = a;
        for i in 1..=16 {
            let next = bz.evaluate(i as f32 / 16.0);
            let deriv = 40.0 * bz.normalized_tangent(i as f32 / 16.0);
            draw_line(next.x, next.y, next.x + deriv.x, next.y + deriv.y, 2.0, RED);
            prev = next;
        }*/

        let ts = quadform(qa.magnitude(), qb.magnitude(), qc.magnitude());
        //draw_text(&format!("{:?}", (qa, qb, qc)), 20.0, 20.0, 30.0, DARKGRAY);
        //draw_text(&format!("{:?}", ts), 20.0, 40.0, 30.0, DARKGRAY);

        next_frame().await
    }
}
