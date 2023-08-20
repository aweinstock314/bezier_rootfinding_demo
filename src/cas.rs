use std::ops::{Add, Mul, Neg, Sub};
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
