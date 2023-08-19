use macroquad::prelude::{clear_background, draw_line, next_frame, LIGHTGRAY, BLUE, screen_height, screen_width, GREEN, RED, draw_circle, draw_text, DARKGRAY};
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
	let disc = (b*b - 4.0 * a * c).sqrt();
	[(-b + disc) / (2.0 * a), (-b - disc) / (2.0 * a)]
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
		let a = Vec2::new(0.05 * w, h/2.0);
		let d = Vec2::new(0.9 * w, h/2.0);
		let b = Lerp::lerp(a, d, 0.5) + Vec2::new(0.0, h / 2.0);
		let c = Lerp::lerp(a, d, 0.9) + Vec2::new(0.0, -h / 3.0);
        draw_circle(a.x, a.y, 5.0, BLUE);
        draw_circle(b.x, b.y, 5.0, BLUE);
        draw_circle(c.x, c.y, 5.0, BLUE);
        draw_circle(d.x, d.y, 5.0, BLUE);
		let bz = CubicBezier2 { start: a, ctrl0: b, ctrl1: c, end: d };

		let (qa, qb, qc) = ((-3.0*a + 3.0 * b - 3.0 * c + 3.0 * d), (6.0*a - 4.0*b + 2.0*c), (-3.0*a + b));
		let mut prev = a;
		let mut ks = Vec::new();
		for i in 1..=32 {
			let t = i as f32 / 32.0;
			let next = bz.evaluate(t);
			draw_circle(next.x, next.y, 2.5, GREEN);
			//let deriv = 40.0 * bz.normalized_tangent(i as f32 / 16.0);
			draw_line(prev.x, prev.y, next.x, next.y, 2.0, GREEN);
			let mid = Lerp::lerp(prev, next, 0.5);
			//draw_text(&format!("{:?}", bz.evaluate_derivative(i as f32/16.0)), mid.x, mid.y - 20.0, 10.0, DARKGRAY);
			//let f1 = t * t * qa + t * qb + qc;
			let sz = Vec2::new(w, h);
			let (a, b, c, d) = (a / sz, b / sz, c / sz, d / sz);
			let f1 = 3.0 * (1.0-t)*(1.0 - t)*(b - a) + 6.0*(1.0-t)*t*(c - b) + 3.0 * t * t * (d - c);
			let f2 = 6.0 * (1.0 - t) * (c - 2.0 * b + a) + 6.0 * t * (d - 2.0 * c + b);
			let k = (f1.magnitude_squared() * f2.magnitude_squared() - f1.dot(f2).powf(2.0)).sqrt()/f1.magnitude_squared().powf(3.0);
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
