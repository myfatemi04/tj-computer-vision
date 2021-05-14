struct Point {
	x: f64,
	y: f64,
	z: f64
}

impl std::fmt::Display for Point {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		write!(f, "({}, {}, {})", self.x, self.y, self.z)
	}
}

struct Tetrahedron<'object> {
	a: &'object Point,
	b: &'object Point,
	c: &'object Point,
	d: &'object Point
}

impl<'object> Tetrahedron<'object> {
	fn new(a: &'object Point, b: &'object Point, c: &'object Point, d: &'object Point) -> Tetrahedron<'object> {
		Tetrahedron {
			a, b, c, d
		}
	}
}

impl std::fmt::Display for Tetrahedron<'_> {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		write!(f, "Tetrahedron {{{} {} {} {}}}", self.a, self.b, self.c, self.d)
	}
}

struct Object<'object> {
	components: Vec<Tetrahedron<'object>>
}

impl<'object> Object<'object> {
	fn empty() -> Object<'object> {
		Object {
			components: Vec::new()
		}
	}

	fn add_component(&mut self, component: Tetrahedron<'object>) {
		self.components.push(component);
	}
}

impl std::fmt::Display for Object<'_> {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		write!(f, "Object {{components: [")?;

		for component in &self.components {
			write!(f, "{} ", component)?;
		}

		write!(f, "]")
	}
}

fn main() {
	let d = Point { x: 0.0, y: 0.0, z: 3.0 };
	let e = Point { x: 1.0, y: 0.0, z: 2.0 };
	let f = Point { x: 2.0, y: 0.0, z: 1.0 };
	let g = Point { x: 3.0, y: 0.0, z: 0.0 };

	let a = Tetrahedron::new(&d, &e, &f, &g);
	let b = Tetrahedron::new(&d, &e, &f, &g);

	let mut o = Object::empty();
	o.add_component(a);
	o.add_component(b);

	println!("{}", o);
}
