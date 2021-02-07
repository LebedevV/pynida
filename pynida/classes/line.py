import numpy as np

from pynida.classes.point import Point
from pynida.classes.points import Points


class Line:

	@staticmethod
	def cosine_distance(lhs: 'Line', rhs: 'Line'):
		return lhs.cosine_distance(rhs)

	def __init__(self, p1: Point, p2: Point):
		# Equation: Ax + By + C = 0
		# https://en.wikipedia.org/wiki/Linear_equation#Two-point_form
		self.A = p1.y - p2.y
		self.B = p2.x - p1.x
		self.C = p1.x*p2.y - p2.x*p1.y

	# Equation: y = kx + beta
	@property
	def coeff_k(self):
		if self.is_vertical():
			raise AttributeError
		else:
			return -self.A / self.B

	# Equation: y = kx + beta
	@property
	def coeff_beta(self):
		if self.is_vertical():
			raise AttributeError
		else:
			return -self.C / self.B

	def is_vertical(self):
		return self.B == 0

	def is_horizontal(self):
		return self.A == 0

	def distances_to_points(self, points: Points):
		numerator = abs(self.A*points.x + self.B*points.y + self.C)
		denominator = np.sqrt(self.A**2 + self.B**2)
		return numerator / denominator

	def orthogonal_projection(self, points: Points):
		if self.is_vertical():
			return points.get_copy_with_x_changed_to(-self.C/self.A)
		elif self.is_horizontal():
			return points.get_copy_with_y_changed_to(self.coeff_beta)
		else:
			gamma = points.y + points.x/self.coeff_k
			c_x_l = (gamma - self.coeff_beta) / (self.coeff_k + 1/self.coeff_k)
			c_y_l = c_x_l*self.coeff_k + self.coeff_beta
			return Points.from_two_arrays(c_x_l, c_y_l)

	def cosine_distance(self, other: "Line"):
		# A - angle between lines
		# cosA = (p, q) / (|p| * |q|)
		# angle between lines is equal to angle between normal vectors to the lines
		# p and q - normal vectors to the lines
		p = Point(self.A, self.B)
		q = Point(other.A, other.B)

		p_abs = Point.get_distance_between(p, Point(0, 0))
		q_abs = Point.get_distance_between(q, Point(0, 0))

		dot_product = p.x*q.x + p.y*q.y
		cosA = dot_product / p_abs / q_abs
		return cosA

