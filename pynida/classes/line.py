import numpy as np

from pynida.classes.point import Point
from pynida.classes.points import Points


class Line:
	def __init__(self, p1: Point, p2: Point):
		# Equation: Ax + By + C = 0
		# https://en.wikipedia.org/wiki/Linear_equation#Two-point_form
		self.A = p1.y - p2.y
		self.B = p2.x - p1.x
		self.C = p1.x*p2.y - p2.x*p1.y

	# Equation: y = kx + beta
	@property
	def coeff_k(self):
		if self.B != 0:
			return -self.A / self.B
		else:
			raise AttributeError

	# Equation: y = kx + beta
	@property
	def coeff_beta(self):
		if self.B != 0:
			return -self.C / self.B
		else:
			raise AttributeError

	def distances_to_points(self, points: Points):
		numerator = abs(self.A*points.x + self.B*points.y + self.C)
		denominator = np.sqrt(self.A**2 + self.B**2)
		return numerator / denominator

	def orthogonal_projection(self, points: Points):
		if self.B == 0:
			return points.get_copy_with_x_changed_to(-self.C/self.A)
		elif self.coeff_k == 0:
			return points.get_copy_with_y_changed_to(self.coeff_beta)
		else:
			gamma = points.y + points.x/self.coeff_k
			c_x_l = (gamma - self.coeff_beta) / (self.coeff_k + 1/self.coeff_k)
			c_y_l = c_x_l*self.coeff_k + self.coeff_beta
			return Points.from_two_arrays(c_x_l, c_y_l)
