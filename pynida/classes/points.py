import numpy as np

from pynida.classes.point import Point

class Points:

	@staticmethod
	def from_two_arrays(X, Y):
		return Points([Point(x, y) for x, y in zip(X, Y)])

	def __init__(self, points=[]):
		self.points = list(points)

	def __getitem__(self, item):
		return self.points[item]

	@property
	def x(self):
		return np.array([p.x for p in self.points])

	@property
	def y(self):
		return np.array([p.y for p in self.points])

	def append(self, point: Point):
		self.points.append(point)

	def __iter__(self):
		return iter(self.points)

	def get_copy_with_x_changed_to(self, value):
		return Points([Point(value, p.y) for p in self])

	def get_copy_with_y_changed_to(self, value):
		return Points([Point(p.x, value) for p in self])

	def get_distances_to(self, p: Point):
		return np.array([Point.get_distance_between(point, p) for point in self.points])
