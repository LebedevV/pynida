import numpy as np

class Point:
	def __init__(self, x=0, y=0):
		self.x = float(x)
		self.y = float(y)

	def __eq__(self, other: 'Point'):
		return self.x == other.x and self.y == other.y

	@staticmethod
	def from_tuple(t):
		return Point(t[0], t[1])

	@staticmethod
	def get_distance_between(lhs: 'Point', rhs: 'Point') -> float:
		delta_x = lhs.x - rhs.x
		delta_y = lhs.y - rhs.y
		return np.sqrt(delta_x**2 + delta_y**2)

	@staticmethod
	def get_middle_point(lhs: 'Point', rhs: 'Point') -> 'Point':
		x = (lhs.x + rhs.x) / 2
		y = (lhs.y + rhs.y) / 2
		return Point(x, y)
