import math
import numpy as np

"""
This script has the purpose of modeling a real slope situation.
It is divided into two parts. Model refers to the handling of information and maths.
View is used for displaying the values calculated.

Desmos link for a visual representation of the problem.
https://www.desmos.com/calculator/p9z9sikh4e
"""

""" TODO:
	- Calc FS Fellenius
	- Calc FS Bishop
	- Exemplo de Talude do Artigo, calcular o FS dele.

	- Graph FS
"""

class Model:
	"""
	Model class
	-----------
	It functions as the brain of the program.
	It performs multiple geometric calculations that serve the purpose of the project.
	Together with it there are the SoilSlope() and Circle() class that make up the model for the slope.
	
	Review of the main methods:
	__init__ -> defines the slope, circle and begins the model.
	initModel -> Organizes all the processes needed for the program and its results.
	"""
	def __init__(self):
		self.soilSlope = SoilSlope()
		self.circle = Circle(self.soilSlope)

		self.initModel()

	def initModel(self):
		self.points = self.intersec()
		self.c_points = self.splitgeometry()
		self.polys = self.mk_polys()
		##############
		self.polys_A = self.calc_areas()

	def intersec(self):
		"""
			Calculates the points of intersection of the slope and the circle.
			It uses second degree equations to find the intersections.
			It returns the points in order from left to right.
		"""
		t = math.tan(self.soilSlope.alp)
		h = self.soilSlope.h
		l = self.soilSlope.slope_len
		R,xc,yc = self.circle.geometry

		a = 1 + t**2
		b1 = -2*xc
		b = b1 -2*t*yc
		c = xc**2 + yc**2 - R**2

		delta = b**2 - 4*a*c

		if delta > 0:
			def f(x):
				if 0 < x < l:
					return (x,t*x)
				elif x < 0:
					return (xc-math.sqrt(R**2-yc**2),0)
				else:
					return (xc+math.sqrt(R**2+2*h*yc-h**2-yc**2),h)

			x1 = (-b + math.sqrt(delta))/(2*a)
			x2 = (-b - math.sqrt(delta))/(2*a)
			p1 = f(x1)
			p2 = f(x2)

			pL = p1 if p1[0] == min(p1[0],p2[0]) else p2
			pR = p2 if pL[0] == p1[0] else p1

			return pL,pR

	def splitgeometry(self):
		"""
			It splits the circle into equal parts based on the number of slices given.
			It uses a trigonometrical approach.
			Together there is the total_angle method that, as the name says, calculates the total angle of the intersected points.
			It returns a list of tuples containing the points.
			One thing might be removed, in the definition of f() there is a rounding done with the map() function.
		"""
		a = math.tan(self.soilSlope.alp)
		pL,pR = self.points
		R,xc,yc = self.circle.geometry
		ns = self.soilSlope.num_slice

		v_pL = np.array(pL)
		v_pR = np.array(pR)
		v_c  = np.array([xc,yc])

		tot_a = self.total_angle(pL,pR,xc,yc)		
		alp = tot_a/ns 
		gam = math.atan(abs((v_c[1]-v_pR[1])/(v_c[0]-v_pR[0])))

		f = lambda n: map(lambda x: round(x,2),R*np.array([math.cos(-(gam+n*alp)),math.sin(-(gam+n*alp))])+v_c)
		return [tuple(f(n)) for n in range(ns+1)]

	def total_angle(self,p1,p2,xc,yc):
		dist = lambda p1,p2: math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
		ct = (xc,yc)
		a,b,c = dist(p1,p2), dist(p1,ct), dist(p2,ct) 
		tot_angle = math.acos(-(a**2-(b**2+c**2))/(2*b*c))

		return tot_angle

	def mk_polys(self):
		"""
			This method creates the polygons whose areas are going to be calculated.
			It takes the list of points in the circle, reflects them into the corresponding part of the surface.
			Together with it there is the pair_points method, which takes the previous points and orders them counter-clockwise.
			It returns an array of 4x1 arrays with the points .
		"""
		c_parts = self.c_points
		pts_x,pts_y = zip(*c_parts)
		a = self.soilSlope.alp
		h = self.soilSlope.h
		l = self.soilSlope.slope_len

		def f(x):
			if 0<= x<= l:
				return round(a*x,2)
			elif x>l:
				return h
			else:
				return 0

		up_c_parts =[(x,f(x)) for x in pts_x]
		
		full_points = [c_parts,up_c_parts]
		
		return self.pair_points(full_points)

	def pair_points(self,points):
		polys = []
		low, up = 0, 1
		for i in range(0,len(points[0])):
			try:
				polys.append([points[0][i],points[1][i],points[1][i+1],points[0][i+1]])
			except IndexError:
				pass

		return np.array(polys)

	def calc_areas(self):
		"""
			It calculates the areas of the polygons.
			It uses the shoelace formula for calculating the area.
			For this purpose, it takes the 4x1 arrays, orders than and takes their dertimant/2, summing the areas.
			It returns an array containging the areas of each of the polygons.
		"""
		p = self.polys

		areas = []
		
		for poly in p:
			n = len(poly)
			a = 0
			for i in range(n):
				x1,y1 = poly[i]
				x2,y2 = poly[(i+1)%n]
				a += (x1*y2-x2*y1)/2
			areas.append(a)

		return areas

"""
	Class representing the soil and the slope.
	userInput gathers input from the user if False is passed.
	contains magic numbers that should be removed.
"""
class SoilSlope:
	def __init__(self,shortcut=True):
		if shortcut:
			self.c,self.phi,self.gam,self.alp,self.h,self.num_slice,self.slope_len =[20,math.radians(25),18.5,math.radians(30),15,5,15/math.tan(math.radians(30))]
		else:
			self.userInput()

	def userInput(self):
		while True:
			try:
				self.c   = float(input('Soil cohesion(kN/m²): '))  
				self.phi = float(input('Attrition angle(º): '))
				self.gam = float(input('Specific Weight(kN/m³): '))
				self.alp = float(input('Slope angle(º): '))
				self.h   = float(input('Slope Height(m): '))
				self.num_slice  = int(input('Number of sections: '))
				self.slope_len	= self.h/math.tan(self.alp)
				break
			except ValueError:
				print('Please, insert correct values')
			 
"""
	Class that defines a circle given a SoilSlope.
	contains magic numbers that should be removed.
"""
class Circle:
	def __init__(self,parent: SoilSlope):
		self.xc = 0.5*parent.h/math.tan(parent.alp)
		self.yc = 1.67*parent.h 
		self.R = 1.5*math.sqrt(self.xc**2+self.yc**2)

		self.geometry = (self.R,self.xc,self.yc)

"""
	TODO:
	Class that makes it possible to view the contents of the model.
	Possibilites: matplotlib or wxpython, for a GUI.
"""
class ModelView:
	def __init__(self):
		pass

def main():
	model = Model()

if __name__ == '__main__':
	main()