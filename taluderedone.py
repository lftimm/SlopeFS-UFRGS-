import math
import sys
import numpy as np

# https://www.desmos.com/calculator/p9z9sikh4e

# Main model class,
# Responsible for handling data and doing math
class Model:
	def __init__(self):
		self.soilSlope = SoilSlope()
		self.circle = Circle(self.soilSlope)

		self.initModel()

	def initModel(self):
		self.points = self.intersec()
		self.splitgeometry()

	def intersec(self):
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
		a = math.tan(self.soilSlope.alp)
		pL,pR = self.points
		R,xc,yc = self.circle.geometry
		ns = self.soilSlope.num_slice

		v_pL = np.array(pL)
		v_pR = np.array(pR)
		v_c  = np.array([xc,yc])

		tot_angle = self.total_angle(pL,pR,xc,yc)		
		omg = tot_angle/ns # Slice of the angle and angular speed
		

		#print(vec_pL,vec_pR,vec_c)

		r = lambda t: R*np.array([math.cos(omg*t),math.sin(omg*t)])

		r1 = v_pL - v_c
		r2 = v_pR - v_c
		
		t0 = math.acos(r1[0]/R)/omg
		tf = math.acos(r2[0]/R)/omg
		
		#r_sections = np.array([r(t0+t) for t in range(dt)])
		#p_sections = np.array([(v_c-r) for r in r_sections])

		#print(p_sections)

	def total_angle(self,p1,p2,xc,yc):
		dist = lambda p1,p2: math.sqrt((p1[0]-p2[0])**2+(p1[1]+p2[1])**2)
		ct = (xc,yc)
		a,b,c = dist(p1,p2), dist(p1,ct), dist(p2,ct) 
		tot_angle = math.acos(-(a**2-(b**2+c**2))/(2*b*c))

		return tot_angle

	def __str__(self):
		string = ''
		string += f'Circle -> R: {self.circle.geometry[0]}, rx: {self.circle.geometry[1]}, ry: {self.circle.geometry[2]}\n'
		string += f'Slope -v\n'
		string += f'	c: {self.soilSlope.c}, phi: {self.soilSlope.phi}, gam: {self.soilSlope.gam}\n'
		string += f'	alp: {self.soilSlope.alp}, h: {self.soilSlope.h}, nf: {self.soilSlope.num_slice}\n'
		string += f'Intersect: {self.points}'
		return string


# Soil and Slope class,
# Not to be called outside the Model Class
# userInput gathers input from the user if False is passed
class SoilSlope:
	def __init__(self,shortcut=True):
		if shortcut:
			self.c,self.phi,self.gam,self.alp,self.h,self.num_slice,self.slope_len =[20,math.radians(25),18.5,math.radians(30),15,5,15/math.tan(math.radians(30))]
		else:
			self.userInput()

	# Might become obsolete if i use a GUI
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
			 

# Creates circle
# Circle is created with basis o the soilSlope given
# Includes magic numbers
class Circle:
	def __init__(self,parent: SoilSlope):
		self.xc = 0.5*parent.h/math.tan(parent.alp)
		self.yc = 1.67*parent.h 
		self.R = 1.5*math.sqrt(self.xc**2+self.yc**2)

		self.geometry = (self.R,self.xc,self.yc)

	def arc_length(self,angle):
		return angle*self.R

# Vizualization
# Maybe using matplotlib only, but might use wxPython together
class ModelView:
	def __init__(self):
		pass

def main():
	model = Model()

if __name__ == '__main__':
	main()