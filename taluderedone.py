import matplotlib.pyplot as plt
import numpy as np
import wx
import math
import typing 

# Main Controller Class
class Controller(wx.Frame):
	def __init__(self):
		super().__init__(None,title='Talude')
		
		self.Center()
		self.Show()

# Main model class,
# Responsible for data handling
class Model:
	def __init__(self):
		self.soilSlope = SoilSlope()
		self.circle = Circle(self.soilSlope)
		self.inter = self.intersec()

		ModelView(self)

	def intersec(self):
		pass

# Soil class,
# Not to be called outside the Model Class
# userInput gathers input if False is passed
class Soil:
	def __init__(self,shortcut=True):
		if shortcut:
			self.c,self.phi,self.gam,self.alp,self.h = [20,
														math.radians(25),
														18.5,
														math.radians(30),
														15]
		else:
			self.c,self.phi,self.gam,self.alp,self.h = self.userInput()

	def userInput(self):

		pass

# Creates circle
class Circle:
	def __init__(self,soil: Soil):
		pass

class ModelView:
	def __init__(self):
		pass

def main():
	model = Controller()

if __name__ == '__main__':
	main()