import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
import math

# Define the size of the plot for cases where it isn't previously specified
x = y = range(100)


# Main field class, initializes one field with it's mean and standard variation
class Field:
	def __init__(self,vr,ln,xi=x,yi=y):
		self.field = self.fieldgen(vr,ln,xi,yi)
		self.mean = np.mean(self.field)
		self.stdvar = np.std(self.field)

	# Generates the field with size xi and yi
	def fieldgen(self,vr,ln,xi=x,yi=y):
		model = gs.Gaussian(dim=2, var=vr, len_scale=ln)
		srf = gs.SRF(model)
		field = srf((xi,yi), mesh_type='structured') 
		return np.array(field)
	
	# Prints te results
	def results(self):
		print({'Mean':round(self.mean,3),'Standard Variation':round(self.stdvar,3)})

# Creates a list of n plots, a list of the plots' means and variations and plots them
class FieldPlot:
	def __init__(self):
		self.plot_sizes = list(range(10,1000,10))
		self.fields = [Field(10,10,n,n) for n in self.plot_sizes]
		self.means = [field.mean for field in self.fields] 
		self.stdvrs = [field.stdvar for field in self.fields]

	def pltM(self):
		plt.plot(self.plot_sizes,self.means)
		plt.plot(self.plot_sizes,[sum(self.means)/len(self.means)]*len(self.plot_sizes))
		plt.show()

	def pltS(self):
		plt.plot(self.plot_sizes,self.stdvrs)
		plt.plot(self.plot_sizes,[sum(self.stdvrs)/len(self.stdvrs)]*len(self.plot_sizes))
		plt.show()

def main():
	plot = Field(10,10)

if __name__ == '__main__':
	main()
