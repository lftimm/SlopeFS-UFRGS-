import math
import numpy as np
import matplotlib.pyplot as plt

PI = math.pi

def main():
	C = np.array([0,5])
	pL = np.array([-1,5])
	pR = np.array([1,5])
	num_slice = 5
	
	dist = lambda p1,p2: math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
	R = dist(pL,C)

	a,b,c = dist(pL,pR), dist(pL,C), dist(pR,C) 
	tot_angle = math.acos(-(a**2-(b**2+c**2))/(2*b*c))

	alp = tot_angle/num_slice
	gam = math.atan(abs((C[1]-pR[1])/(C[0]-pR[0])))

	f = lambda n: R*np.array([math.cos(-(gam+n*alp)),math.sin(-(gam+n*alp))])+C
	pos = [tuple(f(n)) for n in range(num_slice+1)]
	print(pos)
	

if __name__ == '__main__':
	main()