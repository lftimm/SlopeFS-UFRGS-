import matplotlib.pyplot as plt
import scipy.optimize as sp
from math import *
import numpy as np

#função objetivo
def f(x):
	f = calcula_FS(x)
	return f
	
#equação transcendental de FS de Bishop:
def fs(x):
	global R,rx,ry,α,h,nf
	xya,xyb,intercecao_valida = intersecao_circulo_com_talude(R,rx,ry,α,h)
	θf,vxf,vyf,Δxf = geometria_fatias_c(R,rx,ry,α,h,nf,xya,xyb)
	F1 = 0.0
	F2 = 0.0
	for i in range(0,nf):
		A = area_poligono(vxf[i],vyf[i])
		W = A*γ
		F1 += W*sin(θf[i]*(pi/180.0))
		F2 += ( Δxf[i]*c + W*tan(ϕ) ) / ( cos( θf[i]*(pi/180.0) ) + sin( θf[i]*(pi/180.0) )*tan(ϕ)/x   )
	fs = x - F2/F1
	return fs	

#equação transcendental de FS de Janbu:
def fj(x):
	global xa,xb,yi,α,h,nf
	ya,yb,intercecao_valida = intersecao_poligonal_com_talude(xa,xb,yi,α,h)
	θf,vxf,vyf,Δxf = geometria_fatias_p(xa,xb,yi,α,h,nf)
	F1 = 0.0
	F2 = 0.0
	for i in range(0,nf):
		A = area_poligono(vxf[i],vyf[i])
		W = A*γ
		F1 += W*tan(θf[i]*(pi/180.0))
		F2 += ( Δxf[i]*c + W*tan(ϕ) )  / ( ((cos( θf[i]*(pi/180.0) ))**2)*(1.0+(tan(θf[i]*(pi/180.0)))*(tan(ϕ)/x)) )  
	fj = x - F2/F1
	return fj
	
def intersecao_circulo_com_talude(R,rx,ry,α,h):
	
	# ________________________________________
	#
	# ALGORITMO DE INTERSEÇÃO CIRCULAR
	# ________________________________________
	#
	# > R        : raio do círculo de ruptura
	# > rx       : coordenada x do centro do círculo
	# > ry       : coordenada y do centro do círculo
	# > α        : inclinação do talude
	# > h        : altura do talude
	# < xya      : [x,y] da inteseção à esquerda
    # < xyb      : [x,y] da inteseção à direita
   	# < verifica : True ou False se interseção existe
	# ________________________________________
	
	xya = []
	xyb = []
#linha horizontal esquerda
	rad = (R**2-ry**2)
	if ( rad >= 0.0 ):
		intercept = rx-sqrt(rad)
		if ( intercept <= 0.0 ):
			xya = [intercept,0.0]
#linha inclinada central
	aa = (1.0 + tan(α)**2)
	bb = -2.0 * ( rx + tan(α)*ry )
	cc = rx**2 + ry**2 - R**2
	Δ = bb**2 - 4.0*aa*cc
	lx = h/tan(α)
	if ( Δ > 0.0 ):
		xx = ( -bb - sqrt(Δ) ) / ( 2.0 * aa )
		if ( xx >= 0.0 and xx <= lx ):
			xya = [ xx , xx*tan(α) ]
		xx = ( -bb + sqrt(Δ) ) / ( 2.0 * aa )
		if ( xx >= 0.0 and xx <= lx ):
			xyb = [ xx , xx*tan(α) ]
#linha horizontal direita
	rad = (R**2-(h-ry)**2)
	if ( rad >= 0.0 ):
		intercept = rx+sqrt(rad)
		if ( intercept >= lx ):
			xyb = [intercept,h]
#verifica se a interseção por dois pontos é válida
	if ( xya and xyb):
		intercecao_valida = True
		if ( xyb[1] > ry ) :
			intercecao_valida = False
	else:
		intercecao_valida = False
	return xya,xyb,intercecao_valida

def interpolacao_coordenadas_talude(xp,α,h):
	
	# ________________________________________
	#
	# ALGORITMO DE DETERMINAÇÃO DE COORDENADAS
	# NO TALUDE
	# ________________________________________
	#
	# > xp       : coordenada x 
	# > α        : inclinação do talude
	# > h        : altura do talude
	# < yp       : coordenada y
	# ________________________________________
	
	lx = h/tan(α)
#determina o y correspondente ao valor x pertencente ao talude
	if ( xp <= 0.0 ):
		yp = 0.0
	elif ( xp > 0.0 and xp < lx ):
		yp = xp*tan(α)
	else:
		yp = h
	return yp
	
def intersecao_poligonal_com_talude(xa,xb,yi,α,h):
	
	# ________________________________________
	#
	# ALGORITMO DE INTERSEÇÃO POLIGONAL
	# ________________________________________
	#
	# > xa       : x da inteseção à esquerda
	# > xb       : x da inteseção à direita
	# > yi       : lista de coordenadas y dos pontos entre xa e xb
	# > α        : inclinação do talude
	# > h        : altura do talude
	# < ya       : y da inteseção à esquerda
	# < yb       : y da inteseção à direita	
   	# < verifica : True ou False se interseção existe
	# ________________________________________
	
#determina o y correspondente ao valor x das extremidades da poligonal	
	verifica = True
	ya = interpolacao_coordenadas_talude(xa,α,h)
	yb = interpolacao_coordenadas_talude(xb,α,h)
	lyi = len(yi)
	dx = (xb-xa)/(lyi+1.0)
	xi = xa
#verifica se poligonal não intercepta o talude
	for i in range(0,lyi):
		xi += dx
		if ( interpolacao_coordenadas_talude(xi,α,h) < yi[i] ):
			verifica = False
	return ya,yb,verifica
	
def geometria_fatias_c(R,rx,ry,α,h,nf,xya,xyb):
	
	# ________________________________________
	#
	# ALGORITMO PARA OBTER GEOMETRIA DE FATIAS
	# CIRCULARES
	# ________________________________________
	#
	# > R        : raio do círculo de ruptura
	# > rx       : coordenada x do centro do círculo
	# > ry       : coordenada y do centro do círculo
	# > α        : inclinação do talude
	# > h        : altura do talude
	# > nf       : número de fatias
	# > xya      : [x,y] da inteseção à esquerda
    # > xyb      : [x,y] da inteseção à direita
   	# < θf       : lista com ângulo central de cada fatia medido em graus
   	# < vxf      : lista ordenada de coordenadas x de cada vertice de cada fatia
    # < vyf      : lista ordenada de coordenadas x de cada vertice de cada fatia
    # < Δxf      : lista contendo a largura de cada fatia
	# ________________________________________

#determina os ângulo θ1 de xya e θ2 de xyb no círculo
	cosθ = (xya[0]-rx)/R
	sinθ = (xya[1]-ry)/R
	if ( sinθ >= 0.0 ):
		θ1 = acos(cosθ) * (180.0/pi)
	else:
		θ1 = ( 2.0*pi-acos(cosθ) ) * (180.0/pi)
	cosθ = (xyb[0]-rx)/R
	sinθ = (xyb[1]-ry)/R
	if ( sinθ >= 0.0 ):
		θ2 = acos(cosθ) * (180.0/pi)
	else:
		θ2 = ( 2.0*pi-acos(cosθ) ) * (180.0/pi)	
#ângulos entre 0.0 e 360.0 graus
	if ( θ2 <= θ1 ):
		Δθ = ( 360.0 - (θ1-θ2) ) / (nf)
	else :
		Δθ = ( θ2-θ1 ) / (nf)
#lista de valores a serem calculados
	θpf = [θ1]
	θf  = []
	vxf = []
	vyf = []
	Δxf = []
	for i in range(0,nf):
		vxf.append([])
		vyf.append([])
	lx = h/tan(α)
#divide a linha de ruptura em ângulos constantes entre xya e xyb
	for i in range(1,nf+1):
		θpf.append ( θpf[i-1] + Δθ )
#determina a geometria de cada fatia
	for i in range(0,nf):
		dtr = (pi/180.0)
		θf.append( 0.5 * ( θpf[i] + θpf[i+1] ) ) #ângulo central de cada fatia
		xpt1 = rx+R*cos( θpf[i]*dtr )
		vxf[i].append ( xpt1 )
		vyf[i].append ( ry+R*sin( θpf[i]*dtr ) )
		xpt2 = rx+R*cos( θpf[i+1]*dtr )
		vxf[i].append ( xpt2 )
		vyf[i].append ( ry+R*sin( θpf[i+1]*dtr ) )
		if ( xpt2 <= 0.0 ):
			vxf[i].append ( xpt2 )
			vyf[i].append ( 0.0 )
			vxf[i].append ( xpt1 )
			vyf[i].append ( 0.0 )
		elif ( xpt1 >= lx ):
			vxf[i].append ( xpt2 )
			vyf[i].append ( h )
			vxf[i].append ( xpt1 )
			vyf[i].append ( h )
		else:
			if ( xpt1 >= 0.0 and xpt2 <= lx ):
				vxf[i].append ( xpt2 )
				vyf[i].append ( tan(α)*xpt2 )
				vxf[i].append ( xpt1 )
				vyf[i].append ( tan(α)*xpt1 )
			elif ( xpt1 >= 0.0 and xpt2 >= lx ):			
				vxf[i].append ( xpt2 )
				vyf[i].append ( h )
				vxf[i].append ( lx )
				vyf[i].append ( h )			
				vxf[i].append ( xpt1 )
				vyf[i].append ( tan(α)*xpt1 )
			elif ( xpt1 <= 0.0 and xpt2 <= lx):
				vxf[i].append ( xpt2 )
				vyf[i].append ( tan(α)*xpt2 )
				vxf[i].append ( 0.0 )
				vyf[i].append ( 0.0 )										
				vxf[i].append ( xpt1 )
				vyf[i].append ( 0.0 )					
			elif ( xpt1 <= 0.0 and xpt2 >= lx):
				vxf[i].append ( xpt2 )
				vyf[i].append ( h )
				vxf[i].append ( lx )
				vyf[i].append ( h )
				vxf[i].append ( 0.0 )
				vyf[i].append ( 0.0 )												
				vxf[i].append ( xpt1 )
				vyf[i].append ( 0.0 )
#a largura de cada fatia é sempre a variação em x do ponto 1 para o 2				
		Δxf.append( xpt2-xpt1 )
#converte os ângulos com relação ao centro com 0.0 graus sendo o +270.0	
	for i in range(0,nf):
		if ( θf[i] <= 270.0 ):
			θf[i] = -(270.0 - θf[i])
		else:
			θf[i] = θf[i]-270.0
	return θf,vxf,vyf,Δxf

def geometria_fatias_p(xa,xb,yi,α,h,nf):
	
	# ________________________________________
	#
	# ALGORITMO PARA OBTER GEOMETRIA DE FATIAS
	# POLIGONAIS
	# ________________________________________
	#
	# > xa       : x da inteseção à esquerda
	# > xb       : x da inteseção à direita
	# > yi       : lista de coordenadas y dos pontos entre xa e xb
	# > α        : inclinação do talude
	# > h        : altura do talude
	# > nf       : número de fatias
   	# < θf       : lista com ângulo central de cada fatia medido em graus
   	# < vxf      : lista ordenada de coordenadas x de cada vertice de cada fatia
    # < vyf      : lista ordenada de coordenadas x de cada vertice de cada fatia
    # < Δxf      : lista contendo a largura de cada fatia
	# ________________________________________

#cria listas de coordenadas da poligonal
	vx = [xa]
	vy = [interpolacao_coordenadas_talude(xa,α,h)]
	lyi = len(yi)
	dx = (xb-xa)/(lyi+1.0)
	xi = xa
	for i in range(0,lyi):
		xi += dx
		vx.append(xi)
		vy.append(yi[i])
	vx.append(xb)
	vy.append(interpolacao_coordenadas_talude(xb,α,h))
#lista de valores a serem calculados
	θf  = []
	vxf = []
	vyf = []
	Δxf = []
	for i in range(0,nf):
		vxf.append([])
		vyf.append([])
	lx = h/tan(α)
#determina a geometria de cada fatia
	for i in range(0,nf):
		xpt1 = vx[i]
		vxf[i].append ( xpt1 )
		vyf[i].append ( vy[i] )
		xpt2 = vx[i+1]
		vxf[i].append ( xpt2 )
		vyf[i].append ( vy[i+1] )
#inclinação do ponto 1 para o 2				
		θ12 = atan( ( vy[i+1]-vy[i] ) / ( xpt2-xpt1 ) )
		θf.append( θ12*(180.0/pi) )		
		if ( xpt2 <= 0.0 ):
			vxf[i].append ( xpt2 )
			vyf[i].append ( 0.0 )
			vxf[i].append ( xpt1 )
			vyf[i].append ( 0.0 )
		elif ( xpt1 >= lx ):
			vxf[i].append ( xpt2 )
			vyf[i].append ( h )
			vxf[i].append ( xpt1 )
			vyf[i].append ( h )
		else:
			if ( xpt1 >= 0.0 and xpt2 <= lx ):
				vxf[i].append ( xpt2 )
				vyf[i].append ( tan(α)*xpt2 )
				vxf[i].append ( xpt1 )
				vyf[i].append ( tan(α)*xpt1 )
			elif ( xpt1 >= 0.0 and xpt2 >= lx ):			
				vxf[i].append ( xpt2 )
				vyf[i].append ( h )
				vxf[i].append ( lx )
				vyf[i].append ( h )			
				vxf[i].append ( xpt1 )
				vyf[i].append ( tan(α)*xpt1 )
			elif ( xpt1 <= 0.0 and xpt2 <= lx):
				vxf[i].append ( xpt2 )
				vyf[i].append ( tan(α)*xpt2 )
				vxf[i].append ( 0.0 )
				vyf[i].append ( 0.0 )										
				vxf[i].append ( xpt1 )
				vyf[i].append ( 0.0 )					
			elif ( xpt1 <= 0.0 and xpt2 >= lx):
				vxf[i].append ( xpt2 )
				vyf[i].append ( h )
				vxf[i].append ( lx )
				vyf[i].append ( h )
				vxf[i].append ( 0.0 )
				vyf[i].append ( 0.0 )												
				vxf[i].append ( xpt1 )
				vyf[i].append ( 0.0 )
#a largura de cada fatia é sempre a variação em x do ponto 1 para o 2				
		Δxf.append( xpt2-xpt1 )
	return θf,vxf,vyf,Δxf
	
def plota_talude_c(R,rx,ry,α,h,nf,titulo):
	
	# ________________________________________
	#
	# ALGORITMO PARA PLOTAR RESULTADOS (CIRCULAR)
	# ________________________________________
	#
	# > R        : raio do círculo de ruptura
	# > rx       : coordenada x do centro do círculo
	# > ry       : coordenada y do centro do círculo
	# > α        : inclinação do talude
	# > h        : altura do talude
	# > nf       : número de fatias
	# > titulo   : string para título do gráfico
	# ________________________________________

#determina as fatias e se existe a interseção
	xya,xyb,intercecao_valida = intersecao_circulo_com_talude(R,rx,ry,α,h)
	θf,vxf,vyf,Δxf = geometria_fatias_c(R,rx,ry,α,h,nf,xya,xyb)
	if ( intercecao_valida ):
#extremos do arco da circunferência
		cx = []
		cy = []
		cx.append(xya[0])
		cy.append(xya[1])
		cx.append(rx)
		cy.append(ry)
		cx.append(xyb[0])
		cy.append(xyb[1])
#medidas da geometria do problema	
		lx = h/tan(α)
		ly = h
		L  = sqrt(lx**2+ly**2)
		l1x,l1y = [-L,0.0] , [0.0,0.0]
		l2x,l2y = [0.0,lx] , [0.0,h]
		l3x,l3y = [lx,lx+L] , [h,h]
		soilx = [-L,lx+L,lx+L,lx,0.0,-L]
		soily = [-L,-L,h,h,0.0,0.0]
#plota em 600x600 pixels
		plt.figure( 1,figsize=(6,6) )
#contorno do solo		
		plt.plot(l1x,l1y,c='k',lw=1.0)
		plt.plot(l2x,l2y,c='k',lw=1.0)
		plt.plot(l3x,l3y,c='k',lw=1.0)
#círculo
		plt.plot(cx,cy,'--',c='r',lw=1.0)
		plt.plot(xya[0],xya[1],'o',c='r',ms=3.0)
		plt.plot(xyb[0],xyb[1],'o',c='r',ms=3.0)
		plt.plot(rx,ry,'+',c='r',ms=8.0)
#polígono representando o solo		
		plt.fill(soilx,soily,c='g',alpha=0.5)
#fatias
		for i in range(0,nf):
			plt.fill(vxf[i],vyf[i],c='m',alpha=0.3,lw=0.5)
#ajusta a plotagem e mostra na tela
		plt.xlim(-L, lx+L)
		plt.ylim(-L, lx+L)
		plt.title(titulo)
		info = '( '+str(int(round(rx,0)))+' , '+str(int(round(ry,0)))+' )'+'  R = '+str(int(round(R,0)))
		plt.annotate(text=info,xy=(0.5*rx,1.1*ry),fontsize=7.0,c='r')
		plt.show()
	else:
		print('Plotagem Inválida!')

def plota_talude_p(vxf,vyf,α,h,nf,titulo):
	
	# ________________________________________
	#
	# ALGORITMO PARA PLOTAR RESULTADOS (POLIGONAL)
	# ________________________________________
	#
  	# > vxf      : lista ordenada de coordenadas x de cada vertice de cada fatia
    # > vyf      : lista ordenada de coordenadas x de cada vertice de cada fatia
	# > α        : inclinação do talude
	# > h        : altura do talude
	# > nf       : número de fatias
	# > titulo   : string para título do gráfico
	# ________________________________________

#medidas da geometria do problema	
	lx = h/tan(α)
	ly = h
	L  = sqrt(lx**2+ly**2)
	l1x,l1y = [-L,0.0] , [0.0,0.0]
	l2x,l2y = [0.0,lx] , [0.0,h]
	l3x,l3y = [lx,lx+L] , [h,h]
	soilx = [-L,lx+L,lx+L,lx,0.0,-L]
	soily = [-L,-L,h,h,0.0,0.0]
#plota em 600x600 pixels
	plt.figure( 1,figsize=(6,6) )
#contorno do solo		
	plt.plot(l1x,l1y,c='k',lw=1.0)
	plt.plot(l2x,l2y,c='k',lw=1.0)
	plt.plot(l3x,l3y,c='k',lw=1.0)
#polígono representando o solo		
	plt.fill(soilx,soily,c='g',alpha=0.5)
#fatias
	for i in range(0,nf):
		plt.fill(vxf[i],vyf[i],c='m',alpha=0.3,lw=0.5)
#ajusta a plotagem e mostra na tela
	plt.xlim(-L, lx+L)
	plt.ylim(-L, lx+L)
	plt.title(titulo)
	plt.show()
	
def area_poligono(xp,yp):
	
	# ________________________________________
	#
	# ALGORITMO PARA CALCULAR ÁREA DE POLÍGONO
	# ________________________________________
	#
	# > xp  : lista de coordenadas x ordenadas do polígono
	# > yp  : lista de coordenadas y ordenadas do polígono
	# < area: area de um polígono simples
	# ________________________________________	
	
#número de pontos
	n = len(xp)
	area = 0.0
#calcula a área pela fórmula de Shoelace	
	j = n-1
	for i in range(0,n):
		area += (xp[j] + xp[i]) * (yp[j] - yp[i])
		j = i
	area = 0.5*abs(area) 
	return area

def calcula_FS(x):

	# ________________________________________
	#
	# ALGORITMO PARA CALCULAR O FATOR DE 
	# SEGURANÇA PELO MÉTODO DE FELLENIUS
	# PARA UM DADO CÍRCULO DE RUPTURA
	# ________________________________________
	#
	# > x   : lista contendo [R, rx, ry]
	# < FS  : fator de segurança
	# ________________________________________	

#dados globais do problema	
	global c,ϕ,γ,α,h,nf
	global R,rx,ry
	global metodo
	global xa,xb,yi
	
	if ( metodo == 'Fellenius' or metodo == 'Bishop' ):
#dados do círculo de ruptura
		R  =  x[0]
		rx =  x[1]
		ry =  x[2]
#verifica a interseção
		xya,xyb,intercecao_valida = intersecao_circulo_com_talude(R,rx,ry,α,h)
#o círculo cria uma interseção válida com o talude
		if ( intercecao_valida ):
			if ( metodo == 'Fellenius' ):
				θf,vxf,vyf,Δxf = geometria_fatias_c(R,rx,ry,α,h,nf,xya,xyb)
				F1 = 0.0
				F2 = 0.0
				F3 = 0.0
				for i in range(0,nf):
					A = area_poligono(vxf[i],vyf[i])
					W = A*γ
					F1 += Δxf[i] / cos( θf[i]*(pi/180.0) )
					F2 += W*cos(θf[i]*(pi/180.0))
					F3 += W*sin(θf[i]*(pi/180.0))
				F1 *= c
				F2 *= tan(ϕ)
				FS = ( F1 + F2 ) / F3
			elif ( metodo == 'Bishop' ):
				FS = sp.newton(fs,2.0,maxiter=100,tol=1E-6)
			return FS
		else:
#o círculo não cria uma interseção válida com o talude
			FS = 9999.9
			return FS
	
	elif ( metodo == 'Janbu' ):
#dados do polígono de ruptura
		xa =  x[0]
		yi = []
		for i in range(1,nf):
			yi.append( x[i] )
		xb =  x[nf]
#verifica a interseção
		ya,yb,intercecao_valida = intersecao_poligonal_com_talude(xa,xb,yi,α,h)
#a poligonal cria uma interseção válida com o talude
		if ( intercecao_valida ):
			FS = sp.newton(fj,2.0,maxiter=100,tol=1E-6)
			return FS
		else:
#o círculo não cria uma interseção válida com o talude
			FS = 9999.9
			return FS
		
#_______________________________________________________________________

c = 20.0                   #[kN/m²] coesão do solo
ϕ = 30.0 * (pi/180.0)      #[rad]   ângulo de atrito do solo
γ = 18.5                   #[kN/m³] peso específico do solo
α = 45.0 * (pi/180.0)      #[rad]   inclinação do talude
h = 15.0                   #[m]     altura do talude

nf = 50                    #número de fatias
metodo = 'Fellenius'           #Método utilizado: Fellenius / Bishop / Janbu / Log Spiral
otimizacao = 'Gradientes'  #Método de otimização Gradientes / Global
plotagem0 = True           #True plota solução inicial / False apenas Fmin

#_______________________________________________________________________


if ( metodo == 'Fellenius' or metodo == 'Bishop' ):
#solução inicial FS_0
	#R  = sqrt((h*cos(α))**2+h**2)
	rx =  (1.0/2.0) * h/tan(α)
	ry =  (4.0/3.0) * h
	R  =  (2.0/2.0) * sqrt(rx**2+ry**2)
	x0 = [ R, rx, ry ]
	FS = calcula_FS(x0)
	if ( plotagem0 ):
		plota_talude_c(R,rx,ry,α,h,nf,'FS Inicial = '+str(round(FS,2))+' ('+metodo+')' )
#solução mínima FS_MIN
	if ( otimizacao == 'Gradientes' ):
		ineq_cons = {'type': 'ineq',
				      'fun': lambda x: np.array([ ( abs(tan(α)*x[1]-x[2])/( sqrt((tan(α))**2+1.0) ) ) ,
				                                  ( x[0]-1.0 )            ] ) }
		
		res = sp.minimize(f, x0, method='SLSQP', jac={'3-point'},  constraints=[ineq_cons], 
						options={'ftol': 1E-12, 'disp': True, 'maxiter':1000,'iprint': 2, 'eps': 1E-8})
		#res = sp.minimize(f, x0, method='BFGS', jac={'3-point'},
		#				options={'gtol': 1E-12, 'disp': False, 'maxiter':1000, 'eps': 1E-8})
	elif ( otimizacao == 'Global' ):
		res = sp.basinhopping(f, x0, niter=50, disp=True)
	xmin  = res.x
	FSMIN = res.fun
#plotagem Fmin
	plota_talude_c(xmin[0],xmin[1],xmin[2],α,h,nf,'FS Mínimo = '+str(round(FSMIN,2))+' ('+metodo+')' )

if ( metodo == 'Janbu' ):
#solução inicial FS_0
	rx =  (1.0/2.0) * h/tan(α)
	ry =  (4.0/3.0) * h
	R  =  (2.0/2.0) * sqrt(rx**2+ry**2)
	xya,xyb,intercecao_valida = intersecao_circulo_com_talude(R,rx,ry,α,h)
#poligonal inicial é um círculo
	xa = xya[0]
	xb = xyb[0]
	dx =  (xb-xa)/(nf+1.0)
	xi = xa
	x0 = [xa]
	for i in range(1,nf):
		xi += dx
		x0.append( ry-sqrt(R**2-(xi-rx)**2)  )
	x0.append(xb)
	FS = calcula_FS(x0)
	if ( plotagem0 ):
		ya,yb,intercecao_valida = intersecao_poligonal_com_talude(xa,xb,x0[1:nf],α,h)
		θf,vxf,vyf,Δxf = geometria_fatias_p(xa,xb,x0[1:nf],α,h,nf)
		plota_talude_p(vxf,vyf,α,h,nf,'FS Inicial = '+str(round(FS,2))+' ('+metodo+')')
#solução mínima FS_MIN
	if ( otimizacao == 'Gradientes' ):
		#ineq_cons = {'type': 'ineq',
		#		      'fun': lambda x: np.array([ ( abs(tan(α)*x[1]-x[2])/( sqrt((tan(α))**2+1.0) ) ) ,
		#		                                  ( x[0]-1.0 )            ] ) }
		res = sp.minimize(f, x0, method='SLSQP', jac={'3-point'},#  constraints=[ineq_cons], 
						options={'ftol': 1E-12, 'disp': True, 'maxiter':1000,'iprint': 2, 'eps': 1E-8})		
#		res = sp.minimize(f, x0, method='SLSQP', jac={'3-point'},
#						options={'ftol': 1E-12, 'disp': True, 'maxiter':500,'iprint': 2, 'eps': 1E-8})
	elif ( otimizacao == 'Global' ):
		res = sp.basinhopping(f, x0, niter=50, disp=True )
	xmin  = res.x
	FSMIN = res.fun
#plotagem Fmin
	ya,yb,intercecao_valida = intersecao_poligonal_com_talude(xmin[0],xmin[nf],xmin[1:nf],α,h)
	θf,vxf,vyf,Δxf = geometria_fatias_p(xmin[0],xmin[nf],xmin[1:nf],α,h,nf)	
	plota_talude_p(vxf,vyf,α,h,nf,'FS Min = '+str(round(FSMIN,2))+' ('+metodo+')')
						
#_______________________________________________________________________


