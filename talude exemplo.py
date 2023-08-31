#_______________________________________________________________________
#
#                      [ENG01060 - AVALIAÇÃO 2]
#
#    Programação de Métodos Numéricos Aplicados à Engenharia Civil
#
#       Universidade Federal do Rio Grande do Sul >> [UFRGS]
#       Departamento de Engenharia Civil          >> [DECIV]
#
#                                       Professor: Renato Vaz Linn
#
#  CARTÃO                           	 NOME                    
# [000001]   [                   Nome Sobrenome                   ] 
# [000002]   [                   Nome Sobrenome                   ] 
# [000003]   [                   Nome Sobrenome                   ]  
# [000004]   [                   Nome Sobrenome                   ]  
# [000005]   [                   Nome Sobrenome                   ] 
#
#                                                                              
#                           *///////////////////////////////////            
#                      .////////////////////////////////*///////            
#                  *////////////////////////////////*///////////            
#              *////////////////////////////////////////////////            
#         .////////////////////////////////*////////////////////            
#     *////////////////////////////////*////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************//////////////////////////            
#    *********************************///////////////////////*              
#    *********************************///////////////////.                  
#    *********************************///////////////                       
#    *********************************//////////*                           
#    *********************************//////.                               
#    *********************************//                                    
#                                            
#                            AVALIAÇÃO:
#
# [SOLUÇÃO INICIAL] - A solução inicial de FS para x0 está correta e
#                     bem definida/apresentada (1 pt)
#  [CÁLCULO DE FS]  - O cálculo de FS (Fellenius ou Bishop) está
#                     corretamente calculado e implementado (3 pts)
# [CÁLCULO DE FSMIN]- O programa consegue encontar o valor de FS min
#                     através de um algoritmo numérico de otimização
#                     (4 pts)
#[VERIFICAÇÃO NORMA]- O código veririca a adequação do talude aos 
#                     critérios de estabilidade mínimos exigidos
#                     pela norma (1 pt)
#     [CLAREZA]     - O código é claro de entender, apresentando
#                     todas informações para sua compreensão de forma
#                     adequada (1 pt)
#_______________________________________________________________________  

import matplotlib.pyplot as plt
from math import *

#função objetivo para ser invocada por algoritmo de otimização
def f(x):
	f = calcula_FS(x)
	return f

def intersecao_circulo_com_talude(R,rx,ry,α,h):
	
	# ________________________________________
	#
	# ALGORITMO DE INTERSEÇÃO
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
	else:
		intercecao_valida = False
	return xya,xyb,intercecao_valida

def geometria_fatias(R,rx,ry,α,h,nf,xya,xyb):
	
	# ________________________________________
	#
	# ALGORITMO PARA OBTER GEOMETRIA DE FATIAS
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

def plota_talude(R,rx,ry,α,h,nf,titulo):
	
	# ________________________________________
	#
	# ALGORITMO PARA PLOTAR RESULTADOS
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
	θf,vxf,vyf,Δxf = geometria_fatias(R,rx,ry,α,h,nf,xya,xyb)
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
	# SEGURANÇA PARA UM DADO CÍRCULO DE RUPTURA
	# ________________________________________
	#
	# > x   : lista contendo [R, rx, ry]
	# < FS  : fator de segurança
	# ________________________________________	

#dados globais do problema	
	global c,ϕ,γ,α,h,nf
	global R,rx,ry
#dados do círculo de ruptura
	R  =  x[0]
	rx =  x[1]
	ry =  x[2]
#verifica a interseção
	xya,xyb,intercecao_valida = intersecao_circulo_com_talude(R,rx,ry,α,h)
#o círculo cria uma interseção válida com o talude
	if ( intercecao_valida ):
		θf,vxf,vyf,Δxf = geometria_fatias(R,rx,ry,α,h,nf,xya,xyb)
		#exemplo de como calcular a área de cada fatia:
		#for i in range(0,nf):
#			A = area_poligono(vxf[i],vyf[i])
		#----aqui você deve calcular FS e retornar o valor como saída----
		#return ( FS )
	else:
#o círculo não cria uma interseção válida com o talude
		FS = 9999.9
		return FS
				
#_______________________________________________________________________
#						   
# CÓDIGO PRINCIPAL
#
c = 20.0                   #[kN/m²] coesão do solo
ϕ = 25.0 * (pi/180.0)      #[rad]   ângulo de atrito do solo
γ = 18.5                   #[kN/m³] peso específico do solo
α = 30.0 * (pi/180.0)      #[rad]   inclinação do talude
h = 15.0                   #[m]     altura do talude

rx =  12.99                #[m] x do centro do círculo inicial
ry =  25.00                #[m] y do centro do círculo inicial
R  =  42.2                 #[m] Raio círculo inicial

nf = 5                     #número de fatias

#solução inicial FS_0
x0 = [ R, rx, ry ]
FS = calcula_FS(x0)
plota_talude(R,rx,ry,α,h,nf,'Título do Gráfico' )

#encontrar  x* = [R*, rx*, ry*] que fornece FSmin = min(FS)
#plotar
#verificar critérios da norma comparando-a com FSmin


