import math as m
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def snell_fresnel_somatorio(n1, n2, a_inc, n_interfaces, c_onda):
    
    c_onda = c_onda*10**(-9)
    v_luz = 3*10**8
    frequencia = v_luz/c_onda
    t = 1
    omega = 2*m.pi*frequencia
    d_interface = 10**(-6)
    E0 = 1
    
    ## n1 para n2
    seno_inc1 = m.sin(a_inc*m.pi/180)
    cos_inc1 = m.cos(a_inc*m.pi/180)
    a_refracao1 = m.asin((seno_inc1*n1)/n2)*(180/m.pi)
    seno_refracao1 = m.sin(a_refracao1*m.pi/180)
    cos_refracao1 = m.cos(a_refracao1*m.pi/180)
    caminho1 = d_interface/cos_refracao1
    indice1 = n1
    c_onda_ajus1 = c_onda/indice1
    fase1 = abs(c_onda_ajus1 - caminho1)*2*m.pi/c_onda_ajus1
    
    r_paralelo1 = ((n1*cos_refracao1-n2*cos_inc1)/(n1*cos_refracao1+n2*cos_inc1))**2
    t_paralelo1 = ((2*n1*cos_inc1)/(n1*cos_refracao1+n2*cos_inc1))**2
    r_paralelo_porc1 = r_paralelo1/(r_paralelo1+t_paralelo1)
    t_paralelo_porc1 = t_paralelo1/(r_paralelo1+t_paralelo1)
    
    r_perpendicular1 = (n1*cos_inc1-n2*cos_refracao1)/(n1*cos_inc1+n2*cos_refracao1)
    t_perpendicular1 = 2*n1*cos_inc1/(n1*cos_inc1+n2*cos_refracao1)
    r_perpendicular_porc1 = r_perpendicular1/(r_perpendicular1+t_perpendicular1)
    t_perpendicular_porc1 = t_perpendicular1/(r_perpendicular1+t_perpendicular1)
    
    ## n2 para n1
    
    seno_inc2 = m.sin(a_refracao1*m.pi/180)
    cos_inc2 = m.cos(a_refracao1*m.pi/180)
    a_refracao2 = m.asin((seno_inc2*n2)/n1)*(180/m.pi)
    seno_refracao2 = m.sin(a_refracao2*m.pi/180)
    cos_refracao2 = m.cos(a_refracao2*m.pi/180)
    caminho2 = d_interface/cos_refracao2
    indice2 = n2
    c_onda_ajus2 = c_onda/indice2
    fase2 = abs(c_onda_ajus2 - caminho2)*2*m.pi/c_onda_ajus2
    
    r_paralelo2 = ((n2*cos_refracao2-n1*cos_refracao1)/(n2*cos_refracao2+n1*cos_refracao1))**2
    t_paralelo2 = ((2*n2*cos_refracao1)/(n2*cos_refracao2+n1*cos_refracao1))**2
    r_paralelo_porc2 = r_paralelo2/(r_paralelo2+t_paralelo2)
    t_paralelo_porc2 = t_paralelo2/(r_paralelo2+t_paralelo2)
    
    r_perpendicular2 = ((n2*cos_refracao1-n1*cos_refracao2)/(n2*cos_refracao1+n1*cos_refracao2))**2
    t_perpendicular2 = (2*n2*cos_refracao1/(n2*cos_refracao1+n1*cos_refracao2))**2
    r_perpendicular_porc2 = r_perpendicular2/(r_perpendicular2+t_perpendicular2)
    t_perpendicular_porc2 = t_perpendicular2/(r_perpendicular2+t_perpendicular2)
    
    i = 0
    SomaE = 0
    
    for i in range(n_interfaces):
        
        #if (i%2) == 0: 
        if n2>n1:
            fases = ((i/2)*fase1 + (i/2)*fase2)-m.pi
            faseF = fases + (i/2)*fase1 + (i/2)*fase2
            E =1
            E = E0*(t_perpendicular_porc1*(i-1))#(t_perpendicular_porc2*(i-1))(r_perpendicular_porc2)  
        else:
            fases = (((i+1)/2)*fase1 + ((i-1)/2)*fase2)-m.pi
            faseF = fases + ((i+1)/2)*fase1 + ((i-1)/2)*fase2
            E = E0*(t_perpendicular_porc1*(i-1))#(t_perpendicular_porc2*(i-1))(r_perpendicular_porc2)
        
        n1, n2 = n2, n1
        
        E_completo = E**2*m.cos((2*m.pi*0 + faseF)*m.pi/180)
        SomaE = SomaE + E_completo
    
    return SomaE**2
data = []
i = 0
for i in range(1000):
    resultado = snell_fresnel_somatorio(1.45, 1.4505, 30, 5000, 1000 + i)
    data.append(resultado)

data = pd.DataFrame(data, columns = ['Intensidade'])
data['comprimento de onda (nm)'] = 0
    
for j in range(len(data)):
    data['comprimento de onda (nm)'].loc[j] = j + 1500
data.plot(x = 'comprimento de onda (nm)', color = 'red', figsize = (10,6))
plt.ylabel('Intensidade', size = 14)
plt.yticks(size = 14)
plt.show()