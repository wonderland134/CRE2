import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class rate_const():
	
	def __init__(self):
		T = sp.Symbol('T')
		R = 8.314			#J/mol K
		
		
		
		t = T/1000
		#reference 1, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0)
		self.H_O2 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
		self.S_O2 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 2, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0)
		self.H_H2 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
		self.S_H2 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		self.delG_O2 = self.H_O2 - T*self.S_O2/1000			#kJ/mol
		self.delG_H2 = self.H_H2 - T*self.S_H2/1000
		self.delG_CO = -0.0888*T - 111.0465				
		self.delG_CO2 = -0.0016*T - 394.0838
		self.delG_H2O = 0.0526*T - 245.0239
		
		#For Ce2O3 + 0.5O2 -> 2CeO2
		self.delG_S1_FR = -247.8 - T*(29.4)/1000				#kJ/mol O
		
		#For Ce2O3 + 0.5O2 <-> 2CeO2
		self.k1f_FR = 4.47*10**4 * sp.exp(-25.0*1000/(R*T))		#m^3/mol
		self.R_G_1 = self.delG_S1_FR
		
		#for H2+2CeO2 <-> Ce2O3 + H2O
		self.k2f_FR = 3.16*10**6 * sp.exp(-45.0*1000/(R*T))
		self.R_G_2 = self.delG_H2O - (self.delG_S1_FR + 0.5*self.delG_O2) - self.delG_H2
		
		#for CO+2CeO2 <-> Ce2O3 + CO2
		self.k3f_FR = 1.58*10**7 * sp.exp(-58.0*1000/(R*T))
		self.R_G_3 = self.delG_CO2 - self.delG_CO - self.delG_S1_FR
	
	def calc(self):
		T = sp.Symbol('T')
		R = 8.314 		#J/mol K
		T_array = np.arange(100, 600+1, 1)
		k1_array = np.zeros_like(T_array)
		k2_array = np.zeros_like(T_array)
		k3_array = np.zeros_like(T_array)
		
		k1 = self.k1f_FR
		RG1 = self.R_G_1
		K1 = sp.exp(-RG1*1000/(R*T))
		
		k2 = self.k2f_FR
		RG2 = self.R_G_2
		K2 = sp.exp(-RG2*1000/(R*T))
		
		k3 = self.k3f_FR
		RG3 = self.R_G_3
		K3 = sp.exp(-RG3*1000/(R*T))
		
		for i in range(len(T_array)):
			k1_array[i] = k1.subs(T, T_array[i])
			k2_array[i] = k2.subs(T, T_array[i])/K2.subs(T, T_array[i])
			k3_array[i] = k3.subs(T, T_array[i])/K3.subs(T, T_array[i])
		
		plt.close()
		plt.plot(T_array, k1_array, label = 'R1 rate const')
		plt.plot(T_array, k2_array, label = 'R2 rate const')
		plt.plot(T_array, k3_array, label = 'R3 rate const')
		plt.grid()
		plt.legend()
		plt.xlabel('Temperature in K')
		plt.ylabel('rate constant in m^3/(mol*s)')
		plt.title('rate constant comparison')
		plt.show()
		
if __name__ == '__main__':
	test = rate_const()
	test.calc()
