import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class Design_project():
	
	def __init__(self):
		T = sp.Symbol('T')
		R = 8.314			#J/mol K
		t = T/1000
		
		#site density
		omega_1_FR = 35			#mol/m^3
		omega_2_FR = 70
		omega_1_AG = 22.5
		omega_2_AG = 60
		
		#reference 1, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0)
		self.H_O2 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
		self.S_O2 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 2, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0)
		self.H_H2 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
		self.S_H2 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 3, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224)
		self.H_CO2 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H + (-393.51)
		self.S_CO2 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 4, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (25.56759, 6.096130, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665, -110.5271)
		self.H_CO = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H + (-110.53)
		self.S_CO = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 5, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967, -241.8264)
		self.H_H2O = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H + (-241.826)
		self.S_H2O = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		#reference 6, H in kJ/mol, S in J/mol*K
		A, B, C, D, E, F, G, H = (-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, 158.7163, -74.8731)
		self.H_CH4 = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H + (-74.6)
		self.S_CH4 = A*sp.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
		
		
		self.delG_O2 = self.H_O2 - T*self.S_O2/1000			#kJ/mol
		self.delG_H2 = self.H_H2 - T*self.S_H2/1000
		self.delG_CO2 = self.H_CO2 - T*self.S_CO2/1000				
		self.delG_CO = self.H_CO - T*self.S_CO/1000
		self.delG_H2O = self.H_H2O - T*self.S_H2O/1000
		self.delG_CH4 = self.H_CH4 - T*self.S_CH4/1000
		
		#For Ce2O3 + 0.5O2 -> 2CeO2
		self.delG_S1_FR = -247.8 - T*(29.4)/1000				#kJ/mol O
		self.delG_S2_FR = -252.8 - T*(29.4)/1000
		self.delG_S1_AG = -252.8 - T*(-0.6)/1000
		self.delG_S2_AG = -243.8 - T*(-0.6)/1000
		
		#For Ce2O3 + 0.5O2 <-> 2CeO2
		self.k1f_FR = 4.47*10**4 * sp.exp(-25.0*1000/(R*T))		#m^3/mol
		self.k4f_FR = 1.78*10**6 * sp.exp(-40.0*1000/(R*T))
		self.k1f_AG = 1.47*10**4 * sp.exp(-25.0*1000/(R*T))
		self.k4f_AG = 1.78*10**6 * sp.exp(-40.0*1000/(R*T))
		self.R_G_1 = self.delG_S1_FR
		
		#for H2 + 2CeO2 <-> Ce2O3 + H2O
		self.k2f_FR = 3.16*10**6 * sp.exp(-45.0*1000/(R*T))
		self.k5f_FR = 1.58*10**7 * sp.exp(-62.5*1000/(R*T))
		self.k2f_AG = 3.16*10**4 * sp.exp(-45.0*1000/(R*T))
		self.k5f_AG = 1.00*10**6 * sp.exp(-62.5*1000/(R*T))
		self.R_G_2 = self.delG_H2O - (self.delG_S1_FR + 0.5*self.delG_O2) - self.delG_H2
		
		#for CO + 2CeO2 <-> Ce2O3 + CO2
		self.k3f_FR = 1.58*10**7 * sp.exp(-58.0*1000/(R*T))
		self.k6f_FR = 4.68*10**6 * sp.exp(-80.0*1000/(R*T))
		self.k3f_AG = 2.00*10**6 * sp.exp(-58.0*1000/(R*T))
		self.k6f_AG = 5.62*10**5 * sp.exp(-80.0*1000/(R*T))
		self.R_G_3 = self.delG_CO2 - self.delG_CO - self.delG_S1_FR
		
		#for CO + H2O <-> CO2 + H2
		self.k_WGS_FR = 2.82*10**10 * sp.exp(-80*1000/(R*T))
		self.k_WGS_AG = 1.82*10**9 * sp.exp(-80*1000/(R*T))
		self.R_G_WGS = self.delG_H2 + self.delG_CO2 - self.delG_H2O - self.delG_CO
		
		#for CH4 + H2O <-> CO + 3H2
		self.k_SR_FR = None
		self.k_SR_AG = 1.15*10**6 * sp.exp(-85.2*1000/(R*T))
		self.R_G_SR = 3*self.delG_H2 + self.delG_CO - self.delG_H2O - self.delG_CH4
		
	def rate_const_comp(self):
		T = sp.Symbol('T')
		R = 8.314 		#J/mol K
		T_array = np.arange(100, 600+0.1, 0.1)+273
		k1_array = np.zeros_like(T_array)
		k2_array = np.zeros_like(T_array)
		k3_array = np.zeros_like(T_array)
		k4_array = np.zeros_like(T_array)
		
		k1 = self.k1f_FR
		RG1 = self.R_G_1
		K1 = sp.exp(-RG1*1000/(R*T))
		
		k2 = self.k2f_FR
		RG2 = self.R_G_2
		K2 = sp.exp(-RG2*1000/(R*T))
		
		k3 = self.k3f_FR
		RG3 = self.R_G_3
		K3 = sp.exp(-RG3*1000/(R*T))
		
		k4 = self.k_WGS_FR
		RG4 = self.R_G_WGS
		K4 = sp.exp(-RG4*1000/(R*T))
		
		for i in range(len(T_array)):
			k1_array[i] = sp.log(k1.subs(T, T_array[i]))
			k2_array[i] = sp.log(k2.subs(T, T_array[i])/K2.subs(T, T_array[i]))
			k3_array[i] = sp.log(k3.subs(T, T_array[i])/K3.subs(T, T_array[i]))
			k4_array[i] = sp.log(k4.subs(T, T_array[i]))
		
		plt.close()
		plt.plot(T_array, k1_array, label = 'by oxygen')
		plt.plot(T_array, k2_array, label = 'by water')
		plt.plot(T_array, k3_array, label = 'by carbon dioxide')
		plt.plot(T_array, k4_array, label = 'WGS')
		plt.grid()
		plt.legend()
		plt.xlabel('Temperature in K')
		plt.ylabel('ln(k) in m^3/(mol*s)')
		plt.title('Rate constant comparison')
		plt.show()
		
if __name__ == '__main__':
	test = Design_project()
	test.rate_const_comp()
	
	
'''
reference 1 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1&Type=JANAFG&Table=on
reference 2 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1&Type=JANAFG&Plot=on
reference 3 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1
reference 4 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Units=SI&Mask=7
reference 5 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=1#Refs
reference 6 : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Units=SI&Mask=1
'''
