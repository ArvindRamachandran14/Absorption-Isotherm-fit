import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import Model, ODR, Data, RealData
import Isotherm_data
import math
from scipy.optimize import curve_fit
import json

class Isotherm_Fit():

	def __init__(self):

		self.Temperatures = np.asarray(Isotherm_data.Temperatures) # Experimental Temperatures, np array
		self.Pressures_dict = Isotherm_data.pressures # Experimental Pressures^0.5, dictionary
		self.H_conc_dict = Isotherm_data.H_conc # Experimental Concentrations, dictionary

		self.k_B = 8.61733034E-5 #ev/K 
		self.P_0 = 760 #Torr

		self.H_M_ratio_dict = {} # Create H_M_ratio_dict
		self.Pressures_fit = {}

		self.H_M_ratio_data = []
		self.Pressures_data = []
		self.Temperatures_data = []
		
		self.read_input_file()
		self.intialize()


	def read_input_file(self):

		self.f = open('Input_file.json')
		self.data = json.load(self.f)["Fitting parameters"][0]  	


	def intialize(self):

		if self.data['Fit'] == "Partial":
			self.Train_indices = self.data["Partial_fit_isotherm_indices"]

		elif self.data['Fit'] == "Full":
			self.Train_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

		self.H_conc_limit = float(self.data["Fit_limit"])

		for i in range(len(self.H_conc_dict)):
			self.H_M_ratio_dict[i] = []
			for j in range(len(self.H_conc_dict[i])):
				y = self.H_conc_dict[i][j]
				if y <= self.H_conc_limit:
					x = y/(1-y)
					self.H_M_ratio_dict[i].append(x)

		for i in self.Train_indices:
			for j in range(len(self.H_M_ratio_dict[i])):
				self.H_M_ratio_data.append(self.H_M_ratio_dict[i][j])
				self.Pressures_data.append(self.Pressures_dict[i][j])
				self.Temperatures_data.append(self.Temperatures[i])

		self.H_M_ratio_data = np.asarray(self.H_M_ratio_data)
		self.Pressures_data = np.asarray(self.Pressures_data)
		self.Temperatures_data = np.asarray(self.Temperatures_data)

		if self.data["A"] == "Fixed":
			self.A=self.data['A_value']

		if self.data["Temp_dependence"] == "False":
			self.B=0

		elif self.data["B"] == "Fixed":
			self.B=self.data['B_value']

	def plot_raw_data(self):

		fig = plt.figure()
		ax = fig.add_subplot(111)
		for i in range(len(self.Pressures_dict)):
			ax.plot(self.H_M_ratio_dict[i], self.Pressures_dict[i][:len(self.H_M_ratio_dict[i])], 'ko')

		ax.set_xlabel(r'Hydrogen to metal ratio $(\mathrm{\frac{H}{M}})$', fontsize=16)
		ax.set_ylabel('$\mathregular{\sqrt{Pressure (Torr)}}$', fontsize=16)
		plt.tight_layout()
		plt.show()


	def ODR_fit_func(self, parm, var):

		order = self.data["Enthalpy_order"]
		H = 0
		index=0

		while index < order:
			H+=(index+1)*parm[index]*(var[0]**index)
			index+=1

		if self.data["Temp_dependence"] == "True":
			R = parm[4] + parm[5]*var[1] + parm[6]*(var[1]**2)
		else:
			R = parm[4] 

		P = (1/(self.k_B*var[1]))*H + np.log(var[0]/(R-var[0])) + R/(R-var[0])


		if self.data["A"] == "Fixed":
			P+=self.A
		else:
			P+=parm[7]

		if self.data["Temp_dependence"] == "True":
			if self.data["B"] == "Fixed":
				P+=self.B*var[1]
			else:
				P+=parm[8]*var[1]

		'''
		while order>0:
			H+=order*parm[order-1]*(var[0]**(order-1))
			order-=1
		

		i = self.data["Enthalpy_order"]

		if self.data["Temp_dependence"] == "True":
			R = parm[i] + parm[i+1]*var[1] + parm[i+2]*(var[1]**2)
			i+=3
		else:
			R = parm[i] 
			i+=1
	

		P = (1/(self.k_B*var[1]))*H + np.log(var[0]/(R-var[0])) + R/(R-var[0])

		if self.data["A"] == "Fixed":
			P+=self.A
		else:
			P+=parm[i]
			i+=1

		if self.data["Temp_dependence"] == "True":
			if self.data["B"] == "Fixed":
				P+=self.B*var[1]
			else:
				P+=parm[i]*var[1]
	
		'''
		P = (self.P_0**0.5)*np.exp(P)
		

		return P 	

	def ODR_fit(self):

		parameter_initialization = [0, 0, 0, 0, 1, 0, 0, 0, 0]

		'''
		parameter_initialization = []
	
		for i in range(self.data["Enthalpy_order"]):
			parameter_initialization.append(0)
		
		parameter_initialization.append(1)

		if self.data["Temp_dependence"] == "True":
			parameter_initialization.append(0)
			parameter_initialization.append(0)

		if self.data["A"] == "Free":
			parameter_initialization.append(0)
		if self.data["B"] == "Free": 
			parameter_initialization.append(0)
		'''
		self.odr_model = Model(self.ODR_fit_func)
		self.mydata = Data([self.H_M_ratio_data, self.Temperatures_data], self.Pressures_data)
		self.myodr = ODR(self.mydata, self.odr_model, beta0=np.asarray(parameter_initialization), maxit=10000000)
		self.myoutput = self.myodr.run()

		print(self.myoutput.beta)

		self.assign_fitting_constants()

		order = self.data["Enthalpy_order"]

		for i in range(len(self.H_M_ratio_dict)):

			self.Pressures_fit[i] = []
			T = self.Temperatures[i]
			x = np.asarray(self.H_M_ratio_dict[i])

			if self.data["Temp_dependence"] == "True":
				R = self.myoutput.beta[order] + self.myoutput.beta[order+1]*T + self.myoutput.beta[order+2]*(T**2)	
			else:
				R = self.myoutput.beta[order]

			H = self.E + 2*self.Alpha*x + 3*self.Beta*(x**2) + 4*self.Gamma*(x**3)

			self.Pressures_fit[i] = self.P_0**0.5*np.exp(self.A + self.B*T + (1/(self.k_B*T))*H + np.log(x/(R-x)) + R/(R-x))

		self.plot_fit()

	def assign_fitting_constants(self):

		index = 0
		self.E = self.myoutput.beta[0]
		index+=1

		order = self.data["Enthalpy_order"]

		if order > 1:
			self.Alpha = self.myoutput.beta[index]
			index+=1
			if order > 2:
				self.Beta = self.myoutput.beta[index]
				index+=1
				if order > 3:
					self.Gamma = self.myoutput.beta[index]
					index+=1
				else:
					self.Gamma = 0
			else:
				self.Beta=0
				self.Gamma=0
		else:
			self.Alpha =0
			self.Beta=0
			self.Gamma=0
		
		self.R_0 = self.myoutput.beta[index]
		index+=1

		if self.data["Temp_dependence"] == "True":
			self.R_1 = self.myoutput.beta[index]
			self.R_2 = self.myoutput.beta[index+1]
			index+=2

		if self.data["A"] == "Free":
			self.A = self.myoutput.beta[index]
			index+=1

		if self.data["B"] == "Free":
			self.B = self.myoutput.beta[index]


	def plot_fit(self):

		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.plot(self.H_M_ratio_dict[0], self.Pressures_dict[0][:len(self.H_M_ratio_dict[0])], 'ko', label= 'Veleckis et al.')
		plt.plot(self.H_M_ratio_dict[0], self.Pressures_fit[0], 'k-', label= 'Model derived')
		for i in range(1, len(self.Pressures_dict)):
			plt.plot(self.H_M_ratio_dict[i], self.Pressures_dict[i][:len(self.H_M_ratio_dict[i])], 'ko')
			plt.plot(self.H_M_ratio_dict[i], self.Pressures_fit[i], 'k-')

		ax.set_xlabel(r'Hydrogen to metal ratio $(\mathrm{\frac{H}{M}})$', fontsize=16)
		ax.set_ylabel('$\mathregular{\sqrt{Pressure (Torr)}}$', fontsize=16)
		ax.legend(fontsize=12)
		plt.tight_layout()
		plt.show()


	def display_fitting_constants(self):

		T = self.Temperatures
		print("E {0:.5f}".format(self.E))
		print("Alpha {0:.5f}".format(self.Alpha))
		print("Beta {0:.5f}".format(self.Beta))
		print("Gamma {0:.5f}".format(self.Gamma))
		
		if self.data["Temp_dependence"] == "True":
			print("R_0 {0:.5f}".format(self.R_0))
			print("R_1 {0:.5f}".format(self.R_1))
			print("R_2 {0:.10f}".format(self.R_2))
			R = [self.R_0]*len(T) + self.R_1*T + self.R_2*(T**2)
		else:
			print("R_0 {0:.5f}".format(self.R_0))
			R = [self.R_0]*len(T)

		print("A {0:.5f}".format(self.A))
		print("B {0:.5f}".format(self.B))
		print(R)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.plot(T, R, 'ko')
		ax.set_xticks(np.arange(600, 1000, step=50)) 
		ax.tick_params(axis='y', labelsize=16)
		ax.tick_params(axis='x', labelsize=16)
		ax.set_xlabel(r'Temperature $(^{\circ} K)$', fontsize=16)
		ax.set_ylabel('R', fontsize=16)
		plt.tight_layout()
		plt.show()

	def plot_enthalpy(self, referenced): 

		x = np.linspace(0, 0.7, 10)
		#print(self.E, self.Alpha, self.Beta)
		H_exp = [0.01010101, 0.020408163, 0.030927835, 0.041666667, 0.052631579, 0.063829787, 0.075268817, 0.086956522, 0.098901099, 0.111111111, 0.123595506, 0.136363636, 0.149425287, 0.162790698, 0.176470588, 0.19047619, 0.204819277, 0.219512195, 0.234567901, 0.25, 0.265822785, 0.282051282, 0.298701299, 0.315789474, 0.333333333, 0.351351351, 0.369863014, 0.388888889, 0.408450704, 0.428571429, 0.449275362, 0.470588235, 0.492537313, 0.515151515, 0.538461538, 0.5625, 0.587301587, 0.612903226, 0.639344262]
		H_exp = np.asarray(H_exp)
		Enthalpy_veleckis = [-0.371620704, -0.375089742, -0.376824261, -0.37855878, -0.380726929, -0.382027819, -0.385930486, -0.388965895, -0.392001303, -0.39677123, -0.399806639, -0.402408417, -0.405010196, -0.408045604, -0.411081012, -0.413682791, -0.416718199, -0.420187237, -0.422355386, -0.424957165, -0.428426203, -0.431461611, -0.435364279, -0.438399687, -0.441435096, -0.443603244, -0.446638653, -0.448806802, -0.452709469, -0.456178508, -0.459213916, -0.463116584, -0.466151992, -0.46962103, -0.471355549, -0.473523698, -0.474824587, -0.475691847, -0.476125476]
		Enthalpy_veleckis=np.asarray(Enthalpy_veleckis)

		x_DFT = np.asarray([(i+1)/16 for i in range(8)])
		Enthalpy_DFT = [-0.3858700000000046, -0.39626916383213784, -0.3900254937238336, -0.41404223028105536, -0.41129553623267656, -0.41353207729712754, -0.42765350224292015, -0.4269243554026301]
		Enthalpy_DFT = np.asarray(Enthalpy_DFT)
		Enthalpy_DFT = Enthalpy_DFT - Enthalpy_DFT[0]
		Enthalpy_DFT = Enthalpy_DFT -0.01

		fig = plt.figure()
		ax = fig.add_subplot(111)

		H = self.E + 2*self.Alpha*x + 3*self.Beta*(x**2) + 4*self.Gamma*(x**3)

		if referenced:

			H-=self.E
			Enthalpy_veleckis-=Enthalpy_veleckis[0]

		ax.plot(x, H, 'k--', label='Model derived')
		ax.plot(H_exp, Enthalpy_veleckis, 'ko', label= 'Veleckis et al.')

		#ax.errorbar(x_DFT, np.asarray(Enthalpy_DFT), xerr=0, yerr=0.01, fmt='o', color='k', label='DFT calculation')
		ax.set_xlabel(r'Hydrogen to metal ratio $(\mathrm{\frac{H}{M}})$', fontsize=16)
		ax.set_ylabel(r'$\Delta H_{\mathrm{x}}^{\circ} - \Delta H^{\circ} \ (\frac{eV}{\mathrm{H}\:atom})$', fontsize=16)
		ax.legend(loc=1)
		plt.tight_layout()
		plt.show()


fit_object = Isotherm_Fit()

#fit_object.plot_raw_data()

fit_object.ODR_fit()

#fit_object.test_func(np.asarray([1,2,3]))

#fit_object.display_fitting_constants()

fit_object.plot_enthalpy(True)