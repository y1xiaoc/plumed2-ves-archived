import numpy as np

total_files=70

for i in range(total_files):
	file_name="fes_" + str(i) + ".dat"
	matrix=np.genfromtxt(file_name)
	minimum=np.amin(matrix[50:90,1])
	maximum=np.amin(matrix[90:130,1])
	print(str(i) + " " + str(minimum) + " " + str(maximum) + " " + str(maximum-minimum))
	#print(i,minimum,maximum,maximum-minimum)
