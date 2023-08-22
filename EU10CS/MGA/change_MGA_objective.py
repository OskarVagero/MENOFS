print('hello world')

def change_objective(technology,country,optimisation,slacklevel):
	with open(snakemake.input[0], "r") as f:
		list_of_lines = f.readlines()
	#print(list_of_lines[484])
	list_of_lines[483] = 'var_exist_pcap_z("'+str(country)+'","' + str(technology)+'")' + ' + var_new_pcap_z("'+str(country)+'","' + str(technology)+'")' +"\n"
	list_of_lines[761] = 'Solve Dispatch ' + str(optimisation) + ' objective using LP;' + '\n'
	#list_of_lines[484] = 'g("'+str(technology)+'");'+"\n"
	with open("/cluster/work/projects/ec85/MENOFS/EU10CS/results/" + country + "_" + technology + "_" + optimisation + "_" + slacklevel + "/highres_mga.gms", "w") as f:
		f.writelines(list_of_lines)


print(snakemake.wildcards.technology,snakemake.wildcards.country, snakemake.wildcards.optimisation)

change_objective(snakemake.wildcards.technology,snakemake.wildcards.country, snakemake.wildcards.optimisation, snakemake.wildcards.slacklevel)

