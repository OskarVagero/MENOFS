print('Hello World')

def change_slack(technology,country,optimisation,slacklevel):
	with open(snakemake.input[0], "r") as f:
		list_of_lines = f.readlines()
	list_of_lines[7] = slacklevel + "\n"
	with open("/cluster/work/projects/ec85/MENOFS/EU10CS/results/" + country + "_" + technology + "_" + optimisation + "_" + slacklevel + "/mga_parameters.dd", "w") as f:
		f.writelines(list_of_lines)

print(snakemake.wildcards.technology,snakemake.wildcards.country, snakemake.wildcards.optimisation)

change_slack(snakemake.wildcards.technology,snakemake.wildcards.country,snakemake.wildcards.optimisation,snakemake.wildcards.slacklevel)

