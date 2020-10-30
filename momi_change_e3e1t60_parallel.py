import momi
import numpy
import pathlib
import os
import multiprocessing
import subprocess

NeAlb=1e3
NeBot=1e1
tdecline=60
gt=1

def run_sim(simrep):
	#define model
	model_bot = momi.DemographicModel(N_e=NeAlb, gen_time=gt, muts_per_gen=4.5e-8)
	model_bot.add_leaf("ALB",N=NeAlb,t=110)
	model_bot.add_leaf("CON",N=NeBot)
	model_bot.move_lineages("CON","ALB",t=110.01)
	model_bot.set_size("CON",N=NeAlb,t=tdecline)

	# name of sim set + output files. export simdir to os for running scripts
	curpath=pathlib.Path().parent.absolute()
	it=simrep
	expa=int(numpy.log10(NeAlb))
	expb=int(numpy.log10(NeBot))
	simset="momi_bot_10e"+str(expa)+"10e"+str(expb)+"_it"+str(it)
	simdir=str(curpath)+"/momi_sims/momi_bot_sims/" + simset
	print(simdir)
	os.environ['simdir']=simdir
	print(os.getenv('simdir'))

	#simulating data - 25 chromosomes, 200 rad loci * 300 bp each per chromosome. 
	#not completely realistic - could simulate whole chromosomes, create bed file specifying RAD loci, then extract SFS?
	#Think about recomb rate! 

	recoms_per_gen = 1.25e-8
	bases_per_locus = int(6e5)
	n_loci = 25
	ploidy = 2

	# n_alleles per population (n_individuals = n_alleles / ploidy)
	sampled_n_dict = {"ALB":200, "CON":200}

	# create data directory if it doesn't exist
	pathlib.Path(simdir).mkdir(parents=True, exist_ok=True)

	# simulate 20 "chromosomes", saving each in a separate vcf file
	for chrom in range(1, n_loci+1):
		  model_bot.simulate_vcf(
				f"{simdir}/{chrom}",
				recoms_per_gen=recoms_per_gen,
				length=bases_per_locus,
				chrom_name=f"chr{chrom}",
				ploidy=ploidy,
				random_seed=1234+chrom,
				sampled_n_dict=sampled_n_dict,
				force=True)

	# need a dict mapping samples to populations
	ind2pop = {}
	for pop, n in sampled_n_dict.items():
		for i in range(int(n / ploidy)):
			# in the vcf, samples are named like YRI_0, YRI_1, CHB_0, etc
			ind2pop["{}_{}".format(pop, i)] = pop

	with open(simdir+"/ind2pop.txt", "w") as f:
		for i, p in ind2pop.items():
			print(i, p, sep="\t", file=f)

	# can't comment in commands with %%sh so I will explain these next two steps here
	# first compute allele counts with momi.read_vcf
	# usage: python -m momi.read_vcf $VCF $IND2POP $OUTFILE --bed $BED:
	# then extract combined SFS, split into 100 blocks for jackknifing / bootstrap
	# format: python -m momi.extract_sfs $OUTFILE $NBLOCKS $COUNTS...

	subprocess.run('for chrom in `seq 1 25`; do python -m momi.read_vcf ${simdir}/${chrom}.vcf.gz ${simdir}/ind2pop.txt ${simdir}/${chrom}.snpAlleleCounts.gz --bed ${simdir}/${chrom}.bed;done',shell=True)
	subprocess.run('python -m momi.extract_sfs ${simdir}/sfs.gz 100 ${simdir}/*.snpAlleleCounts.gz',shell=True)

	# read sfs into python!
	sfsfile=simdir+"/sfs.gz"
	sfs = momi.Sfs.load(sfsfile)

	#set model for inference - constant pop size, contemporary samples only!
	model_inf_constant_contemp =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)

	#add data to model
	model_inf_constant_contemp.set_data(sfs)

	#define parameters to infer
	model_inf_constant_contemp.add_size_param("n_constant")

	#add contemp pop
	model_inf_constant_contemp.add_leaf("CON",N="n_constant")

	#infer parameters
	model_inf_constant_contemp.optimize(method="TNC")

	#set model for inference - contemporary samples, 2 different pop sizes assumed!
	model_inf_change_contemp =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)

	#add data to model
	model_inf_change_contemp.set_data(sfs)

	#define parameters to infer - model with size change and unknown time of decline
	model_inf_change_contemp.add_size_param("n_alb")
	model_inf_change_contemp.add_size_param("n_bot")
	model_inf_change_contemp.add_time_param("t_bot",upper=1e2)

	#add contemp pop and size change event

	model_inf_change_contemp.add_leaf("CON",N="n_bot")
	model_inf_change_contemp.set_size("CON", N="n_alb", t="t_bot")

	#infer parameters
	model_inf_change_contemp.optimize(method="TNC")

	#set model for inference - albatross +contemporary samples, constant pop size assumed!
	model_inf_constant_temporal =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)

	#add data to model
	model_inf_constant_temporal.set_data(sfs)

	#define parameters to infer - model with constant size
	model_inf_constant_temporal.add_size_param("n_constant")

	#add contemporary and historic pops
	model_inf_constant_temporal.add_leaf("ALB",N="n_constant",t=110)
	model_inf_constant_temporal.add_leaf("CON",N="n_constant")
	model_inf_constant_temporal.move_lineages("CON","ALB",t=110.01)

	#infer parameters
	model_inf_constant_temporal.optimize(method="TNC")


	#set model for inference - albatross + contemporary samples, 2 different pop sizes assumed!
	model_inf_change_temporal =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)

	#add data to model
	model_inf_change_temporal.set_data(sfs)

	#define parameters to infer - model with population size change and unknown time of decline
	model_inf_change_temporal.add_size_param("n_alb")
	model_inf_change_temporal.add_size_param("n_bot")
	model_inf_change_temporal.add_time_param("t_bot",upper=1e2)

	#add populations and size change even
	model_inf_change_temporal.add_leaf("CON",N="n_bot")
	model_inf_change_temporal.add_leaf("ALB",N="n_alb",t=110)
	model_inf_change_temporal.move_lineages("CON","ALB",t=110.01)
	model_inf_change_temporal.set_size("CON", N="n_alb", t="t_bot")

	#infer parameters
	model_inf_change_temporal.optimize(method="TNC")

	#bootstrapping
	n_bootstraps = 10
	# make copies of the original models to avoid changing them
	model_inf_cons_cont_copy = model_inf_constant_contemp.copy()
	model_inf_change_cont_copy = model_inf_change_contemp.copy()
	model_inf_cons_temp_copy = model_inf_constant_temporal.copy()
	model_inf_change_temp_copy = model_inf_change_temporal.copy()

	bootstrap_cons_cont = []
	bootstrap_change_cont = []
	bootstrap_cons_temp = []
	bootstrap_change_temp = []


	for i in range(n_bootstraps):

		# resample the data
		resampled_sfs = sfs.resample()

		print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for constant + contemp")

		# tell model to use the new dataset
		model_inf_cons_cont_copy.set_data(resampled_sfs)
		# choose new random parameters for submodel, optimize
		model_inf_cons_cont_copy.set_params(randomize=True)
		model_inf_cons_cont_copy.optimize()
		# append results
		bootstrap_cons_cont.append(model_inf_cons_cont_copy.get_params())

		print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for change + contemp")

		# tell model to use the new dataset
		model_inf_change_cont_copy.set_data(resampled_sfs)
		# choose new random parameters for submodel, optimize
		model_inf_change_cont_copy.set_params(randomize=True)
		model_inf_change_cont_copy.optimize()
		# append results
		bootstrap_change_cont.append(model_inf_change_cont_copy.get_params())

		print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for constant + temporal")

		# tell model to use the new dataset
		model_inf_cons_temp_copy.set_data(resampled_sfs)
		# choose new random parameters for submodel, optimize
		model_inf_cons_temp_copy.set_params(randomize=True)
		model_inf_cons_temp_copy.optimize()
		# append results
		bootstrap_cons_temp.append(model_inf_cons_temp_copy.get_params())
	
		print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for change + temporal")

		# tell model to use the new dataset
		model_inf_change_temp_copy.set_data(resampled_sfs)
		# choose new random parameters for submodel, optimize
		model_inf_change_temp_copy.set_params(randomize=True)
		model_inf_change_temp_copy.optimize()
		# append results
		bootstrap_change_temp.append(model_inf_change_temp_copy.get_params())
	
	#output/compute estimates, log likelihood, and AIC
	cn=model_inf_constant_contemp.get_params().get('n_constant')
	nparam1=len(model_inf_constant_contemp.get_params())
	lik1=model_inf_constant_contemp.log_likelihood()
	AIC1=2*nparam1-2*lik1
	cnh=model_inf_change_contemp.get_params().get('n_alb')
	cnc=model_inf_change_contemp.get_params().get('n_bot')
	ct=model_inf_change_contemp.get_params().get('t_bot')
	nparam2=len(model_inf_change_contemp.get_params())
	lik2=model_inf_change_contemp.log_likelihood()
	AIC2=2*nparam2-2*lik2
	tn=model_inf_constant_temporal.get_params().get('n_constant')
	nparam3=len(model_inf_constant_temporal.get_params())
	lik3=model_inf_constant_temporal.log_likelihood()
	AIC3=2*nparam3-2*lik3
	tnh=model_inf_change_temporal.get_params().get('n_alb')
	tnc=model_inf_change_temporal.get_params().get('n_bot')
	tt=model_inf_change_temporal.get_params().get('t_bot')
	nparam4=len(model_inf_change_temporal.get_params())
	lik4=model_inf_change_temporal.log_likelihood()
	AIC4=2*nparam4-2*lik4


	#write estimates and model characteristics
	estimate=simdir+"/estimates.csv"
	f = open(estimate,"a")
	f.write("Model,Data,Nh,Nc,T,lnL,AIC"+"\n")
	f.write("Constant,Contemporary,"+str(cn)+","+str(cn)+",NA,"+str(lik1)+","+str(AIC1)+"\n")
	f.write("Change,Contemporary,"+str(cnh)+","+str(cnc)+","+str(ct)+","+str(lik2)+","+str(AIC2)+"\n")
	f.write("Constant,Temporal,"+str(tn)+","+str(tn)+",NA,"+str(lik3)+","+str(AIC3)+"\n")
	f.write("Change,Temporal,"+str(tnh)+","+str(tnc)+","+str(tt)+","+str(lik4)+","+str(AIC4))
	f.close()


	#write bootstrap results
	boot=simdir+"/bootstraps.csv"
	f = open(boot,"a")
	f.write("Model,Data,Param,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10"+"\n")
	f.write("Constant,Contemporary,Nconstant")
	for i in range(len(bootstrap_cons_cont)):
		f.write(',')
		f.write(str(bootstrap_cons_cont[i].get('n_constant')))
	f.write('\n')
	f.write("Change,Contemporary,Nhistoric")
	for i in range(len(bootstrap_change_cont)):
		f.write(',')
		f.write(str(bootstrap_change_cont[i].get('n_alb')))
	f.write('\n')
	f.write("Change,Contemporary,Ncontemporary")
	for i in range(len(bootstrap_change_cont)):
		f.write(',')
		f.write(str(bootstrap_change_cont[i].get('n_bot')))
	f.write('\n')
	f.write("Change,Contemporary,Tbottleneck")
	for i in range(len(bootstrap_change_cont)):
		f.write(',')
		f.write(str(bootstrap_change_cont[i].get('t_bot')))
	f.write('\n')
	f.write("Constant,Temporal,Nconstant")
	for i in range(len(bootstrap_cons_temp)):
		f.write(',')
		f.write(str(bootstrap_cons_temp[i].get('n_constant')))
	f.write('\n')
	f.write("Change,Temporal,Nhistoric")
	for i in range(len(bootstrap_change_temp)):
		f.write(',')
		f.write(str(bootstrap_change_temp[i].get('n_alb')))
	f.write('\n')
	f.write("Change,Temporal,Ncontemporary")
	for i in range(len(bootstrap_change_temp)):
		f.write(',')
		f.write(str(bootstrap_change_temp[i].get('n_bot')))
	f.write('\n')
	f.write("Change,Temporal,Tbottleneck")
	for i in range(len(bootstrap_change_temp)):
		f.write(',')
		f.write(str(bootstrap_change_temp[i].get('t_bot')))
	f.close()

a_pool = multiprocessing.Pool()
result = a_pool.map(run_sim, range(10))
