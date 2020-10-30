
# coding: utf-8

# In[105]:


import momi
import numpy
import pathlib
import os


# In[106]:


# here is a stable population model with historical sampling
# to actually generate data for this I had to add split event slightly before sampling

NeConstant=1e4

model_stab = momi.DemographicModel(N_e=NeConstant, gen_time=1,
                              muts_per_gen=4.5e-8)
model_stab.add_leaf("ALB",N=NeConstant,t=110)
model_stab.add_leaf("CON",N=NeConstant)
model_stab.move_lineages("CON","ALB",t=110.01)


# In[108]:


# name of sim set + output files. export simdir to os for running scripts
curpath=pathlib.Path().parent.absolute()

it=0

exp=int(numpy.log10(NeConstant))

simset="momi_stable_10e"+str(exp)+"_it"+str(it)

simdir=str(curpath)+"/momi_sims/momi_stable_sims/" + simset

print(simdir)

os.environ['simdir']=simdir

print(os.getenv('simdir'))


# In[109]:


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
      model_stab.simulate_vcf(
            f"{simdir}/{chrom}",
            recoms_per_gen=recoms_per_gen,
            length=bases_per_locus,
            chrom_name=f"chr{chrom}",
            ploidy=ploidy,
            random_seed=1234+chrom,
            sampled_n_dict=sampled_n_dict,
            force=True)


# In[110]:


# need a dict mapping samples to populations
ind2pop = {}
for pop, n in sampled_n_dict.items():
    for i in range(int(n / ploidy)):
        # in the vcf, samples are named like YRI_0, YRI_1, CHB_0, etc
        ind2pop["{}_{}".format(pop, i)] = pop

with open(simdir+"/ind2pop.txt", "w") as f:
    for i, p in ind2pop.items():
        print(i, p, sep="\t", file=f)


# In[111]:


# can't comment in commands with %%sh so I will explain these next two steps here
# first compute allele counts with momi.read_vcf
# usage: python -m momi.read_vcf $VCF $IND2POP $OUTFILE --bed $BED:
# then extract combined SFS, split into 100 blocks for jackknifing / bootstrap
# format: python -m momi.extract_sfs $OUTFILE $NBLOCKS $COUNTS...


# In[112]:


get_ipython().run_cell_magic(u'sh', u'', u'for chrom in `seq 1 25`;\ndo\n    python -m momi.read_vcf \\\n           ${simdir}/${chrom}.vcf.gz ${simdir}/ind2pop.txt \\\n           ${simdir}/${chrom}.snpAlleleCounts.gz \\\n           --bed ${simdir}/${chrom}.bed\ndone')


# In[113]:


get_ipython().run_cell_magic(u'sh', u'', u'python -m momi.extract_sfs ${simdir}/sfs.gz 100 ${simdir}/*.snpAlleleCounts.gz')


# In[114]:


# read sfs into python!
sfsfile=simdir+"/sfs.gz"
sfs = momi.Sfs.load(sfsfile)


# In[115]:


#set model for inference - constant pop size, contemporary samples only!
model_inf_constant_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=4.5e-8)


# In[116]:


#add data to model
model_inf_constant_contemp.set_data(sfs)


# In[117]:


#define parameters to infer
model_inf_constant_contemp.add_size_param("n_constant")


# In[118]:


model_inf_constant_contemp.add_leaf("CON",N="n_constant")


# In[119]:


model_inf_constant_contemp.optimize(method="TNC")


# In[120]:


#set model for inference - contemporary samples, 2 different pop sizes assumed!
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=4.5e-8)


# In[121]:


#add data to model
model_inf_change_contemp.set_data(sfs)


# In[122]:


#define parameters to infer - model with size change and unknown time of decline
model_inf_change_contemp.add_size_param("n_alb")
model_inf_change_contemp.add_size_param("n_bot")
model_inf_change_contemp.add_time_param("t_bot",upper=1e2)


# In[123]:


model_inf_change_contemp.add_leaf("CON",N="n_bot")
model_inf_change_contemp.set_size("CON", N="n_alb", t="t_bot")


# In[124]:


model_inf_change_contemp.optimize(method="TNC")


# In[125]:


#set model for inference - albatross +contemporary samples, constant pop size assumed!
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)


# In[126]:


#add data to model
model_inf_constant_temporal.set_data(sfs)


# In[127]:


#define parameters to infer - model with size change and unknown time of decline
model_inf_constant_temporal.add_size_param("n_constant")


# In[128]:


model_inf_constant_temporal.add_leaf("ALB",N="n_constant",t=109)
model_inf_constant_temporal.add_leaf("CON",N="n_constant")
model_inf_constant_temporal.move_lineages("CON","ALB",t=110)


# In[129]:


model_inf_constant_temporal.optimize(method="TNC")


# In[130]:


#set model for inference - albatross + contemporary samples, 2 different pop sizes assumed!
model_inf_change_temporal =  momi.DemographicModel(N_e=NeAlb, gen_time=1, muts_per_gen=4.5e-8)


# In[131]:


#add data to model
model_inf_change_temporal.set_data(sfs)


# In[132]:


#define parameters to infer - model with unknown time of decline
model_inf_change_temporal.add_size_param("n_alb")
model_inf_change_temporal.add_size_param("n_bot")
model_inf_change_temporal.add_time_param("t_bot",upper=1e2)


# In[133]:


model_inf_change_temporal.add_leaf("CON",N="n_bot")
model_inf_change_temporal.add_leaf("ALB",N="n_alb",t=109)
model_inf_change_temporal.move_lineages("CON","ALB",t=110)
model_inf_change_temporal.set_size("CON", N="n_alb", t="t_bot")


# In[134]:


model_inf_change_temporal.optimize(method="TNC")


# In[135]:


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
    


# In[136]:


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


# In[137]:


estimate=simdir+"/estimates.csv"
f = open(estimate,"a")
f.write("Model,Data,Nh,Nc,T,lnL,AIC"+"\n")
f.write("Constant,Contemporary,"+str(cn)+","+str(cn)+",NA,"+str(lik1)+","+str(AIC1)+"\n")
f.write("Change,Contemporary,"+str(cnh)+","+str(cnc)+","+str(ct)+","+str(lik2)+","+str(AIC2)+"\n")
f.write("Constant,Temporal,"+str(tn)+","+str(tn)+",NA,"+str(lik3)+","+str(AIC3)+"\n")
f.write("Change,Temporal,"+str(tnh)+","+str(tnc)+","+str(tt)+","+str(lik4)+","+str(AIC4))
f.close()


# In[138]:


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

