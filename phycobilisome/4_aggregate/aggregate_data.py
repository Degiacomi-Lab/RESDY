
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

reject_minor = True
reject_major = True


#################################
### GATHER CHAINS INFORMATION ###
#################################

# get information on problematic calculations
problem_minor = np.loadtxt("../pka/problem_minor.dat", dtype=str)
problem_major = np.loadtxt("../pka/problem_major.dat", dtype=str)

# get information on identical chains
chaindata = {}
fin = open("chains_raw.dat", "r")
fout = open("chains_summary.dat", "w")
while True:

    key = fin.readline().strip('\n').strip('\t')
    if not key:
        break

    value = np.array(fin.readline().split(" ")[::3])

    if reject_major:
        test = np.isin(value, problem_major)
        value = value[~test]

    if reject_minor:
        test = np.isin(value, problem_minor)
        value = value[~test]


    chaindata[key] = value
    line = fin.readline()
    if not line:
        break

    fout.write("%s\t%s\n"%(key, " ".join(value)))

fout.close()

######################
### AGGREGATE PKAs ###
######################

sasa = np.loadtxt("../sasa/sasa.dat", dtype=str)
sasa_index = sasa[:, 0:2]
sasa_vals = sasa[:, 2].astype(float)


######################
### AGGREGATE PKAs ###
######################

# read pKa data
# line format example: HN	LYS	32	10.20"
pkas = np.loadtxt("../result.dat", dtype=str)
pka_index = pkas[:, ::2]
pka_vals = pkas[:, -1].astype(float)

# aggregate data according to chain information
result_index = []
result_values = []
fout = open("results_aggregated.dat", "w")
for k in list(chaindata):

    nb_prots = len(chaindata[k])
    lineout="\n%s (%i)"%(k, nb_prots)
    fout.write("%s\n"%lineout)
    print(lineout)

    test = np.isin(sasa_index[:, 0], chaindata[k])   
    resid_sasa = sasa_index[test, 1].astype(int)
    values_sasa = sasa_vals[test]
   
    test = np.isin(pka_index[:, 0], chaindata[k])
    resid_pka = pka_index[test, 1].astype(int)
    values_pka = pka_vals[test]
    for r in np.unique(resid_pka):

        # find pKa of relevant residue
        test_pka = resid_pka == r
        lys_pkas = values_pka[test_pka]

        # find SASA of relevant residue
        test_sasa = resid_sasa == r
        lys_sasa = values_sasa[test_sasa]

        result_index.append([k, r])
        result_values.append([np.mean(lys_pkas), np.std(lys_pkas), len(lys_pkas), np.mean(lys_sasa), np.std(lys_sasa), len(lys_sasa)])

        lineout = "K%i. pKa: %4.2f pm %4.2f (%i). SASA: %4.2f pm %4.2f (%i)"%\
                  (r, np.mean(lys_pkas), np.std(lys_pkas),len(lys_pkas),np.mean(lys_sasa), np.std(lys_sasa), len(lys_sasa))
        fout.write("%s\n"%lineout)
        print(lineout)


fout.close()

result_index = np.array(result_index)
result_values = np.array(result_values)


######################
### GET TOP RANKED ###
######################

minsasa = np.mean(sasa_vals)-np.std(sasa_vals) # CRITERION FOR SASA ACCEPTANCE!

pos = np.argsort(result_values[:, 0])
result_values_sorted = result_values[pos]
result_index_sorted = result_index[pos]

cutoff = 20
minmeas = 1

fout = open("top20.dat", "w")
lineout = "\n### TOP %s (with >= %s measures and SASA > %i) ###\n"%(cutoff, minmeas, int(minsasa))
fout.write("%s\n"%lineout)
print(lineout)

cnt = 0
for i in range(len(result_values_sorted)):

    if result_values_sorted[i, 2] >= minmeas and result_values_sorted[i, 3] > minsasa:
        cnt += 1
        lineout = "%s (%i), K%s. pKa: %4.2f pm %4.2f (%i). SASA: %4.2f pm %4.2f (%i)"%\
             (result_index_sorted[i, 0], len(chaindata[result_index_sorted[i, 0]]), result_index_sorted[i, 1],\
              result_values_sorted[i, 0], result_values_sorted[i, 1], result_values_sorted[i, 2],\
              result_values_sorted[i, 3], result_values_sorted[i, 4], result_values_sorted[i, 5])
        fout.write("%s\n"%lineout)
        print(lineout)

    if cnt == cutoff:
        break

fout.close()



#######################
### MAKE SOME PLOTS ###
#######################

hist, bin_edges = np.histogram(sasa_vals, bins=50)
bin_centres = (bin_edges[:-1]+bin_edges[1:])/2.0
bin_width = bin_edges[1]-bin_edges[0]


fig = plt.figure(1, figsize=(4,6))

ax1 = fig.add_subplot(2, 1, 1)
rect1 = patches.Rectangle((0, 0), minsasa, np.max(hist)*1.05, facecolor='gray', alpha=0.5)
ax1.add_patch(rect1)
ax1.set_ylim([0, np.max(hist)*1.05])
ax1.bar(bin_centres, hist, width=bin_width, color="maroon", edgecolor="k")
ax1.set_ylabel("count (#)")
ax1.set_xlim(xmin=0)


ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
ax2.plot(result_values[:, 3], result_values[:, 0], "k.")
#ax2.errorbar(result_values[:, 3], result_values[:, 0], xerr=result_values[:, 4], yerr=result_values[:, 1], fmt='o')
ax2.set_xlabel("SASA ($\AA^2$)")
ax2.set_xlim(xmin=0)
ax2.set_ylabel("pKa")
rect2 = patches.Rectangle((0, 0), minsasa, np.max(result_values[:, 0])*1.05,  facecolor='gray', alpha=0.5)
ax2.add_patch(rect2)
ax2.set_ylim([np.min(result_values[:, 0])*0.95, np.max(result_values[:, 0])*1.05])

plt.subplots_adjust(hspace=0)

plt.savefig("lysine_stats.png")
plt.show()

