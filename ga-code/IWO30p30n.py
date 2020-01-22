import random,os,math,copy,STLgen,subprocess
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
MAX_POPULATION_SIZE = 10 # maximum number of plants in a colony(or population)
#STANDARD Deviation
#SD for slat
sigma_fin1 = 0.0005  # final standard deviation
sigma_ini1 = 6  # initial standard deviation
#SD for flap
sigma_fin2 = 0.0005  # final standard deviation
sigma_ini2 = 1.5  # initial standard deviation

Smin = 0  # min seeds produced
Smax = 3  # max seeds produced
n_mi = 3  # modulation index
iter_max = 40  # Maximum number of iterations to be done
# slat config[[hori dist,vert dist],[deflection,0]] + flap config[[hori dist,vert dist],[deflection,0]] + output
CHROMOSOME_SIZE = 2+2+1
#design space
max_d_s=-15
min_d_s=-45
max_o_s=-2.25
min_o_s=-4.25
max_v_s=-2.2
min_v_s=2.2
max_d_f=45
min_d_f=10
max_o_f=4.5
min_o_f=0.6
max_v_f=4.5
min_v_f=1.3
#max and min lines to be read from cfd file
MAX_Line=10000
MIN_Line=9500

aoa=math.radians(16)
#filepaths
cfdstl="/home/fmg/OpenFOAM/fmg-6/run/30P30N/cfd_openfoam/wingMotion_snappyHexMesh/constant/triSurface/air.stl" #add path for simulation stl file
cu_file=r"/home/fmg/OpenFOAM/fmg-6/run/30P30N/curves/" #add path for curves location
cfd_call="/home/fmg/OpenFOAM/fmg-6/run/30P30N/cfd_openfoam/Allrun" # add path for cfd allrun code
gcost="/home/fmg/OpenFOAM/fmg-6/run/30P30N/simvalue/CP%i/forces.dat" #add path for point where data has to be taken
delete="/home/fmg/OpenFOAM/fmg-6/run/30P30N/ga-code/delete" #pleasse add path of delete file here

#importing initial profile
Init_Flap_Profile=pd.read_csv("flap1.csv")#flap
Init_Flap_Profile=Init_Flap_Profile.values
Init_ME_Profile=pd.read_csv("me.csv") #main element
Init_ME_Profile=Init_ME_Profile.values/4
trailing_edge_me=np.array([0.874007,0.030943])
leading_edge_me=np.array([0.043807,-0.016668])
Init_Slat_Profile=pd.read_csv("slat1.csv") #slat
Init_Slat_Profile=Init_Slat_Profile.values
c_=np.array([[0.018772,0.005079]])
Init_Slat_profile=np.append(Init_Slat_Profile,c_,axis=0)

# class that generates chromosomes
def create_file(pathname, filename, openmode):
    filepath = os.path.join(pathname, filename)
    if not os.path.exists(pathname):
        os.makedirs(pathname)
    file1 = open(pathname+filename, openmode)
    return file1


class Chromosome:
    def __init__(self, chrom_num=0, gen_number=0, mode=" ",):
        # creating initial gene
        self._genes = np.zeros((CHROMOSOME_SIZE, 2), dtype=float)
        self.I = (100*gen_number) + chrom_num
        if(mode == "initialise"):
            flag=False
            while(not flag):
                self._genes[0][0] = min_o_s + (max_o_s-min_o_s) * random.random()
                self._genes[0][1] = min_v_s + ((max_v_s-min_v_s) * random.random())
                self._genes[1][0] = min_d_s + ((max_d_s-min_d_s) * random.random())
                self._genes[2][0] = min_o_f + ((max_o_f-min_o_f) * random.random())
                self._genes[2][1] = min_v_f + ((max_v_f-min_v_f) * random.random())
                self._genes[3][0] = min_d_f + ((max_d_f-min_d_f) * random.random())
                print("SEED#:",chrom_num)
                curve=self.profile_gen()
                Chromosome.genstl(curve,cfdstl)
                time=self.cal_cost()
                print("TIME TAKEN:",time," mins")
                temp=self.get_cost()
                if(temp==False):
                    flag=False
                else:
                    self._genes[-1][1] = temp
                    print("acceptable individual")
                    print("\nCOST: ",self._genes[-1][1])
                    #print("\nCurve Coordinates\n", curve)
                    break

    def cal_cost(self):
        start=time.time()
        subprocess.call([cfd_call,str(self.I)]) #ADD ALLRUN FILENAME HERE
        end=time.time()
        return (end-start)/60

        #for debugging
        # return random.randint(1,6000)

    def get_cost(self):
        cfd_file = open(gcost% self.I, "r")  #ADD FORCECOEFF FILE PATH HERE
        list=cfd_file.readlines()
        #actual code
        forcex=0
        forcey=0
        if(len(list)>=MAX_Line):
            forcex=0
            forcey=0
            for i in range(MIN_Line,MAX_Line):
                forcex+=float(list[i].split('(')[2].split()[0])
                forcey+=float(list[i].split('(')[2].split()[1])

            forcex/=(MAX_Line-MIN_Line)
            forcey/=(MAX_Line-MIN_Line)
            LDcoeffmat=np.array([[math.cos(aoa) ,-math.sin(aoa)],[math.sin(aoa),math.cos(aoa)]])
            LDmat=np.dot(LDcoeffmat,np.array([forcey,forcex]).reshape(2,1))
            LDratio=LDmat[0]/LDmat[1]
            return LDratio.squeeze()

        else:
            print("SIMULATION INCOMPLETE........REGENERATING INDIVIDUAL")
            return False

        #for debugging

        # return random.randint(1,6000)
    def get_genes(self):
        return self._genes
    def profile_gen(self):
        curve = self.transform_profile(self._genes)  # transforming slat

        return curve
    @staticmethod
    def genstl(curve,filename):
        STLgen.STL_Gen(curve[0].T[0],curve[0].T[1],'new_slat.stl')
        STLgen.STL_Gen(curve[1].T[0],curve[1].T[1],'me.stl')
        STLgen.STL_Gen(curve[2].T[0],curve[2].T[1],'new_flap.stl')
        #BELOW LINE TAKES THREE ARGUMENTS FIRST:SLAT FILE, SECOND: AIRFOIL FILE , THIRD: COMBINED FILE PATH along with name
        STLgen.combine('new_slat.stl','me.stl','new_flap.stl',filename)
        print("profile stl generated and saved")
    @staticmethod
    def transform_slat(slat_prof,angle,overhang,vert):
        overhang/=100
        vert/=100
        max=-1000
        j=0
        for i in range(slat_prof.shape[0]):
            if(slat_prof[i][0]>max):
                max=slat_prof[i][0]
                j=i
        fin=copy.deepcopy(slat_prof[i])
        slat_prof.T[0]-=fin[0]
        slat_prof.T[1]-=fin[1]
        slat_prof=np.append(slat_prof,np.zeros((slat_prof.shape[0],1)),axis=1)
        slat_prof=STLgen.rotate(slat_prof,angle+30)
        slat_prof.T[0]+=leading_edge_me[0]+overhang
        slat_prof.T[1]+=leading_edge_me[1]-vert#rotate end
        slat_prof=np.delete(slat_prof,2,1)
        return slat_prof
    @staticmethod
    def transform_flap(flap_prof,angle,overhang,vert):
        overhang/=100
        vert/=100
        min,j=100000,0
        for i in range(flap_prof.shape[0]):
            if(flap_prof[i][0]<min):
                min=flap_prof[i][0]
                j=i
        fin=copy.deepcopy(flap_prof[j])
        flap_prof.T[0]-=fin[0] #rotate start
        flap_prof.T[1]-=fin[1]
        flap_prof=np.append(flap_prof,np.zeros((flap_prof.shape[0],1)),axis=1)
        flap_prof=STLgen.rotate(flap_prof,angle-30)
        flap_prof.T[0]+=trailing_edge_me[0]+overhang
        flap_prof.T[1]+=trailing_edge_me[1]-vert #rotate end
        flap_prof=np.delete(flap_prof,2,1)
        return flap_prof

    @staticmethod
    def transform_profile(genes):
        #slat_transformation
        slat=Chromosome.transform_slat(Init_Slat_Profile,genes[1][0],genes[0][0],genes[0][1])
        flap=Chromosome.transform_flap(Init_Flap_Profile,genes[3][0],genes[2][0],genes[2][1])
        slat/=4
        flap/=4
        fin_prof=[slat,Init_ME_Profile,flap]
        return fin_prof
    def __str__(self):
        return self._genes.__str__()
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# class that create one set of generations
class Population:
    def __init__(self, size, gen_num=0, mode=" "):
        self._chromosomes = []
        i = 0
        gen_num = str(gen_num)
        while i < size:
            self.add_chromosomes(Chromosome(i+1, int(gen_num), mode))
            i += 1

    def add_chromosomes(self, chromosome):
        self._chromosomes.append(chromosome)

    def get_chromosomes(self):
        pop_chroms_2d_array = np.zeros(
            (len(self._chromosomes), CHROMOSOME_SIZE, 2), dtype=float)
        #pop_chroms_2d_array = np.around(pop_chroms_2d_array, 6)
        for i in range(len(self._chromosomes)):
            pop_chroms_2d_array[i] = self._chromosomes[i].get_genes()
        return pop_chroms_2d_array
# class that helps in evolving and mutating the genes of the chromosomes


class GeneticAlgorithm:
    @staticmethod
    def truncated_normal(mean, stddev, minval, maxval):
        return np.clip(np.random.normal(mean, stddev), minval, maxval)
    @staticmethod
    def reproduce(pop, iter):
        new_pop = copy.deepcopy(pop)
        worst_cost = pop._chromosomes[-1].get_genes()[CHROMOSOME_SIZE - 1][1]
        best_cost = pop._chromosomes[0].get_genes()[CHROMOSOME_SIZE - 1][1]
        sigma_iter1 = GeneticAlgorithm.std_deviation(iter, iter_max , sigma_ini1 , sigma_fin1) #for shape
        sigma_iter2 = GeneticAlgorithm.std_deviation(iter, iter_max , sigma_ini2 , sigma_fin2) #for position
        if(best_cost != worst_cost):
            seed_num = 0
            for i in range(MAX_POPULATION_SIZE): #limiting the number of individuals that can reproduce
                ratio = (pop._chromosomes[i]._genes[- 1][1] - worst_cost) / (best_cost - worst_cost)
                # number of seeds chromosome can produce on the basis of rank
                S = Smin + (Smax - Smin) * ratio
                #print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS",int(S))
                for j in range(int(S)):
                    seed_num += 1
                    seed = Chromosome(seed_num, iter)
                    flag=False
                    while(not flag):
                        seed._genes[0][0] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[0][0],sigma_iter1,min_o_s,max_o_s)
                        seed._genes[0][1] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[0][1],sigma_iter1,min_v_s,max_v_s)
                        seed._genes[1][0] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[1][0],sigma_iter1,min_d_s,max_d_s)
                        seed._genes[2][0] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[2][0],sigma_iter2,min_o_f,max_o_f)
                        seed._genes[2][1] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[2][1],sigma_iter2,min_v_f,max_v_f)
                        seed._genes[3][0] = GeneticAlgorithm.truncated_normal(pop._chromosomes[i].get_genes()[3][0],sigma_iter2,min_d_f,max_d_f)

                        flag=True
                        print("SEED# ",seed_num)
                        print("Slope Check Passed")
                        curve=seed.profile_gen()
                        Chromosome.genstl(curve,cfdstl)
                        time = seed.cal_cost()
                        print("\nTIME TAKEN:",time," mins")
                        temp=seed.get_cost()
                        if(temp==False):
                            flag= False
                        else:
                            seed._genes[-1][1] = temp
                            #print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS",temp)
                            print("ACCEPTABLE INDIVIDUAL")
                            print("\nCOST: ",seed._genes[-1][1])
                            #print("\nCurve Coordinates\n", curve)
                            new_pop.add_chromosomes(seed)


            GeneticAlgorithm.sort(new_pop)
            for i in range(MAX_POPULATION_SIZE):
                pop._chromosomes[i] = new_pop._chromosomes[i]
        else:
            print("best and worst cost equal can`t reproduce")
            return False,False
        return pop,new_pop

    @staticmethod
    def std_deviation(iter, iter_max,sigma_ini,sigma_fin):
        sigma_iter = (((iter_max - iter)**n_mi) / iter_max **
                      n_mi) * (sigma_ini - sigma_fin) + sigma_fin
        return sigma_iter

    @staticmethod
    def sort(pop):
        pop_chroms_2d_array = pop.get_chromosomes()
        #print("chroms",pop.get_chromosomes())
        sindices = np.argsort(pop_chroms_2d_array[:, :, 1][:, -1], axis=0)
        print("sindices", sindices)
        sorted_chroms = Population(len(pop._chromosomes), 0, "zeros")
        for i in range(0, len(pop._chromosomes)):
            sorted_chroms._chromosomes[i]._genes = pop_chroms_2d_array[sindices[-(i+1)]]
        for i in range(0, len(pop._chromosomes)):
            pop._chromosomes[i] = sorted_chroms._chromosomes[i]

#------------------------------------------------------------------------------------------------------------------------------------#-



#------------------------------------------------------------------------------------------------------------------------------------#-


def _print_population(new_pop, gen_number):
    chroms = new_pop.get_chromosomes()
    print("\n---------------------------------------------------------")
    print("PRINTING AFTER SORTING THE POPULATION")
    print("Generation#", gen_number, "|Fittest chromosome fitness:",chroms[0][-1][1])
    print("-----------------------------------------------------------")
    #dplot1.update(gen_number,chroms[0][-1][1])

    i = 0
    for i in range(MAX_POPULATION_SIZE):
        print("PLANT#", i + 1, " ", "||Fitness:",chroms[i][-1][1], "\n")
        print("Parameters\n",chroms[i])
        #handle.write("PLANT NO%i")
        print("--------------------------------------------------------------")
    i=0
    for i in range(len(new_pop._chromosomes)):
        I = (100*gen_number)+i+1
        curve_file = create_file(cu_file,"curve_%04i.dat" % I, 'w+')  # saving curve coordinates
        for j in range(CHROMOSOME_SIZE):
            curve_file.write(str(new_pop._chromosomes[i]._genes[j])+'\n')
        curve_file.close()
        plt.figure()
        curve=new_pop._chromosomes[i].profile_gen()
        Chromosome.genstl(curve,cu_file+"cstl_%04i.stl"%I)
        plt.plot(curve[0].transpose()[0],curve[0].transpose()[1])
        plt.plot(curve[1].transpose()[0],curve[1].transpose()[1])
        plt.plot(curve[2].transpose()[0],curve[2].transpose()[1])
        plt.axes().set_aspect('equal')
        # plt.xlim(-0.04,0.05)
        # plt.ylim(-0.05,0.07)
        plt.savefig(cu_file+"fig_%04i"%I)
    print("Curve File Saved")
#-------------------------------------------------------------------------





# main
# initialising population
s=time.time()
subprocess.call([delete])
print("Running for.......\nMAX_ITERATION:",iter_max)
print("POPULATION SIZE:",MAX_POPULATION_SIZE)
print("GENERATION#0-----------------------------------------------------------------------------------")
population = Population(MAX_POPULATION_SIZE, 0, "initialise")
#dplot1=plotter.DynamicUpdate()#plotting maximum fitness dynamically
GeneticAlgorithm.sort(population)
#handle=open(r"COST.dat",'w+')
#handle.write("Generation# 0 \n")
_print_population(population, 0)
iter = 1
sigma = [GeneticAlgorithm.std_deviation(0, iter_max,sigma_ini1,sigma_fin1)]
print("**************************************EVOLUTION STARTED*********************************************************")
while iter <= iter_max:
    sigma.append(GeneticAlgorithm.std_deviation(iter, iter_max,sigma_ini1,sigma_fin1))
    print("GENERATION#",iter,"-----------------------------------------------------------------------------------")
    print("sigma shape:",GeneticAlgorithm.std_deviation(iter,iter_max,sigma_ini1,sigma_fin1))
    print("sigma position:",GeneticAlgorithm.std_deviation(iter,iter_max,sigma_ini2,sigma_fin2))
    population,new_pop = GeneticAlgorithm.reproduce(population, iter)
    print("*************************************REPRODUCED*********************************************")
    if(population == False):
        iter+=1
        break
    _print_population(new_pop, iter)
    iter += 1
plt.figure()
plt.plot(np.arange(1, iter+1), sigma)
e=time.time()
print("TOTAL TIME TAKEN: ",((e-s)/60),"mins")
#plt.show()
