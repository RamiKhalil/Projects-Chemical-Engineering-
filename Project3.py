#copyright Rami Yousef Khalil 260558325
import numpy as np
import scipy.sparse as sparse

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from time import sleep
import math
import pdb


a = np.array([1, 2, 3])
dT = 0.01 #temperature increament

#Thermodynamic Data
R = 8.314 # pa.m^3/(mol.k)
Tc = 514 # Critical temperature of EtOH in Kelvins
Pc = 6140000 # Critical pressure in pascals
omega = 0.635 #The acentric factor
kappa = 0.37464+1.54226*omega-0.26992*omega**2 # used in the calculation of a
Tlow = 0.3*514 # starting temperature in kelvins
T = Tlow - dT # subtract dT because we will increment by dT when we enter the big loop
x = np.array([2, 3, 1, 0])
n = 39000
m = 5
matrix = [0] * n
for j in range(n):
    matrix[j] = [0] * m
#matrix = np.zeros(39000, 5) # col. 1,2,3,4 and 5 are T,p,vapor molar volume,liquid molar volume and heat of vap., respectively.
# T=0.3*Tc;

def der(T,h): # da/dT used for enthalpy calculations
#==============afunction=======================
    def a(T):
        Tc=514 # Kelvin
        R=8.314
        Pc=6140000 # pascals
        omega=0.635 # unitless
        Tr=T/Tc
        kappa=0.37464+1.54226*omega-0.26992*omega**2
        alpha=(1+kappa*(1-Tr**0.5))**2;
        a=alpha*0.45724*1/Pc*(R*Tc)**2
        return a
# Calculates the derivate of a function numerically based on the 5-pt Stencil
    return (-a(T+2*h)+8*a(T+h)-8*a(T-h)+a(T-2*h))/(12*h);

yaa=der(250, 0.0001)
print(yaa)
caa= matrix[0][2]
print(caa)



i=0; #initialize a counter for the entries in the matrix
pdb.set_trace()
while T < (Tc-0.22): # Big while loop(at each pass the temperature is incremented)
    T = T + dT # incrementing temperature as discussed before
    matrix[i][0] = T # setting the elements of the first column equal to the incremented temperature
#    T=matrix(i,1) # use T instead of the previos line for the sake of ease while reading the code
    if i == 0: # only the first pass we give a guessing for the pressure
        matrix[i][1] = math.exp(-5) # Pa
    else: # for other passes the pressure is set equal to the previous pressure as an initial guess
        matrix[i][1] = matrix[i-1][1]
#==============================================================================
    p = matrix[i][1] # for the ease of reading the code(same as using T-line 23-)
    Tr = T/Tc #reduced temperature
    alpha = (1+kappa*(1-Tr**0.5))**2 # modifying factor for the interactions parameter at different temperatures
    b = 0.07780*R*Tc/Pc
    a=alpha*0.45724*1/Pc*(R*Tc)**2

    check = 1
    pnew = p

    while check > math.exp(-2): # this function does not need a singularity check because it is a continuous cubic function
        A = a*pnew/(R*T)**2
        B = b*pnew/(R*T)
        #=============Cubic equation of state coeffiecients===================#
        c = [1, -1+B, -3*B**2-2*B, -A*B+B**2+B**3]
        #c(1) = 1
        #c(2) = -(1-B)
        #c(3) = A-3*B**2-2*B
        #c(4) = -(A*B-B**2-B**3)
        zroots = np.roots(c) #using the built-in function in matlab to the find the roots of eos
        # The following part is used to take the real roots of eos only
        #index = zroots.find(np.imag(zroots)==0)
        index = np.where(zroots.imag == 0)
        if len(index)>1:
            z = zroots[index]
        else:
            ztemp = zroots[index].real
            #lengg=len(ztemp)
            #z = int(ztemp)
            pdb.set_trace()
            ztemp = ztemp.astype(np.int64)
            z = np.asscalar(ztemp)
            #z = np.take(ztemp,0)
        #taking real roots only ends here
        vliquid = min(z)*R*T/pnew
        vvapor = max(z)*R*T/pnew
        fugv = pnew * math.exp(max(z) - 1 - math.log10((max(z)-B))-A/(2 * math.sqrt(2)*B)*math.log10(((max(z)+(1+math.sqrt(2))*B)/(max(z)+(1-math.sqrt(2))*B)))) # z(1) is the bigger root
        fugl = pnew * math.exp(min(z) - 1 - math.log10((min(z)-B))-A/(2 * math.sqrt(2)*B)*math.log10(((min(z)+(1+math.sqrt(2))*B)/(min(z)+(1-math.sqrt(2))*B)))) # z(2)is the smaller root
        pnew = pnew-R*T*(math.log(fugv/fugl))/(vvapor-vliquid); #newton raphson to find a better estimate of vapor pressure
        check = abs(1-fugv/fugl)
#The end of the while loop that breaks until fugacities are the same for both phases at each T
#This is the criteria for equilibrium in thermodynamics


    h = 0.001

    daoverdt = der(T,h)
    #print(daoverdt)

    Hv=R*T*(max(z)-1)+(T*daoverdt-a)/(2*math.sqrt(2)*b)*math.log10((max(z)+2.44*B)/(max(z)-0.414*B)) # enthalpy of the vapor phase at equilibriem in joules
    Hl=R*T*(min(z)-1)+(T*daoverdt-a)/(2*math.sqrt(2)*b)*math.log10((min(z)+2.44*B)/(min(z)-0.414*B))# enthalpy of the liquid phase at equilibriem in joules


    matrix[i][1]=pnew
    matrix[i][2]=max(z)*R*T/pnew # zrt/p vapor molar volumein m^3
    matrix[i][3]=min(z)*R*T/pnew # liquid molar volume in m^3
    matrix[i][4]= Hv-Hl  # hvap in Joules
    i=i+1
    pdb.set_trace()
print(matrix)
# modifications used for the second plot to reach the peak and make the
# lines touchhing. note that the numbers used are completey random!
matrix[i-1][1] = matrix[i-2][1]+1000
matrix[i-1][2] = matrix[i-2][2]-2.2*math.exp(-5)
matrix[i-1][3] = matrix[i-2][3]+2*math.exp(-5)



# emptying the rest of the matrix after the last temperature increment
matrix[i:end]=[]
# declaring the arrays used in the plotting function
Temperature=matrix[:][1]
Ps=matrix[:][2]
hvap=matrix[:][5]
vvs=matrix[:][3]
vls=matrix[:][4]



# Plotting
fh=figure(1)
set(fh, 'color','w')

#colordef white
subplot(2,1,1)
[ax]=plotyy(Temperature,Ps/1000,Temperature,hvap/1000)

#set(get(ax(1),'Ylabel'),'String','Pressure (kPa)','fontsize',12,'fontangle','oblique','fontweight','bold')
#set(get(ax(2),'Ylabel'),'String','\Delta H ^v (kJ/mol)','fontsize',12,'fontangle','oblique','fontweight','bold','Rotation',-90,'units','normalized','position',[1.1 0.5])

xlabel('Temperature (K)','fontsize',13,'fontangle','normal','fontweight','bold')

title('Ethanol Vapor Pressure & Enthalpy of Vaporization Vs. Temperature','fontsize',15,'fontangle','normal')

hlegend=legend('Vapor Pressure','Enthalpy of Vaporization');

#set(hlegend,'fontsize',12,'box','off','units','normalized','position',[0.37 0.72 0.001 0.001],'fontangle','normal','orientation','vertical')

#set(ax,'box','off','fontsize',13,'xlim',[150 520],'TickDir','out','XTick',150:(520-150)/10:550,'fontangle','normal');

#grid on

#set(ax(1),'ylim',[0 7e3],'YTick',0:7e3/5:7e3)
#set(ax(2),'ylim',[0 6e1],'YTick',0:6e1/5:6e1)
###########################################
subplot(2,1,2)

semilogx(vls,Ps/1000,vvs,Ps/1000)
#grid on
#set(gca,'box','off','TickDir','out','fontsize',13,'ylim',[0 7e3])

title('Pressure Vs. Molar Volume for Ethanol','fontsize',15,'fontangle','normal')
xlabel('Molar Volume (m^3/mol)','fontsize',13,'fontangle','normal','fontweight','bold')
ylabel('Pressure (kPa)','fontsize',13,'fontangle','normal','fontweight','bold')

hlegend=legend('liquid','Vapor')

#set(hlegend,'fontsize',13,'box','off','units','normalized','position',[0.53 0.4 0.001 0.001],'fontangle','normal','orientation','horizontal')

#print(figure(1),'-dpng','-r600','picture1')
