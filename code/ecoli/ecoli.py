import numpy as np
import matplotlib.pyplot as plot

class Ecoli:
    def __init__(self, count):
        self.ecolicount = count

    @property
    def ecoli_count(self):
        return self.ecolicount


class Glucose:
    def __init__(self, count):
        self.glucosecount = count

    @property
    def glucose_count(self):
        return self.glucosecount

class Galactose:
    def __init__(self, count):
        self.galactosecount = count

    @property
    def galactose_count(self):
        return self.galactosecount


seed_input=[500, 1000]
time = np.linspace(0, 200, 200)
#ode_res = odeint(simulate_ecoli_glucose_consumption, seed_input, time)

ecoli_count = Ecoli(50)
glucose_count = Glucose(1000)
galactose_count = Galactose(500)

dedt = np.empty(201)
dgdt = np.empty(201)
dgadt = np.empty(201)

def death_by_starvation(indx):
    ecoli_count.ecolicount = dedt[indx] = ecoli_count.ecoli_count/2


def ecoli_multiplication(i):
    if glucose_count.glucose_count <= 0 and galactose_count.galactose_count <= 0:
        death_by_starvation(i)
    elif i == 0:
        dedt[i] = ecoli_count.ecolicount
        dgadt[i] = galactose_count.galactose_count
    elif i % 2 == 0:
        dedt[i] = 2 * ecoli_count.ecoli_count
        ecoli_count.ecolicount = dedt[i]
    else:
        dedt[i] = ecoli_count.ecoli_count


def glucose_consumption(i):
    if i % 3 == 0:
        glucose_count.glucosecount = glucose_count.glucosecount + 100

    if glucose_count.glucose_count > 0:
        dgdt[i] = glucose_count.glucose_count - ecoli_count.ecoli_count * 0.12
        glucose_count.glucosecount = dgdt[i] if dgdt[i] >= 0 else 0

        if dgdt[i] < 0:
            galactose_consumption(i, dgdt[i])
            glucose_count.glucosecount = 0
            dgdt[i] = 0
        else:
            glucose_count.glucosecount = dgdt[i]

        dgadt[i] = galactose_count.galactose_count
    else:
        dgdt[i] = 0
        galactose_consumption(i)
        glucose_count.glucosecount = 0



def galactose_consumption(indx, difference=0):
    dgadt[indx] = galactose_count.galactose_count - ecoli_count.ecolicount * 0.1 + difference
    galactose_count.galactosecount = dgadt[indx] = dgadt[indx] if dgadt[indx] >= 0 else 0

for i in range(201):
    plot.clf()
    glucose_consumption(i)
    ecoli_multiplication(i)

    if i == 200:

        plot.plot(time, dgdt[0:200], 'r', label='Glucose-Consumption')
        plot.plot(time, dedt[0:200], 'g', label='Ecoli-Multiplication')
        plot.plot(time, dgadt[0:200], 'y', label='GalactoseConsumption')
        plot.legend()
        plot.show()
