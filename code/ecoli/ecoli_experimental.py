import numpy as np
import matplotlib.pyplot as plot
from math import floor as floor

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

ecoli_count = Ecoli(5)
glucose_count = Glucose(1000)
galactose_count = Galactose(500)

dedt = np.full(201, ecoli_count.ecoli_count)
dgdt = np.full(201, glucose_count.glucose_count)
dgadt = np.full(201, galactose_count.galactose_count)

def normalize_to_zero(value):
    return value if value > 0 else 0


def death_by_starvation(indx):
    ecoli_count.ecolicount = dedt[indx] =  normalize_to_zero(ecoli_count.ecoli_count/2)


def get_indx_to_multiply(baseIndx):
    if baseIndx <= 20:
        return baseIndx
    else:
        return baseIndx - 20


    # Math equation to remember - ex: 101 - isEvent = true; 0 + i%10 would go back to 1-20 series and 10 + i %10 would be for odd values
    # isEven = True if floor(i/10) % 2 == 0 else False
    # return 0 + i % 10 if isEvent else 10 + i % 10


def check_for_sufficient_food(i):
    return glucose_count.glucose_count + galactose_count.galactose_count - ecoli_count.ecoli_count * 0.12 < 0


def partial_multiplication_and_starvation(i):
    multiplying_indx = get_indx_to_multiply(i)
    cell_to_divide = floor(glucose_count.glucose_count + galactose_count.galactose_count / 0.12)
    cells_under_starvation = abs(cell_to_divide - dedt[multiplying_indx])
    dedt[i] = normalize_to_zero(cell_to_divide * 2 - cells_under_starvation)


def ecoli_multiplication(i):
    if glucose_count.glucose_count <= 0 and galactose_count.galactose_count <= 0:
        death_by_starvation(i)

    elif check_for_sufficient_food(i) and dedt[i] > 0:
        partial_multiplication_and_starvation(i)

   # elif i % 2 == 0:
   #     dedt[i] = 2 * ecoli_count.ecoli_count
   #     ecoli_count.ecolicount = dedt[i]

    else:
        multiplying_indx = get_indx_to_multiply(i)
        #dedt[i] = dedt[multiplying_indx] * 2
        dedt[i] = dedt[i - 1] + dedt[i] if i !=0 else dedt[i]
        ecoli_count.ecolicount = dedt[i]


def glucose_consumption(i):
    if i % 30 == 0:
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
