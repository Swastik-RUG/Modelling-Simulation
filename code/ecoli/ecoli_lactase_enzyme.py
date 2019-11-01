from math import floor as floor

import matplotlib.animation as animation
import matplotlib.pyplot as plot
import matplotlib.style as style
import numpy as np
from configparser import ConfigParser
from tkinter import *


class Ecoli:
    def __init__(self, count):
        self.ecolicount = count

    # self.betaGalactosidase = 10
    # self.starvationperiod = 0

    @property
    def ecoli_count(self):
        return self.ecolicount

    # @property
    # def getBetaGalctosidase(self):
    #     return self.betaGalactosidase
    #
    # @property
    # def getStarvationPeriod(self):
    #     return self.starvationperiod


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


class Lactose:
    def __init__(self, count):
        self.lactosecount = count

    @property
    def lactose_count(self):
        return self.lactosecount


class LactaseEnzyme:
    def __init__(self, count):
        self.lactase = count

    @property
    def lactose_count(self):
        return self.lactase


parser = ConfigParser()
parser.read("conf.ini")

master = Tk()
var = StringVar()
Label(master, text="Select OCR language").grid(row=0, sticky=W)
Radiobutton(master, text="Normal", variable=var, value='Normal').grid(row=1, sticky=W)
Radiobutton(master, text="OnlyLactose", variable=var, value='OnlyLactose').grid(row=2, sticky=W)
Radiobutton(master, text="Starvation", variable=var, value='Starvation').grid(row=3, sticky=W)
Radiobutton(master, text="stable_colony", variable=var, value='stable_colony').grid(row=4, sticky=W)
Radiobutton(master, text="custom", variable=var, value='custom').grid(row=5, sticky=W)
Button(master, text="OK", command=master.quit).grid(row=6, sticky=W)
mainloop()


def getSelectionType(selection):
    if selection == 1:
        return "Normal"
    elif selection == 2:
        return "OnlyLactose"
    elif selection == 3:
        return "Starvation"
    elif selection == 4:
        return "custom"
    else:
        return "Normal"


selection = var.get()

selection2 = getSelectionType(parser.get(selection, "ecoli"))

itr = parser.getint(selection, "iterations")
# seed_input = [500, 1000]
time = np.linspace(0, itr + 1, itr + 1)

ecoli_count = Ecoli(parser.getint(selection, "ecoli"))
glucose_count = Glucose(parser.getint(selection, "glucose_count"))
galactose_count = Galactose(parser.getint(selection, "galactose_count"))
lactose_count = Lactose(parser.getint(selection, "lactose_count"))
beta_galactosidase = LactaseEnzyme(parser.getint(selection, "lactase"))
replenish_glucose = parser.getint(selection, "replenish_glucose") if (parser.has_option(selection, "replenish_glucose")) else 0
replenish_galactose = parser.getint(selection, "replenish_galactose") if (parser.has_option(selection, "replenish_galactose")) else 0
replenish_lactose = parser.getint(selection, "replenish_lactose") if (parser.has_option(selection, "replenish_lactose")) else 0
replenish_interval = parser.getint(selection, "replenish_interval") if (parser.has_option(selection, "replenish_interval")) else 0
sample_frame = 2


# window = Tk()
# User_input = Entry()
# User_input.pack()
# #\user_problem  = int(User_input.get())
# window.mainloop()
#
# user_problem  = int(User_input.get())
# print(user_problem)

# itr = 50
# seed_input = [500, 1000]
# time = np.linspace(0, itr + 1, itr + 1)
#
# ecoli_count = Ecoli(1)
# glucose_count = Glucose(100000)
# galactose_count = Galactose(1000)
# lactose_count = Lactose(100000)
# beta_galactosidase = LactaseEnzyme(10000)
# sample_frame = 2


# ecoli_colony = []

# for i in range(sample_frame):
#    ecoli_colony.append(ecoli_count)

def create_seed_array(initial_seed, dummy, frame, extra=1):
    return np.concatenate([np.full(frame, initial_seed), np.full(itr + extra - sample_frame, dummy)], axis=0)


dedt = create_seed_array(ecoli_count.ecoli_count, 0, sample_frame)
dgdt = create_seed_array(glucose_count.glucose_count, 0, 1, 2)
dgadt = create_seed_array(galactose_count.galactose_count, 0, 1, 2)
dldt = create_seed_array(lactose_count.lactose_count, 0, 1, 2)
lactaseDt = create_seed_array(10 * ecoli_count.ecoli_count, 0, 1, 2)
resource_consumption = 0.12
lactase_enzyme_depletion_rate = 0.1


def normalize_to_zero(value):
    return value if value > 0 else 0


def death_by_starvation(indx):
    ecoli_count.ecolicount = dedt[indx] = normalize_to_zero(ecoli_count.ecoli_count / 2)


def get_indx_to_multiply(baseIndx):
    if baseIndx < sample_frame:
        return baseIndx
    else:
        return baseIndx - sample_frame


def get_available_glucose_and_galactose():
    return glucose_count.glucose_count + galactose_count.galactose_count


def partial_multiplication_and_starvation(indx):
    multiplying_indx = get_indx_to_multiply(indx)
    cell_to_divide = floor(glucose_count.glucose_count + galactose_count.galactose_count / 0.12)
    cells_under_starvation = abs(cell_to_divide - dedt[multiplying_indx])
    dedt[indx] = normalize_to_zero(cell_to_divide * 2 - cells_under_starvation)


def check_if_colony_is_dead(itr):
    return ecoli_count.ecolicount <= 0


def check_for_lactose_availability(itr):
    return lactose_count.lactose_count > 0


def calc_possible_lactose_metabolism(cellColony, lactoseAvailable):
    possible_metabolism = floor(beta_galactosidase.lactase / 100)
    return possible_metabolism if (lactoseAvailable >= possible_metabolism) else lactoseAvailable


def multiply_under_adundant_resources(expected_resource_requirement,
                                      expected_ecoli_division,
                                      curr_ecoli_count,
                                      lactose_swap_period,
                                      possible_lactose_metabolism,
                                      itr):
    glucose_consumed = glucose_count.glucose_count - expected_resource_requirement
    deduct_from_galactose = (galactose_count.galactose_count - glucose_consumed) if glucose_consumed < 0 else 0
    glucose_count.glucosecount = glucose_consumed if glucose_consumed >= 0 else 0
    galactose_count.galactosecount = (
            galactose_count.galactose_count - deduct_from_galactose) if glucose_consumed >= 0 else 0

    # lactose_metabolism = possible_lactose_metabolism if deduct_from_galactose < 0 and lactose_swap_period > 2 else 0

    ecoli_count.ecolicount = expected_ecoli_division
    dedt[itr] = ecoli_count.ecoli_count
    # dedt[itr].starvationperiod = 0

    dgdt[itr] = glucose_count.glucose_count
    dgadt[itr] = galactose_count.galactose_count
    dldt[itr] = lactose_count.lactose_count

    curr_lactase_count = lactaseDt[itr:itr + 1]
    # lactase_from_new_cells = (expected_ecoli_division - curr_ecoli_count) * 10
    new_lactase = curr_lactase_count - curr_lactase_count * lactase_enzyme_depletion_rate  # + lactase_from_new_cells
    lactaseDt[
        itr] = new_lactase if new_lactase >= ecoli_count.ecoli_count * 10 else curr_lactase_count  # + lactase_from_new_cells


starvationperiod = 0


def consume_resources_and_multiply(itr):
    lactose_swap_period = 0
    colony_status = check_if_colony_is_dead(itr)
    if colony_status:
        print("The E-coli colony is dead")
    #    elif resource_scarcity <= 0:
    #        lactose_availability = check_for_lactose_availability(itr)
    #        partial_multiplication_and_starvation(itr)

    else:
        multiplying_indx = get_indx_to_multiply(itr)
        curr_ecoli_count = dedt[multiplying_indx]

        colony_slice = dedt[itr - sample_frame:itr] if itr >= sample_frame else dedt[0:itr + 1]
        # curr_total_ecoli_count = sum(e.ecoli_count for e in (ecoli_colony[i-20:i] if i >= 20 else ecoli_colony[0:i]))
        expected_ecoli_division = curr_ecoli_count * 2
        expected_resource_requirement = curr_ecoli_count * resource_consumption
        available_glucose_galactose = get_available_glucose_and_galactose()
        possible_lactose_metabolism = calc_possible_lactose_metabolism(colony_slice, lactose_count.lactose_count)
        total_resources_available = available_glucose_galactose + possible_lactose_metabolism
        insufficient_raw_resources = True if available_glucose_galactose - expected_resource_requirement < 0 else False
        resource_state = total_resources_available - expected_resource_requirement
        # curr_lactase_count = sum(l for l in (lactaseDt[itr - 20:itr] if itr >= 20 else lactaseDt[0:itr + 1]))
        lactaseDt[itr] = beta_galactosidase.lactase

        lactose_swap_period = 0 if insufficient_raw_resources else lactose_swap_period + 1
        if lactose_swap_period > 2:
            glucose_count.glucosecount += possible_lactose_metabolism
            galactose_count.galactosecount += possible_lactose_metabolism
            lactose_count.lactosecount -= possible_lactose_metabolism
            beta_galactosidase.lactase = beta_galactosidase.lactase + ecoli_count.ecoli_count * 50
            lactaseDt[itr] = beta_galactosidase.lactase

        if resource_state >= 0:

            multiply_under_adundant_resources(expected_resource_requirement,
                                              expected_ecoli_division,
                                              curr_ecoli_count,
                                              lactose_swap_period,
                                              possible_lactose_metabolism,
                                              itr)

            # Colony starvation - Assuming colony survives for 2 iterations and kills 30% of the population
        # elif sufficient_raw_resources <= 0 < possible_lactose_metabolism:

        elif resource_state <= 0:
            starved_cells = colony_slice
            # Part of cells that does not divide
            available_raw_resource = get_available_glucose_and_galactose()
            cells_under_starvation = expected_ecoli_division - (available_raw_resource / resource_consumption)
            cells_with_sufficient_resources = expected_ecoli_division - cells_under_starvation
            # bug
            dedt[itr] = ecoli_count.ecoli_count + cells_with_sufficient_resources
            dedt[itr] = (dedt[itr] - dedt[itr] * 0.3)  # if (starvationperiod > 2) else dedt[itr]
            # dedt[itr].starvationperiod += 1 if (dedt[itr].starvationperiod > 2) else -dedt[itr].starvationperiod
            # starvationperiod += -starvationperiod if (starvationperiod > 2) else 1
            ecoli_count.ecolicount = dedt[itr]
            dgdt[itr] = 0  # glucose_count.glucose_count
            glucose_count.glucosecount = 0
            dgadt[itr] = 0  # galactose_count.galactose_count
            galactose_count.galactosecount = 0
            lactose_count.lactosecount = dldt[itr] = lactose_count.lactose_count - possible_lactose_metabolism
            lactaseDt[itr] = beta_galactosidase.lactase
        else:
            dedt[itr] = ecoli_count.ecoli_count
            dgdt[itr] = glucose_count.glucose_count
            dgadt[itr] = galactose_count.galactose_count
            dldt[itr] = lactose_count.lactose_count
            lactaseDt[itr] = beta_galactosidase.lactase


for i in range(itr + 1):
    if replenish_interval != 0 and itr % replenish_interval == 0:
        glucose_count.glucosecount += replenish_glucose
        galactose_count.galactosecount += replenish_galactose
        lactose_count.lactosecount += replenish_lactose

    consume_resources_and_multiply(i)

style.use("fivethirtyeight")
fig = plot.figure("Ecoli bacterial growth simulation:" + selection)
ax = fig.add_subplot(1, 1, 1)


def displayPlot(i):
    #    ax1.clear()
    ax.cla()
    ax.title.set_text('Ecoli bacterial growth simulation:' + selection)
    if parser.getboolean("display_control", "plot_glucose"):
        ax.plot(time[0:i], dgdt[0:i], 'r', label='Glucose-Consumption')
    if parser.getboolean("display_control", "plot_ecoli"):
        ax.plot(time[0:i], dedt[0:i], 'g', label='Ecoli-Multiplication')
    if parser.getboolean("display_control", "plot_galactose"):
        ax.plot(time[0:i], dgadt[0:i], 'y', label='GalactoseConsumption')
    if parser.getboolean("display_control", "plot_lactase"):
        ax.plot(time[0:i], lactaseDt[0:i], 'b', label='Lactase enzyme')
    ax.legend()
    ax.set_xlabel('Time in min')
    ax.set_ylabel('Ecoli Growth')


#  ax1.plot(time[0:i], dldt[0:i], 'c', label='LactoseConsumption')


ani = animation.FuncAnimation(fig, displayPlot, interval=500)
plot.tight_layout()
plot.show()
