from math import floor as floor

import matplotlib.animation as animation
import matplotlib.pyplot as plot
import matplotlib.style as style
import numpy as np
from configparser import ConfigParser
from tkinter import *

# ----------------------------------------------------------------
# Declaration of classes used in the Simulation
# Classes were required as this was not a simple model to solve
# ----------------------------------------------------------------
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


class LactaseSwapPeriod:
    def __init__(self, count):
        self.lactaseSwap = count

    @property
    def lactase_swap_period(self):
        return self.lactaseSwap

# Create a parser object to read and read the conf.ini file
parser = ConfigParser()
parser.read("conf.ini")

# --------------------------------------------------------------------------
# Generate the option selection window and take the input from the user
# --------------------------------------------------------------------------
master = Tk()
var = StringVar()
Label(master, text="Select Simulation to run").grid(row=0, sticky=W)
Radiobutton(master, text="Normal", variable=var, value='Normal').grid(row=1, sticky=W)
Radiobutton(master, text="OnlyLactose", variable=var, value='OnlyLactose').grid(row=2, sticky=W)
Radiobutton(master, text="Starvation", variable=var, value='Starvation').grid(row=3, sticky=W)
Radiobutton(master, text="stable_colony", variable=var, value='stable_colony').grid(row=4, sticky=W)
Radiobutton(master, text="stable_colony_glucose", variable=var, value='stable_colony_glucose').grid(row=5, sticky=W)
Radiobutton(master, text="stable_colony_lactose", variable=var, value='stable_colony_lactose').grid(row=6, sticky=W)
Radiobutton(master, text="custom", variable=var, value='custom').grid(row=7, sticky=W)
Button(master, text="OK", command=master.quit).grid(row=8, sticky=W)
mainloop()

# ------------------------------------------------------------------------------------------------------------------
# Extract the Independent variable parameters from the Configuration file (conf.ini) based on the selection
# ------------------------------------------------------------------------------------------------------------------
selection = var.get()

# Set number of iterations
itr = parser.getint(selection, "iterations")
# Generate the time axis
time = np.linspace(0, itr + 1, itr + 1)
# Set the starting ecoli count
ecoli_count = Ecoli(parser.getint(selection, "ecoli"))
# Set the starting glucose count
glucose_count = Glucose(parser.getint(selection, "glucose_count"))
# Set the starting galactose count
galactose_count = Galactose(parser.getint(selection, "galactose_count"))
# Set the starting Lactose count
lactose_count = Lactose(parser.getint(selection, "lactose_count"))
# Set the starting Lactase count
beta_galactosidase = LactaseEnzyme(parser.getint(selection, "lactase"))
# --------------------------------------------------------------------------------------
# Set the replenishment values if they are set
# --------------------------------------------------------------------------------------
replenish_glucose = parser.getint(selection, "replenish_glucose") if (
    parser.has_option(selection, "replenish_glucose")) else 0
replenish_galactose = parser.getint(selection, "replenish_galactose") if (
    parser.has_option(selection, "replenish_galactose")) else 0
replenish_lactose = parser.getint(selection, "replenish_lactose") if (
    parser.has_option(selection, "replenish_lactose")) else 0
replenish_interval = parser.getint(selection, "replenish_interval") if (
    parser.has_option(selection, "replenish_interval")) else 0

# ----------------------------------------------------------------------------------------------------------------------
# To make our model realistic we have divided the ecoli into sample frames
# Since 1 unit of time is 10min are each E.Coli replicate at 2 iterations
# We treat E.coli as a batch of 2 sets ( 1st element of array and 2nd element)
# This was done to simulate an environment where an Ecoli cell always divides at any given time like a biological ecoli.
# ----------------------------------------------------------------------------------------------------------------------
sample_frame = 2
# Initialize a variable to keep track of Lactose swap period; at 2 the colony swaps to lactose consumption
lactose_swap_period = LactaseSwapPeriod(0)


# Create a numpy array, generic method
def create_seed_array(initial_seed, dummy, frame, extra=1):
    return np.concatenate([np.full(frame, initial_seed), np.full(itr + extra - sample_frame, dummy)], axis=0)

# ----------------------------------------------------------------------------------------------------------------------
# Create seed for Ecoli, Glucose, Galactose, Lactose, Lactase
# ----------------------------------------------------------------------------------------------------------------------
dedt = create_seed_array(ecoli_count.ecoli_count, 0, sample_frame)
dgdt = create_seed_array(glucose_count.glucose_count, 0, 1, 2)
dgadt = create_seed_array(galactose_count.galactose_count, 0, 1, 2)
dldt = create_seed_array(lactose_count.lactose_count, 0, 1, 2)
lactaseDt = create_seed_array(10 * ecoli_count.ecoli_count, 0, 1, 2)
# Rate at which each E.Coli consumes the resources
resource_consumption = 0.12
# Rate at which the lactase enzyme deplete when the E.Coli uses only Glucose and Galactose
lactase_enzyme_depletion_rate = 0.1

# ----------------------------------------------------------------------------------------------------------------------
# An internal helper method to fetch the parent of the sample we are replicating
# Ex: Sample_frame = 2;
#     at iteration 3 -> the E.Coli colony to replicate would be 1 (itr - sample_frame)
#     at iteration 4 -> The E.Coli colony to replicate would be 2
#     at iteration 5 -> Parent would be iteration 3 output.
# ----------------------------------------------------------------------------------------------------------------------
def get_indx_to_multiply(baseIndx):
    if baseIndx < sample_frame:
        return baseIndx
    else:
        return baseIndx - sample_frame


# ----------------------------------------------------------------------------------------------------------------------
# Get the available Glucose and Galactose available in the system
# ----------------------------------------------------------------------------------------------------------------------
def get_available_glucose_and_galactose():
    return glucose_count.glucose_count + galactose_count.galactose_count


# ----------------------------------------------------------------------------------------------------------------------
# Checks if the E.Coli colony is extinct
# ----------------------------------------------------------------------------------------------------------------------
def check_if_colony_is_dead():
    return ecoli_count.ecolicount <= 0


# ----------------------------------------------------------------------------------------------------------------------
# Calculate the possible lactose metabolism using Lactase count
# ----------------------------------------------------------------------------------------------------------------------
def calc_possible_lactose_metabolism(cellColony, lactoseAvailable):
    possible_metabolism = floor(beta_galactosidase.lactase / 100)
    return possible_metabolism if (lactoseAvailable >= possible_metabolism) else lactoseAvailable


# ----------------------------------------------------------------------------------------------------------------------
# Handles the multiplication and resource consumption of the colony when there is sufficient resources
# ----------------------------------------------------------------------------------------------------------------------
def multiply_under_adundant_resources(expected_resource_requirement,
                                      expected_ecoli_division,
                                      itr):
    # Deduct the require resource count from Glucose and Galactose
    glucose_consumed = glucose_count.glucose_count - expected_resource_requirement
    deduct_from_galactose = (galactose_count.galactose_count - glucose_consumed) if glucose_consumed < 0 else 0
    glucose_count.glucosecount = glucose_consumed if glucose_consumed >= 0 else 0
    galactose_count.galactosecount = (galactose_count.galactose_count - deduct_from_galactose) if glucose_consumed >= 0 else 0

    # Multiply the E.coli * 2 as there is abundant resources
    ecoli_count.ecolicount = expected_ecoli_division
    # Store the E.Coli colony status to the numpyArray to plot later
    dedt[itr] = ecoli_count.ecoli_count

    # Record the status of all the resources for plotting
    dgdt[itr] = glucose_count.glucose_count
    dgadt[itr] = galactose_count.galactose_count
    dldt[itr] = lactose_count.lactose_count

    # For every new E.Coli in the system we add 10 lactase enzymes so as to have a minimum number of lactase in the system.
    # This is done to prevent the model from being punished by lack of Lactase enzyme when the colony grows massively
    curr_lactase_count = lactaseDt[itr:itr + 1]
    new_lactase = curr_lactase_count - curr_lactase_count * lactase_enzyme_depletion_rate
    lactaseDt[itr] = new_lactase if new_lactase >= ecoli_count.ecoli_count * 10 else curr_lactase_count


# ----------------------------------------------------------------------------------------------------------------------
# A method that Orchestrates the replication and resource utilization of the model
# ----------------------------------------------------------------------------------------------------------------------
def consume_resources_and_multiply(itr):
    colony_status = check_if_colony_is_dead()
    if colony_status:
        # Colony is dead, set the resources to the count it was prior to extinction of the colony
        print("The E-coli colony is dead")
        dldt[itr] = lactose_count.lactose_count
        dgdt[itr] = glucose_count.glucose_count
        dgadt[itr] = galactose_count.galactose_count
        lactaseDt[itr] = beta_galactosidase.lactase
    else:
        # Get the parent index to multiply
        multiplying_indx = get_indx_to_multiply(itr)
        # Get the E.coli count from the parent index
        curr_ecoli_count = dedt[multiplying_indx]

        # Get a colony slice to compute the Lactase available in the system
        # The colony is bucketed into sample_frame but the environment they live in is the same
        colony_slice = dedt[itr - sample_frame:itr] if itr >= sample_frame else dedt[0:itr + 1]
        # Expected cell division under Glucose environment
        expected_ecoli_division = curr_ecoli_count * 2
        # Required Glucose and Galactose resources for the cell division to occur
        expected_resource_requirement = curr_ecoli_count * resource_consumption
        # Available Glucose and Galactose in the system
        available_glucose_galactose = get_available_glucose_and_galactose()
        # The amount of Lactose that can be metabolised (LM equation in the report)
        possible_lactose_metabolism = calc_possible_lactose_metabolism(colony_slice, lactose_count.lactose_count)
        # Total resources Glucose, Glactose and Lactose combined
        total_resources_available = available_glucose_galactose + possible_lactose_metabolism
        # State of the system, decides if the colony has sufficient resources or not
        resource_state = total_resources_available - expected_resource_requirement
        # Lactase enzyme present in the system, stored for plotting
        lactaseDt[itr] = beta_galactosidase.lactase

        # If the E.Coli colony is in a Lactose consumption phase then metabolize the Lactose to Glucose and Galactose
        if lactose_swap_period.lactase_swap_period > 2:
            glucose_count.glucosecount += possible_lactose_metabolism
            galactose_count.galactosecount += possible_lactose_metabolism
            lactose_count.lactosecount -= possible_lactose_metabolism
            # Reducing the Lactase enzyme
            beta_galactosidase.lactase = (beta_galactosidase.lactase + ecoli_count.ecoli_count * 50) if(lactose_count.lactose_count > 0) else beta_galactosidase.lactase
            # Reset the Lactase enzyme value again for plotting
            lactaseDt[itr] = beta_galactosidase.lactase

        # if Lactose is the only source, decrease the growth by 20%
        expected_ecoli_division -= 0 if available_glucose_galactose > 0 else expected_ecoli_division * 0.20

        # If there is Sufficient resources, including Lactose
        if resource_state >= 0:
            # Check if Glucose and Galactose is sufficient
            # if insufficient, start the Lactose swap counter to 1, else decrease it by 1
            if available_glucose_galactose >= expected_resource_requirement and lactose_swap_period.lactase_swap_period > 0:
                lactose_swap_period.lactaseSwap = lactose_swap_period.lactase_swap_period - 1
            else:
                lactose_swap_period.lactaseSwap = lactose_swap_period.lactase_swap_period + 1

            # Multiply the Colony and Consume the resources
            multiply_under_adundant_resources(expected_resource_requirement,
                                              expected_ecoli_division,
                                              itr)

        # Colony starvation - Assuming colony survives for 2 iterations and kills 30% of the population
        elif resource_state <= 0:
            # Available Glucose and Galctose
            available_raw_resource = get_available_glucose_and_galactose()
            # Part of the Colony under starvation
            cells_under_starvation = expected_ecoli_division - (available_raw_resource / resource_consumption)
            # Part of the Colony who received sufficient resources
            cells_with_sufficient_resources = expected_ecoli_division - cells_under_starvation
            # Divide the cells which had sufficient resources only
            dedt[itr] = ecoli_count.ecoli_count + cells_with_sufficient_resources
            # Kill 30% of the E.Coli
            dedt[itr] = (dedt[itr] - dedt[itr] * 0.3)
            # Set the Ecoli count for plotting
            ecoli_count.ecolicount = dedt[itr]
            # Set other resource values for plotting
            dgdt[itr] = 0
            glucose_count.glucosecount = 0
            dgadt[itr] = 0
            galactose_count.galactosecount = 0
            lactose_count.lactosecount = dldt[itr] = lactose_count.lactose_count - possible_lactose_metabolism
            lactaseDt[itr] = beta_galactosidase.lactase
        else:
            # This else part never occurs, this is just a placeholder.
            # By the time resources are completely empty the colony would be dead
            dedt[itr] = ecoli_count.ecoli_count
            dgdt[itr] = glucose_count.glucose_count
            dgadt[itr] = galactose_count.galactose_count
            dldt[itr] = lactose_count.lactose_count
            lactaseDt[itr] = beta_galactosidase.lactase


# Iterate the model for itr number of times
for i in range(itr + 1):
    # If Replenish interval has been set, replenish the resources at its interval
    if replenish_interval != 0 and itr % replenish_interval == 0:
        glucose_count.glucosecount += replenish_glucose
        galactose_count.galactosecount += replenish_galactose
        lactose_count.lactosecount += replenish_lactose

    # Start the simulation for the iteration
    consume_resources_and_multiply(i)

# Create subplots for the graphs
style.use("fivethirtyeight")
fig = plot.figure("Ecoli bacterial growth simulation:" + selection)
ax = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)


# Helper function used by animation to animate the simulation
def displayPlot(i):
    ax.cla()
    ax2.cla()
    ax3.cla()
    ax4.cla()
    ax.title.set_text('Ecoli bacterial growth simulation:' + selection)
    if parser.getboolean("display_control", "plot_ecoli"):
        ax.plot(time[0:i], dedt[0:i], 'g', label='Ecoli-Multiplication')
    if parser.getboolean("display_control", "plot_glucose"):
        ax2.plot(time[0:i], dgdt[0:i], 'r', label='Glucose-Consumption')
    if parser.getboolean("display_control", "plot_galactose"):
        ax2.plot(time[0:i], dgadt[0:i], 'y', linestyle="-.", label='GalactoseConsumption')
    if parser.getboolean("display_control", "plot_lactase"):
        ax3.plot(time[0:i], lactaseDt[0:i], 'b', label='Lactase enzyme')
    if parser.getboolean("display_control", "plot_lactose"):
        ax4.plot(time[0:i], dldt[0:i], 'c', linestyle="--", label='Lactose Consumption')
    ax.legend()
    ax.set_xlabel('Time in min')
    ax.set_ylabel('Ecoli Growth')
    ax2.legend()
    ax2.set_xlabel('Time in min')
    ax2.set_ylabel('Glucose + Galactose')
    ax3.legend()
    ax3.set_xlabel('Time in min')
    ax3.set_ylabel('Lactase enzyme')
    ax4.legend()
    ax4.set_xlabel('Time in min')
    ax4.set_ylabel('Lactose Consumption')

# Animation support with interval 500, 0.5sec
ani = animation.FuncAnimation(fig, displayPlot, interval=500)
plot.tight_layout()
plot.show()
