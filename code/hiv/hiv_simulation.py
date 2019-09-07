import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def hiv_infection_simulator(cells, time):
    healthy_cells = cells[0]
    infected_cells = cells[1]
    virus_cells = cells[2]

    if time > 1:
        antidote = 10
    else:
        antidote = 0


    kr1 = 1e5  # new healthy cells per year
    kr2 = 0.1  # death rate of healthy cells
    kr3 = 2e-7  # healthy cells converting to infected cells
    kr4 = 0.5  # death rate of infected cells
    kr5 = 5  # death rate of virus
    kr6 = 100  # production of virus by infected cells
    kr7 = 1.0
    healthy_cells_affected =  kr3 * healthy_cells * virus_cells

    # Ordinary differential equations
    # Assuming 100ml antidote kills 80% virus cells but kills the patient at 1000ml dosage
    dhdt = kr1 - kr2 * healthy_cells - healthy_cells_affected
    didt = healthy_cells_affected - kr4 * infected_cells
    dvdt = - healthy_cells_affected - kr5 * virus_cells + kr6 * infected_cells - kr7 * (virus_cells * (antidote / 100 * 0.1))

    return [dhdt, didt, dvdt]


seed_input = [1000, 0, 100]
time = np.linspace(0, 15, 1000)
ode_result = odeint(hiv_infection_simulator, seed_input, time)

healthy_cells = ode_result[:, 0]
infect_cells = ode_result[:, 1]
virus_cells = ode_result[:, 2]

plt.semilogy(time, healthy_cells)
plt.semilogy(time, infect_cells)
plt.semilogy(time, virus_cells)
plt.show()

