# Modelling-Simulation
# Modelling and Simulation - University of Groningen

## Simulation of Escherichia coli (E. Coli) Growth in Presence of VariousCarbon Sources
### Group - 19
### By, Swastik S Nayak (S4151968) & Siddharth Baskaran (S3922782).

**ecoli_simulation.py**  
  Simulates a theoritcal model of the Ecoli division and Glucose & Galactose consumption.
  The Ecoli cell division duration is downscaled from 20min to 2min to obtain a wider range on the graph.  
  A detailed report of the concept, modelling and simulation can be found in our **[report](https://github.com/Swastik-RUG/Modelling-Simulation/blob/dev/report/E_Coli_Colony_simulation_team_19_report.pdf)**  

### Simulations:
 
## Static Snapshots of the Simulation can be found in the [images/ folder](https://github.com/Swastik-RUG/Modelling-Simulation/tree/dev/images)  
 
## Ecoli Division under Normal scenario, without replenishment of food  
  ![](simulations/Ecoli-Normal.gif)  

## Ecoli Division under Starvation scenario, without replenishment of food  
 ![](simulations/Ecoli-Starvation.gif)  

## Ecoli Division under Lactose Only scenario, without replenishment of food  
   ![](simulations/Ecoli-LactoseOnly.gif)  

 ## Ecoli Division under Stable scenario, with replenishment of food  
  ![](simulations/Ecoli-stable.gif)  

 ## Ecoli Division under Stable Glucose and Galactose only scenario, with replenishment of food  
  ![](simulations/Ecoli-stableGlucose.gif)  

 ## Ecoli Division under Stable Lactose scenario, with replenishment of food  
 ![](simulations/Ecoli-stableLactose.gif)  

## Create you Custom Simulations by following these steps

- Clone this repository
- Find the **"conf.ini"** file.
  ![](images/conf_ini.PNG)
- Execute the **"ecoli_simulation.py"** and select **simulation mode "custom"**
  ![](images/CustomSimulationRun.PNG)
- Done

## Conf.ini meaning of independent variables

lactose count: An integer value, which can be used to set the initial lactose molecule count in the system.
- **glucose count:** An integer value, which can be used to set the initial glucose molecule count in the system.
- **galactose count:** An integer value, which can be used to set the initial galactose molecule count in the system.
- **lactase:** An integer value, which can be used to set the initial glucose lactase or beta-galactosidase count in the
system.
- **ecoli:** An integer value, which can be used to set the initial ecoli count in the system.
- **iterations:** Number of iterations the simulation has to be performed.
- **replenish interval:** If the experimenter wants to simulate an environment where the carbon sources are replenished for specific time-frames. If this value is set to 0, the simulation will not replenish any resources, if it
is set to a R, where R ¿ 0; the simulation will replenish the carbon sources in the intervals R,2R,3R.....
- **replenish glucose:** An integer value, which can be set to indicate the glucose count that has to be replenished
at replenish intervals.
- **replenish galactose:** An integer value, which can be set to indicate the galactose count that has to be replenished at replenish intervals.
- **replenish lactose:** An integer value, which can be set to indicate the lactose count that has to be replenished
at replenish intervals.

For any Assistance regarding this project, kindly create an issue and we will resolve the issue at our earliest possible.
