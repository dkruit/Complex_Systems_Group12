# Complex Systems Group 12 Project

## Research Question
In 2017 Tesla built the largest battery in the world in South Australia [1]. The goal of this 100 MW battery is to reduce power blackouts and stabilize the grid by supplying or absorb burst of energy very quickly. The magnitude of blackouts is known to have a power law tail [2]. We want to study how a large battery affects the critcial load of the system and the exponent of the power law.

## Method
We will implement the ORNL-PSerc-Alaska (OPL) model [3,4]. This model was designed to study the complex behavior of the dynamics of series of blackouts. The researches we demonstrated self-organization of the system to a critical point at which the probability distribution of blackout size resembles real world data from the North American Electric Reliability Corporation.

We first aim to implement the model and replicate results in papers [4] and [5]. Then we want to add a large stabilizing battery to the network and compare the network with and without battery.

## Running the simulation
To run the simulation, set the parameters as desired in 'Network.py'.
Important note: in order to reduce the running time, the default parameters are now set to running 2 experiments for 500 days, instead of 10 experiments for 1000 days.
The data and figures will be saved into a map with name 'EXPERIMENT_NAME'.
The following data will be saved per simulation i:
- flines_simulation_i: the total number of lines that failed per event;
- mean_M_linei: the fractional overload M per line per day;
- mean_Mi: the average fractional overload M of all lines per day;
- N_load_shed_i: the fractional load shedded per event.

If PLOT_FIGS = True, the following figures will be saved per simulation i:
- frequency_failedlines_hist_i: a histogram of the number of failed lines per event on log linear scale
- frequency_loadshed_log_i: the frequency of the shedded load on a log log scale
- Nrevents_failedlines_i: nr of events with nr of failed lines on a linear scale
- probability_failedlines_log_i: probability density function of the number of failed lines per event on a log log scale

These plots are all plots without error bars. In order to create plots with error bars, fun the file 'fline_plot.py'.

If NETWORK_FIGS = True, the following figures will be saved per simulation day i (and overwritten every new simulation, since this function is only for visualization):
- day_i_a_updated_p: the network after updating p;
- day_i_b_initial_failure: the network after the initial failures;
- day_i_c_after_redispatch: the network after each redispatchment loop.


## References
1. BBC.com, "Tesla mega-battery in Australia activated", 1 December 2017. Url: [link](https://www.bbc.com/news/world-australia-42190358 "https://www.bbc.com/news/world-australia-42190358")
2. Carreras, B. A., et al. "Evidence for self-organized criticality in electric power system blackouts." Proceedings of the 34th Annual Hawaii International Conference on System Sciences. IEEE, 2001.
3. Dobson, Ian, et al. "An initial model for complex dynamics in electric power system blackouts." Proceedings of the 34th annual Hawaii international conference on system sciences. Vol. 2. 2001.
4. Carreras, Benjamin A., et al. "Dynamics, criticality and self-organization in a model for blackouts in power transmission systems." Proceedings of the 35th Annual Hawaii International Conference on System Sciences. IEEE, 2002.
5. Dobson, Ian, et al. "An initial complex systems analysis of the risks of blackouts in power transmission systems." Power Systems and Communications Infrastructures for the future, Beijing, China (2002).
