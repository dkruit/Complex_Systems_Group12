# Complex Systems Group 12 Project

## Research Question
In 2017 Tesla built the largest battery in the world in South Australia [1]. The goal of this 100 MW battery is to reduce power blackouts and stabilize the grid. The magnitude of blackouts is known to have a power law tail [2]. We want to study how a large battery affects the critcial load of the system and the exponent of the power law.

## Method
We will implement the ORNL-PSerc-Alaska (OPL) model [3,4]. This model was designed to study the complex behavior of the dynamics of series of blackouts. The researches we demonstrated self-organization of the system to a critical point at which the probability distribution of blackout size resembles real world data from the North American Electric Reliability Corporation.

We first aim to implement the model and replicate results in from [4] and [5]. Then we want to add a large stabilizing battery to the network and compare the network with and without battery.


## References
1. BBC.com, "Tesla mega-battery in Australia activated", 1 December 2017. Url: [link](https://www.bbc.com/news/world-australia-42190358 "https://www.bbc.com/news/world-australia-42190358")
2. Carreras, B. A., et al. "Evidence for self-organized criticality in electric power system blackouts." Proceedings of the 34th Annual Hawaii International Conference on System Sciences. IEEE, 2001.
3. Dobson, Ian, et al. "An initial model for complex dynamics in electric power system blackouts." Proceedings of the 34th annual Hawaii international conference on system sciences. Vol. 2. 2001.
4. Carreras, Benjamin A., et al. "Dynamics, criticality and self-organization in a model for blackouts in power transmission systems." Proceedings of the 35th Annual Hawaii International Conference on System Sciences. IEEE, 2002.
5. Dobson, Ian, et al. "An initial complex systems analysis of the risks of blackouts in power transmission systems." Power Systems and Communications Infrastructures for the future, Beijing, China (2002).
