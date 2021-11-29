# Code for neutral network model of DOC compounds and microbes

This is the computer code corresponding to the following publication:

"Long-term stability of marine dissolved organic carbon emerges from a neutral network of compounds and microbes"  
A. Mentges, C. Feenders, C. Deutsch, B. Blasius, T. Dittmar  
Scientific Reports, 9, Article number: 17780 (2019)  
doi: 10.1038/s41598-019-54290-z  
URL: https://www.nature.com/articles/s41598-019-54290-z

The code is organized as follows:
- the main script is "main_ode_DOM_model.m", it contains three small working examples
- "ode_DOM_model" contains the differential equations
- "wrap_DOM_model" is a wrapper function which calls the ODE solver
- "get_consumption_matrix" and "get_excretion_matrix" generate the network of uptake- and release abilities

Note: we used an adapted version of the standard ode45 solver, which returns non-interpolated time steps for the accurate calculation of DOM age (see "ode45_.m").
