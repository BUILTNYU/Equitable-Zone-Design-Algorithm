# Equitable-Zone-Design-Algorithm

This repository gives the Equitable Zone Design Algorithm. The Equitable Zone Design Algorithm is an optimization model that can be used to create an aggregated zone district system based on data collected on basic spatial units (BSUs). It is designed to minimize total heterogeneity while also minimizing data margin of error (MOE) of aggregated zones using a controlling MOE threshold. This algorithm can be helpful for public agencies and mobility providers to make decisions and design services based on reliable and equitable data.

The algorithm has two phases: zone district growing and Tabu search. The zone district growing phase involves growing zone districts from seed BSUs, and adding neighboring BSUs that minimize total heterogeneity. The enclave assignment phase assigns any remaining enclave BSUs to the grown aggregated zones, minimizing total heterogeneity. The Tabu search phase involves making swaps on the grown districts to improve total heterogeneity.

This algorithm framework is inspired by previous research and has been tested on a section of Lower Manhattan, as well as all of NYC with ACS data reliability. The algorithm uses a heuristic approach to solve the non-linear integer programming problem, making it an efficient and effective tool for decision making and service design.

Files includes:
- "zoning_functions_v4_2.py": packages of the Equitable Zone Design Algorithm.
- "NYC Equitable Zoning Generation.ipynb": Jupyter Notebook for NYC Equitable Zoning generation using the Equitable Zone Design Algorithm.
- "2010 Census Tracts.zip": NYC census tract shapefile.
- "nyc_nei.csv": Neighboring relationship of NYC census tracts.
