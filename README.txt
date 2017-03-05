This tool is made as an optimisation of existing solution-model of population dynamics.

The purpose of this tool is to generate dadi models and test them against real data. To proceed  with tool you need to give parameters for initial model. Then tool will try to improve it using simulated annealing (optimize by likelihood).

You can insert specific population changes into model. During simulated annealing it can be added actions, modify their parameters or delete some of them.

To run with demo model (population dynamos of cheetah):
python generate_adaptive_model.py