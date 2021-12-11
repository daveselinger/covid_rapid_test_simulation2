from simulation import SimulationParameters, Simulation

def run_sim():
    parameters = SimulationParameters()
    parameters.startingInfectionRate=0.01
    parameters.variantParameters['omicron'].transmissionRate = 0.21
    simulation=Simulation(parameters)
    susceptible,infected,recovered,deceased=[],[],[],[]
    for i in range(100):
        simulation.tick()
        susceptible.append(simulation.totals.susceptible)
        infected.append(simulation.totals.infected)
        recovered.append(simulation.totals.recovered)
        deceased.append(simulation.totals.deceased)

# def parse_args_and_run():
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-c", "--cfg", required=True, type=pathlib.Path)
    # args, overrides = parser.parse_known_args()
    # return run_sim(overrides, **vars(args))


if __name__ == "__main__":
    # parse_args_and_run()
    run_sim()
