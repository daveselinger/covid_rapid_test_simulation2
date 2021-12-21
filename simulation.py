from dataclasses import dataclass
import math
import random
from actor import Actor, ACTOR_STATUS, InfectionRecord
from util import gaussianRandom
import pandas as pd


#Demographics
ageBrackets = [4,17,29,39,49,64,74,84,110]
usPopulationByAge=[0.047,0.163,0.162,0.136,0.123,0.129,0.101,0.053,0.023]

# approx from Dec 2020. Crude interpolation
#  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7721859/pdf/10654_2020_Article_698.pdf
infectionFatalityRateByAge=[0.00004,0.00004,0.00004,0.00068,0.023,0.0775,0.025,0.085,0.283]

@dataclass
class SimulationParameters :
    # Population is made up of People.
    # People are called Actors.  
    # Actors have Interactions with other Actors.  
    # Interactions may lead to Exposures/Infections.
    # Testing is either Rapid antigen or PCR

    ###  Population Parameters  ##########################################################

    # Starting Population size (Int)
    populationSize = 10000

    # Starting infected rate (Float, 0-1)
    startingInfectionRate = 0.004

    # Starting recovered rate
    # A list with one tuple per variant.  Each tuple has variant, rate (float 0-1), mean days ago, std dev days ago.
    startingRecoveredList = [('beta',0.0010,510, 60), ('delta',0.0020,180,60)]

    # % of positives that will quarantine effectively (after positive test).
    positiveQuarantineRate = 0.9

    # isolation after positive test or self isolation
    positiveTestIsolationInterval = 21

    # Starting vaccinated rate (Float, 0-1)
    startingVaccinationRate = 0.600
    # Mean/STD of how long ago these were vaccinated (in days)
    vaccinationMean = 240
    vaccinationSTD = 60
    # Daily Vacination rate
    vaccinationRate = 0.0035

    # Vaccination delay
    vaccinationDelay = 28

    # Non compliance rate
    # TODO: This parameter is sampled and assigned to actors, but never actively used
    nonCompliantRate = 0.0

    # constructor(props =:) {
    #     Object.assign(this, props)
    #

    ###  Interation Parameters  ##########################################################
    
    # Mean/STD of # of interactions per day
    numInteractions = 2.5
    numInteractionsSTD = 1.0

    # Rate of actors interacting with people external to the simulation (Float, 0-1)
    externalInteractionRate = 0.05
    
    # Rate at which external people are infected
    externalBaseInfected  = 0.01
    
    # TODO: Contact tracing parameters go here

    ###  Infection Parameters  ##########################################################
    
    # Variants (dict)
    # This dict should contain all variants in the simulation.  
    # Use a zero probability for any variants not present at the start
    startingVariantMix = { 'delta': 0.75, 'omicron': 0.05, 'beta':0.17, 'alpha':0.03}

    # The variantParameters dict will get populated later
    variantParameters = {}

    ###  Rapid Testing Parameters  ##########################################################

    # Testing interval: Start out with every 3 days (Float: 0-100.0)
    testingInterval = 3.0

    # Testing rate
    testingRate = 0.6

    # Random rapid testing rate
    testingRateRandom = 0.0

    # False positive % (Float, 0-1)
    falsePositiveRate = 0.01

    # False negative % (Float, 01)
    falseNegative = 0.02

    ###  PCR Testing Parameters  ##########################################################

    # Testing rate PCR
    testingRatePcr = 0.000

    # RandomTesting rate PCR
    testingRateRandomPcr = 0.003

    # Testing interval for PCR: Start out with every 3 days (Float: 0-100.0)
    testingIntervalPcr = 3.0

    # False positive % for PCR (Float, 0-1)
    falsePositiveRatePcr = 0.01

    # False negative % for PCR (Float, 01)
    falseNegativePcr = 0.02

    # The delay from PCR test to results. Isolation is delayed by this ammount
    daysToPcrResults = 1.5


    # Days to detectable (Float)
    # NOT CURRENTLY USED
    daysToDetectable = 2.5
    daysToDetectableSTD = 0.5

    #Demographics
    ageBrackets = [4,17,29,39,49,64,74,84,110]

    # US population ratios
    populationByAge=[0.047,0.163,0.162,0.136,0.123,0.129,0.101,0.053,0.023]

    def __init__(self):
        # Create a dictionary of variant parameters
        for v in self.startingVariantMix:
            self.variantParameters[v] = VariantParameters(v)
            
        # Customize the parameters for the different variants
        self.variantParameters['alpha'].transmissionRate /= 2
        self.variantParameters['omicron'].transmissionRate *= 2
        self.variantParameters['delta'].recoveredResistance['delta'] = 0.95
        self.variantParameters['omicron'].recoveredResistance['omicron'] = 0.95



class VariantParameters :

    def __init__(self, name = 'Default'):
        self.name = name

        ###  Interation Parameters  ##########################################################

        # Mean/STD of transmission per interaction
        # TODO:  Only mean is used currently
        self.transmissionRate = 0.2
        self.transmissionRateSTD = 0.1
        
        ###  Infection Parameters  ##########################################################
        
        # The rate of people who are infected but do not show symptoms.
        # TODO:  This is sampled in both Actor and Infection, but never used to control behavior
        self.asymptomaticRate = 0.2

        # % of symptomatic that will self isolate if symptomatic.
        # TODO: This parameter is sampled and assigned to actors, but never actively used
        self.selfIsolationRate = 0.0

        # Days to contagious (int)
        self.daysToContagious = 2.5
        self.daysToContagiousSTD = 0.5

        # Days to Recovery
        self.daysToRecovery = 10
        self.daysToRecoverySTD = 4

        # Days to symptos (Float)
        self.daysToSymptoms = 5.5
        self.daysToSymptomsSTD = 2

        # days_to_pcr detectable. Sampled for this actor.
        self.daysToPcrDetectable = 2
        self.daysToPcrDetectableSTD = 0.5

        # The duration of virus shedding when antigen detection is positive. Sampled for this actor.
        self.durationDaysOfPcrDetection = 24
        self.durationDaysOfPcrDetectionSTD = 7

        # days until the rapid test will detect. Sampled for this actor.
        self.daysToAntigenDetectable = 3
        self.daysToAntigenDetectableSTD = 1

        # The duration of virus shedding when antigen detection is positive. Sampled for this actor.
        self.durationDaysOfAntigenDetection = 10
        self.durationDaysOfAntigenDetectionSTD = 3

        # Recovered Resistance (%, as a probability?)
        self.recoveredResistance = {'alpha': 0.95, 'beta': 0.95, 'delta': 0.9, 'omicron': 0.8}

        # approx from Dec 2020. Crude interpolation
        #  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7721859/pdf/10654_2020_Article_698.pdf
        self.infectionFatalityRateByAge=[0.00004,0.00004,0.00004,0.00068,0.023,0.0775,0.025,0.085,0.283]

        ###  Vaccination Parameters  ##########################################################

        # Vaccination effectiveness
        self.vaccinationEfficacy = 0.95
        
    

# Activity is a risk modifier. 1.0 is normal, 0.0 is safe, >1.0 is risky
class ACTIVITY:
    NORMAL= 1.0
    SAFE=   0.1


class RunStatistics:
    susceptible = 0
    infected = 0
    recovered = 0
    deceased = 0
    testsConducted = 0
    daysLost = 0


class Simulation:
    def __init__(self, simulationParameters):
        self.simulationParameters = simulationParameters
        self.actors = []
        self.totals = RunStatistics()
        self.simClock = 0

        rows = math.floor(math.sqrt(self.simulationParameters.populationSize))
        for i in range(self.simulationParameters.populationSize):
            a = Actor(self)
            a.xPosition = i % rows
            a.yPosition = math.floor(i / rows)
            a.id = i
            a.isTesting = random.random() < self.simulationParameters.testingRate
            a.isTestingPcr = random.random() < self.simulationParameters.testingRatePcr
            a.isNonCompliant = random.random() < self.simulationParameters.nonCompliantRate
            a.ageBracket =  random.choices(range(len(self.simulationParameters.ageBrackets)),
                self.simulationParameters.ageBrackets)[0]
            self.actors.append(a)

        # Initial infected subpopulation
        exposed_list = random.sample(range(len(self.actors)), 
                                     int(max(1, self.simulationParameters.startingInfectionRate * self.simulationParameters.populationSize)))
        for idx in exposed_list:
           # Choose variant randomly according to starting mix
            variant = random.choices(list(self.simulationParameters.variantParameters.values()), list(self.simulationParameters.startingVariantMix.values()))[0]
            self.actors[idx].infect(variant, -1)    # Initial exposures get dummy ID of -1
            self.totals.infected += 1

        # Initial recovered subpopulation
        for variant, startingRecoveredRate, recoveredDaysMean, recoveredDaysSTD in self.simulationParameters.startingRecoveredList:
            recovered_list = random.sample(range(len(self.actors)), 
                                           int(max(1, startingRecoveredRate * self.simulationParameters.populationSize)))
            for idx in recovered_list:
                recoveredDays = 0 - gaussianRandom(recoveredDaysMean, recoveredDaysSTD)
                if recoveredDays > -2:
                    recoveredDays = -2
                self.actors[idx].infections.append(InfectionRecord(-1, idx, variant, recoveredDays))
                self.totals.recovered += 1

        # Initial vaccinated subpopulation
        # TODO:  Currently these are random and independent of the starting infected population.  Maybe they should be negatively correlated.
        vaccinated_list = random.sample(range(len(self.actors)), 
                                        int(max(1, self.simulationParameters.startingVaccinationRate * self.simulationParameters.populationSize)))
        for idx in vaccinated_list:
            vaccinatedDays = gaussianRandom(self.simulationParameters.vaccinationMean, 
                                            self.simulationParameters.vaccinationSTD)
            if vaccinatedDays < 2:
                vaccinatedDays = 2
            self.actors[idx].vaccinate(vaccinatedDays)

        # The remaining susceptible, after we've created the initially infected
        self.totals.susceptible = self.simulationParameters.populationSize - self.totals.infected - self.totals.recovered
        
    # *
    # Models transmission from an infected individual to a susceptible
    # individual.
    # duration is in days (default is 15 minutes 0.0104)
    #
    # TODO: model full viral load dynamics and exposure duration

    def hasBeenExposed(self, susceptible, infected=None, duration=0.0104, activity=ACTIVITY.NORMAL, variant=None):
        if (((infected is not None) and (infected.status != ACTOR_STATUS.INFECTIOUS))
                or (susceptible.status != ACTOR_STATUS.SUSCEPTIBLE and susceptible.status != ACTOR_STATUS.RECOVERED)
        ):
            return False

        if ((infected is not None) and (infected.isolated)):
            return False
            
        if infected is None:
            infected_protection = 1.0
        else:
            infected_protection = infected.protection
            
        if variant is None:
            variant = infected.myInfection.variant

        if (random.random() < variant.transmissionRate
                * infected_protection
                * susceptible.vaccinationProtection(variant.name)
                * susceptible.reinfectionProtection(variant.name)
                * activity
                * duration / 0.0104
        ):
            return True

        return False

    # Handles the disease progression in all actors

    def tickDisease(self, days=1):
        newTotals = RunStatistics()

        # This handles disease progression
        # TODO: cost model
        for actor in self.actors:
            actor.tick(days)
            # Update totals
            if (actor.status == ACTOR_STATUS.RECOVERED):
                newTotals.recovered += 1
            elif (actor.status == ACTOR_STATUS.SUSCEPTIBLE):
                newTotals.susceptible += 1
            elif (actor.status == ACTOR_STATUS.INFECTIOUS):
                newTotals.infected += 1
            elif (actor.status == ACTOR_STATUS.EXPOSED):
                newTotals.infected += 1
            elif (actor.status == ACTOR_STATUS.DECEASED):
                newTotals.deceased += 1

            newTotals.testsConducted += actor.testsConducted
            newTotals.daysLost += actor.daysIsolated

        self.totals = newTotals

    # Check for exposure in either direction and infect the susceptible actor
    # if exposure occured. This can be called for collision detect based interactions.
    # @returns:Boolean - Whether there was a resulting infection `TRUE` or not `FALSE`.

    def checkExposure(self, actor, other, duration=0.0104, activity=ACTIVITY.NORMAL):
        if (self.hasBeenExposed(other, actor)):
            # print((actor.id, 'infects', other.id)
            other.infect(actor.myInfection.variant, actor.id)
            return True

        if (self.hasBeenExposed(actor, other)):
            # print((other.id, 'infects', actor.id)
            actor.infect(other.myInfection.variant, other.id)
            return True

        return False

    def externalInteraction(self):
        ''' Helper function to determine external interactions 
            If desired, both rate and variants could vary over time
        '''
        external_infected_rate = self.simulationParameters.externalBaseInfected    # This rate could vary with time using self.simClock
        if random.random() < self.simulationParameters.externalInteractionRate * external_infected_rate:
            # The variant weighting could vary with time
            variant = random.choices(list(self.simulationParameters.variantParameters.values()), list(self.simulationParameters.startingVariantMix.values()))[0]
            return True, variant
        return False, None
    
    #  Generate daily interactions based on simulation parameters.
    #  This is not used if interactions are based on collision detection.

    def tickInteractions(self, days=1.0):
        for actor in self.actors:
            if (actor.status == ACTOR_STATUS.INFECTIOUS and not actor.isolated):
                # Determine if we infect based on # of interactions and % of day passed
                if (random.random() < days):
                    interactions = int(gaussianRandom(self.simulationParameters.numInteractions,
                                                  self.simulationParameters.numInteractionsSTD))
                    if interactions < 0:
                        interactions = 0
                    encounter_list = random.sample(range(len(self.actors)), int(interactions))
                    for idx in encounter_list:
                        self.checkExposure(self.actors[idx], actor)
            elif ((actor.status == ACTOR_STATUS.SUSCEPTIBLE) or (actor.status == ACTOR_STATUS.SUSCEPTIBLE)) and (not actor.isolated):
                # Determine if actor with have an encounter with an external source of infection
                interacts, variant = self.externalInteraction()
                # If we've had an external encounter, use hasBeenExposed() to check for transmission
                if interacts and self.hasBeenExposed(actor, variant=variant):
                    actor.infect(variant, -1)    # External exposures get dummy ID of -1

    #   This is the outer tick. To be overriden by subclasses.
    #   Should implement policies such as social distancing,
    #   isolation or testing.
    #   This base model just picks random actors to infect

    def tick(self, days=1):
        self.simClock += days

        self.tickInteractions(days)
        self.tickRapidTesting(days)
        self.tickPcrTesting(days)
        self.tickVaccination(days)
        self.tickDisease(days)

    # Implement daily rapid testing policy

    def tickRapidTesting(self, days=1.0):
        # Perform rapid testing
        for actor in self.actors:
            if (((actor.isTesting and (actor.testTime is None
                                       or actor.testTime >= self.simulationParameters.testingInterval))
                 or (random.random() < self.simulationParameters.testingRateRandom / days))
                    and not actor.isolated
                    and actor.rapidTest()):
                # TODO: Can sample these as well.
                if (random.random() < self.simulationParameters.positiveQuarantineRate):
                    actor.isolateFor(self.simulationParameters.positiveTestIsolationInterval)
                    # print(('Isolated ', actor.id)
                else:
                    pass
                    # print(('Isolation non compliance', actor.id)

            # TODO: Some actors become sick and never become "unsick" so they isolate forever.
            if (actor.isSymptomatic and actor.willSelfIsolate and not actor.isolated):
                actor.isolateFor(self.simulationParameters.positiveTestIsolationInterval)

    # Implement pcr testing policy

    def tickPcrTesting(self, days=1.0):
        # Perform pcr testing
        for actor in self.actors:
            if (
                    ((random.random() < self.simulationParameters.testingRateRandomPcr / days) or
                     (actor.isTestingPcr and (actor.testTimePcr is None or
                                              actor.testTimePcr >= self.simulationParameters.testingIntervalPcr)))
                    and not actor.isolated and actor.pcrTest()):
                # TODO: Can sample these as well.
                if (random.random() < self.simulationParameters.positiveQuarantineRate):
                    # print('Isolated PCR', actor.id)
                    actor.isolateFor(self.simulationParameters.positiveTestIsolationInterval,
                                     self.simulationParameters.daysToPcrResults)
                else:
                    pass
                    # print(('Isolation non compliance', actor.id)

            # TODO: Some actors become sick and never become "unsick" so they isolate forever.
            if (actor.isSymptomatic and actor.willSelfIsolate and not actor.isolated):
                actor.isolateFor(self.simulationParameters.positiveTestIsolationInterval)

    def tickVaccination(self, days=1.0):
        # Perform random vaccination
        # TODO: find a better way to do
        for actor in self.actors:
            if (not actor.isVaccinated and
                    random.random() < self.simulationParameters.vaccinationRate * days):
                actor.vaccinate()
                
    def infectionsDF(self):
        '''Return a pandas dataframe with the infection spread data'''
        vars = []
        for a in self.actors:
            for i in a.infections:
                vars.append((i.from_id, i.to_id, i.variant_name, i.time))

        df = pd.DataFrame(vars, columns=['from_id', 'to_id', 'variant_name', 'time'])
        df.sort_values(['time','from_id'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        return df
