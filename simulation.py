from dataclasses import dataclass
import math
import random
from actor import Actor, ACTOR_STATUS
from util import gaussianRandom


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

    # Mean/STD of # of interactions per day
    numInteractions = 2.5
    numInteractionsSTD = 1.0

    # % of positives that will quarantine effectively (after positive test).
    positiveQuarantineRate = 0.9

    # isolation after positive test or self isolation
    positiveTestIsolationInterval = 21

    # Daily Vacination rate
    vaccinationRate = 0.0035

    # Vaccination delay
    vaccinationDelay = 28

    # Vaccination effectiveness
    vaccinationEfficacy = 0.95

    # Non compliance rate
    # TODO: This parameter is sampled and assigned to actors, but never actively used
    nonCompliantRate = 0.0

    # constructor(props =:) {
    #     Object.assign(this, props)
    #

    ###  Interation Parameters  ##########################################################

    # Mean/STD of transmission per interaction
    # TODO:  Only mean is used currently
    transmissionRate = 0.05
    transmissionRateSTD = 0.1
    
    # TODO: Contact tracing parameters go here

    ###  Infection Parameters  ##########################################################
    
    # The rate of people who are infected but do not show symptoms.
    # TODO:  This is sampled in both Actor and Infection, but never used to control behavior
    asymptomaticRate = 0.2

    # % of symptomatic that will self isolate if symptomatic.
    # TODO: This parameter is sampled and assigned to actors, but never actively used
    selfIsolationRate = 0.0

    # Days to contagious (int)
    daysToContagious = 2.5
    daysToContagiousSTD = 0.5

    # Days to Recovery
    daysToRecovery = 10
    daysToRecoverySTD = 4

    # Days to symptos (Float)
    daysToSymptoms = 5.5
    daysToSymptomsSTD = 2

    # days_to_pcr detectable. Sampled for this actor.
    daysToPcrDetectable = 2
    daysToPcrDetectableSTD = 0.5

    # The duration of virus shedding when antigen detection is positive. Sampled for this actor.
    durationDaysOfPcrDetection = 24
    durationDaysOfPcrDetectionSTD = 7

    # days until the rapid test will detect. Sampled for this actor.
    daysToAntigenDetectable = 3
    daysToAntigenDetectableSTD = 1

    # The duration of virus shedding when antigen detection is positive. Sampled for this actor.
    durationDaysOfAntigenDetection = 10
    durationDaysOfAntigenDetectionSTD = 3

    # Mortality rate (Float, 0-1)
    mortalityRate = 0.01

    # Recovered Resistance (%, as a probability?)
    # NOT CURRENTLY USED
    recoveredResistance = 0.98

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

# Activity is a risk modifier. 1.0 is normal, 0.0 is safe, >1.0 is risky
@dataclass
class ACTIVITY:
    NORMAL= 1.0
    SAFE=   0.1


@dataclass
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
        self.daysElapsed = 0

        rows = math.floor(math.sqrt(self.simulationParameters.populationSize))
        for i in range(self.simulationParameters.populationSize):
            a = Actor(self)
            a.xPosition = i % rows
            a.yPosition = math.floor(i / rows)
            a.id = i
            a.isTesting = random.random() < self.simulationParameters.testingRate
            a.isTestingPcr = random.random() < self.simulationParameters.testingRatePcr
            a.isNonCompliant = random.random() < self.simulationParameters.nonCompliantRate
            self.actors.append(a)

        # Initial exposed subpopulation
        for cnt in range(int(max(1,
                                 self.simulationParameters.startingInfectionRate * self.simulationParameters.populationSize))):
            idx = math.floor(random.random() * len(self.actors))
            self.actors[idx].infect()
            self.totals.infected += 1

        # The remaining susceptible, after we've created the initially infected
        self.totals.susceptible = self.simulationParameters.populationSize - self.totals.infected

    # *
    # Models transmission from an infected individual to a susceptible
    # individual.
    # duration is in days (default is 15 minutes 0.0104)
    #
    # TODO: model full viral load dynamics and exposure duration

    def hasBeenExposed(self, susceptible, infected, duration=0.0104, activity=ACTIVITY.NORMAL):
        if (infected.status != ACTOR_STATUS.INFECTIOUS
                or susceptible.status != ACTOR_STATUS.SUSCEPTIBLE
        ):
            return False

        if (susceptible.isVaccinatedProtected):
            return False

        if (infected.isolated):
            return False

        if (random.random() < self.simulationParameters.transmissionRate
                * infected.protection
                * activity
                * duration / 0.0104
        ):
            susceptible.infect()
            return True

        return False

    # Handles the disease progression in all actors

    def tickDisease(self, days=1):
        self.daysElapsed += days

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
            other.infect()
            return True

        if (self.hasBeenExposed(actor, other)):
            # print((other.id, 'infects', actor.id)
            actor.infect()
            return True

        return False

    #  Generate daily interactions based on simulation parameters.
    #  This is not used if interactions are based on collision detection.

    def tickInteractions(self, days=1.0):
        for actor in self.actors:
            if (actor.status == ACTOR_STATUS.INFECTIOUS and not actor.isolated):
                # Determine if we infect based on # of interactions and % of day passed
                if (random.random() < days):
                    interactions = gaussianRandom(self.simulationParameters.numInteractions,
                                                  self.simulationParameters.numInteractionsSTD)
                    for encounter in range(int(interactions)):
                        # Pick a random nearby actor to expose
                        other = self.actors[math.floor(random.random() * len(self.actors))]
                        self.checkExposure(other, actor)

    #   This is the outer tick. To be overriden by subclasses.
    #   Should implement policies such as social distancing,
    #   isolation or testing.
    #   This base model just picks random actors to infect

    def tick(self, days=1):
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
