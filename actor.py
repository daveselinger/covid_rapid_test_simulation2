from dataclasses import dataclass
import random
from infection import Infection
from util import gaussianRandom
from enum import Enum



class ACTOR_STATUS (Enum):
    SUSCEPTIBLE= 0
    EXPOSED=     1
    INFECTIOUS=  2
    RECOVERED=   3
    DECEASED=    4

# Protection is 1.0 for no protection, 0 for full protection
@dataclass
class ACTOR_PROTECTION:
    NONE= 1.0
    MASK= 0.1

class Actor:

    def __init__(self, simulation):
        self.simulation = simulation

        # Initialize the day of the next test using uniform random.
        #  self.testTimePcr = math.floor(random.random() * self.simulationParameters.testingIntervalPcr)
        #  self.testTime = math.floor(random.random() * self.simulationParameters.testingInterval)

        # The parameters of the overall simulation for this actor
        self.simulationParameters = simulation.simulationParameters

        # Blue-healthy/Susceptible.
        self.status = ACTOR_STATUS.SUSCEPTIBLE

        # Isolated/Not isolated (impacted by rapid testing).
        self.isolated = False

        # Number of days of isolation remaining.
        self.isolatedRemain = 0

        # Number of days to delay before beginning isolation
        self.isolateAfterRemain = 0

        # The number of days this actor was infected.
        self.infectedTime = None

        # The horizontal position of this actor.
        self.xPosition = 0

        # The vertical position of this actor.
        self.yPosition = 0

        # ID of this actor.
        self.id = 0

        # By default actor has no protection
        self.protection = ACTOR_PROTECTION.NONE

        # Days since most recent rapid test
        self.testTime = None

        # Days since most recent PCR test
        self.testTimePcr = None

        # Days isolated
        self.daysIsolated = 0

        # Whether or not this actor is showing symptoms once infected.
        self.isAsymptomatic = False

        # Whether or not this actor is showing symptoms or not.
        self.isSymptomatic = False

        # Number of days remaining from vaccination to protection.
        self.vaccinationRemain = 0

        # actor has been  vaccinated.
        self.isVaccinated = False

        # actor has vaccine protection.
        self.isVaccinatedProtected = False

        # Whether or not this actor will self-isolate when symptoms appear.
        self.willSelfIsolate = True

        # Nmber of tests conducted.
        self.testsConducted = 0

        # Nmber of Pcr tests conducted.
        self.testsConductedPcr = 0

        # Whether or not this actor will comply with orders.
        self.isNonCompliant = False

        # Is actor part of the rapid test testing program
        self.isTesting = False

        # Is actor part of the Pcr routine testing program
        self.isTestingPcr = False

        self.myInfection = None

    # Infect the individual. Starts as EXPOSED.

    def infect(self, variant):
        self.myInfection = Infection(self, variant)

        self.infectedTime = 0
        self.status = ACTOR_STATUS.EXPOSED
        self.isAsymptomatic = random.random() < variant.asymptomaticRate
        self.willSelfIsolate = random.random() < variant.selfIsolationRate

    def vaccinate(self):
        self.isVaccinated = True
        self.vaccinationRemain = gaussianRandom(self.simulationParameters.vaccinationDelay)

    # Perform rapid test on actor

    def rapidTest(self):
        self.testsConducted += 1
        self.testTime = 0

        if (self.status == ACTOR_STATUS.EXPOSED
                or self.status == ACTOR_STATUS.INFECTIOUS):
            if (self.myInfection.detectRapidTest()
                    and random.random() > self.simulationParameters.falseNegative):
                # print(('Tested positive', self.id)
                return True

        if (random.random() < self.simulationParameters.falsePositiveRate):
            # print(('False positive ', self.id)
            return True

        return False

    # Perform PCR test on actor

    def pcrTest(self):
        self.testsConductedPcr += 1
        self.testTimePcr = 0

        if (self.status == ACTOR_STATUS.EXPOSED or self.status == ACTOR_STATUS.INFECTIOUS):
            if (self.myInfection.detectPcrTest()
                    and random.random() > self.simulationParameters.falseNegativePcr):
                # print('Tested PCR positive', self.id)
                return True

        if (random.random() < self.simulationParameters.falsePositiveRatePcr):
            # print(('False positive ', self.id)
            return True

        return False

    #  Isolate actor for a number of days.
    #  @param:number days - The number of days to isolate (int).

    def isolateFor(self, days, after=0):
        if (after == 0):
            self.isolated = True
        else:
            self.isolateAfterRemain = after

        self.isolatedRemain = days

    #
    # Perform updates to actor for each cycle.

    def tick(self, days=1.0):
        # First progress the status based on lifecycle
        if (self.myInfection is not None and self.status != ACTOR_STATUS.RECOVERED):
            if (self.status == ACTOR_STATUS.EXPOSED):
                if (self.myInfection.isContagious()):
                    self.status = ACTOR_STATUS.INFECTIOUS

            elif (self.status == ACTOR_STATUS.INFECTIOUS):
                if (not self.myInfection.isContagious()):
                    if (random.random() < self.myInfection.variant.mortalityRate):
                        self.status = ACTOR_STATUS.DECEASED
                    else:
                        self.status = ACTOR_STATUS.RECOVERED

            self.isSymptomatic = self.myInfection.isSymptomatic()

        if (self.isolateAfterRemain > 0):
            self.isolateAfterRemain -= days
            if (self.isolateAfterRemain <= 0):
                self.isolated = True

        # Advance the clock
        if self.infectedTime is not None:
            self.infectedTime += days

        if (self.myInfection is not None):
            self.myInfection.tick(days)

        if (self.isolated):
            self.daysIsolated += days
            self.isolatedRemain -= days
            if (self.isolatedRemain <= 0):
                self.isolated = False

        if (self.testTime is not None):
            self.testTime += days

        if (self.testTimePcr is not None):
            self.testTimePcr += days

        if (self.vaccinationRemain > 0):
            self.vaccinationRemain -= days

            # This can only happen once.
            if (self.vaccinationRemain <= 0):
                if (random.random() < self.simulationParameters.vaccinationEfficacy):
                    self.isVaccinatedProtected = True
        random.random() < self.simulationParameters.vaccinationEfficacy
