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
class ACTOR_PROTECTION:
    NONE= 1.0
    MASK= 0.1

@dataclass
class InfectionRecord:
# Small record to log past infection(s)
    from_id: int
    to_id: int
    variant_name: str
    time: int

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

        # actor has been  vaccinated.
        self.isVaccinated = False
        # Number of days from vaccination to full level of protection
        self.vaccinationDelay = 0
        # Global time when vaccinated
        self.vaccinationClock = float('-inf')

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

        # The age of the actor (brackets)
        self.ageBracket = 0 

        # currently active infection
        self.myInfection = None
        
        # list of infections
        self.infections = []

    # Infect the individual. Starts as EXPOSED.

    def infect(self, variant, exposer_id):
        self.myInfection = Infection(self, variant)

        self.infectedTime = 0
        self.status = ACTOR_STATUS.EXPOSED
        self.isAsymptomatic = random.random() < variant.asymptomaticRate
        self.willSelfIsolate = random.random() < variant.selfIsolationRate
        
        # Log infection record
        self.infections.append(InfectionRecord(exposer_id, self.id, variant.name, self.simulation.simClock))

    def vaccinate(self, days_ago = None):
        self.isVaccinated = True
        if days_ago is None:
            self.vaccinationClock = self.simulation.simClock
        else:
            self.vaccinationClock = self.simulation.simClock - days_ago
        self.vaccinationDelay = gaussianRandom(self.simulationParameters.vaccinationDelay)
        
    def vaccinationProtection(self, variant):
        ''' vaccinationProtection() returns a multiplier of exposure risk
            A return value of 1.0 means no protection
            A return value of 0.5 means half as likely to get infected
        '''
        
        # TODO: implement waning efficacy over time
        if self.isVaccinated == False:
            # Unvaccinated, so no protection
            return 1.0
        elif (self.isVaccinated == True) and (self.vaccinationClock + self.vaccinationDelay > self.simulation.simClock):
            # Not yet fully vaccinated
            # TODO:  Current behaviro is no protection, but could model a ramp up
            return 1.0
        else:
            # Fully vaccinated, so return the protection to the exposure variant
            full = self.simulationParameters.variantParameters[variant].vaccinationEfficacy
            current = full
            # TODO:  Model waning behavior.  Need a reasonable curve.  Sample code below is exponential decay.
            # current = full * 0.999 ** (self.simulation.simClock - self.vaccinationClock - self.vaccinationDelay)
            return 1.0 - current

        
    def reinfectionProtection(self, variant):
        ''' reinfectionProtection() returns a multiplier of exposure risk
            A return value of 1.0 means no protection
            A return value of 0.5 means half as likely to get infected
        '''
        
        # TODO: implement waning efficacy over time
        if len(self.infections) == 0:
            # Never infected, so no protection
            #print('*')
            return 1.0
        else:
            # Previously infected, so return the max of the reinfection protections of previous variants against the new variant
            efficacies = []
            for i in self.infections:
                efficacies.append(self.simulationParameters.variantParameters[i.variant_name].recoveredResistance[variant])
            # TODO:  Model waning behavior.  Need a reasonable curve.
            # Can use (self.simulation.simClock - i.time), which is the number of days since infection i
            return 1.0 - max(efficacies)

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
                    if self.myInfection.isFatal:
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
