import random
from util import gaussianRandom


class Infection:
    # The actor who is infected
    myActor = None

    # The day the infection started. Date of exposure.
    # TODO: Start doing this off a timestamp instead of a tick

    startDay = None

    # Number of days infected
    infectedTime = 0

    # The days until contagious. Sampled for this actor.
    daysToContagious = -1

    # Whether this infection will be asymptomatic. Sampled for this actor.
    asymptomatic = False

    # If symptomatic, the days until symptomatic. Sampled for this actor.
    daysToSymptomatic = -1

    # days_to_pcr detectable. Sampled for this actor.
    daysToPcrDetectable = -1

    # days until the rapid test will detect. Sampled for this actor.
    daysToAntigenDetectable = -1

    # The last day of contagious. Sampled for this actor.
    daysToNotContagious = -1

    # The last day of symptoms. Sampled for this actor.
    daysToNotSymptomatic = -1

    # The last day of detectable virus shedding using PCR test. Sampled for this actor.
    daysToPcrNotDetectable = -1

    # The last day of detectable virus shedding using antigen test. Sampled for this actor.
    daysToAntigenNotDetectable = -1

    def __init__(self, actor):
        self.myActor = actor
        self.infectedTime = 0
        self.asymptomatic = (random.random() < self.myActor.simulationParameters.asymptomaticRate)

        # TODO: Sample these values in the future!
        params = self.myActor.simulationParameters
        self.daysToContagious = gaussianRandom(params.daysToContagious, params.daysToContagiousSTD)
        self.daysToNotContagious = self.daysToContagious + gaussianRandom(
            params.daysToRecovery - params.daysToContagious, params.daysToRecoverySTD)
        self.daysToSymptomatic = gaussianRandom(params.daysToSymptoms, params.daysToSymptomsSTD)
        self.daysToNotSymptomatic = self.daysToSymptomatic + gaussianRandom(params.daysToRecovery,
                                                                            params.daysToRecoverySTD)
        self.daysToPcrDetectable = gaussianRandom(params.daysToPcrDetectable, params.daysToPcrDetectableSTD)
        self.daysToPcrNotDetectable = self.daysToPcrDetectable + gaussianRandom(params.durationDaysOfPcrDetection,
                                                                                params.durationDaysOfPcrDetectionSTD)
        # this assures that the sampled daysToAntigenDetectible is after the sampled daysToPcrDetectible
        self.daysToAntigenDetectable = self.daysToPcrDetectable + gaussianRandom(
            params.daysToAntigenDetectable - params.daysToPcrDetectable, params.daysToAntigenDetectableSTD)
        self.daysToAntigenNotDetectable = self.daysToAntigenDetectable + gaussianRandom(
            params.durationDaysOfAntigenDetection, params.durationDaysOfAntigenDetectionSTD)

    def tick(self, days=1.0):
        self.infectedTime += days

    def duration(self):
        return self.infectedTime

    def detectPcrTest(self):
        return (self.duration() > self.daysToPcrDetectable and
                self.duration() < self.daysToPcrNotDetectable)

    def detectRapidTest(self):
        return (self.duration() > self.daysToAntigenDetectable and
                self.duration() < self.daysToAntigenNotDetectable)

    def isSymptomatic(self):
        if (self.asymptomatic):
            return False

        return (self.duration() > self.daysToSymptomatic and
                self.duration() < self.daysToNotSymptomatic)

    def isContagious(self):
        return (self.duration() > self.daysToContagious and
                self.duration() < self.daysToNotContagious)
