import random
from util import gaussianRandom


class Infection:

    def __init__(self, actor, variant = None):
        # The actor who is infected
        self.myActor = actor

        # The day the infection started. Date of exposure.
        # TODO: Start doing this off a timestamp instead of a tick
        startDay = None

        # Number of days infected
        self.infectedTime = 0

        # TODO: Sample these values in the future!
        params = self.myActor.simulationParameters
        
        # Assign variant
        if variant is not None:
            self.variant = variant
        else:
            rnd = random.random()
            for v, pct in params.startingVariantMix.items():
                rnd -= pct
                if rnd <= 0:
                    self.variant = v
                    break

        # Whether this infection will be asymptomatic. Sampled for this actor.
        self.asymptomatic = (random.random() < self.myActor.simulationParameters.asymptomaticRate)
        # The days until contagious. Sampled for this actor.
        self.daysToContagious = gaussianRandom(params.daysToContagious, params.daysToContagiousSTD)
        self.daysToNotContagious = self.daysToContagious + gaussianRandom(
            params.daysToRecovery - params.daysToContagious, params.daysToRecoverySTD)
        # If symptomatic, the days until symptomatic. Sampled for this actor.
        self.daysToSymptomatic = gaussianRandom(params.daysToSymptoms, params.daysToSymptomsSTD)
        self.daysToNotSymptomatic = self.daysToSymptomatic + gaussianRandom(params.daysToRecovery,
                                                                            params.daysToRecoverySTD)
        # days_to_pcr detectable. Sampled for this actor.
        self.daysToPcrDetectable = gaussianRandom(params.daysToPcrDetectable, params.daysToPcrDetectableSTD)
        self.daysToPcrNotDetectable = self.daysToPcrDetectable + gaussianRandom(params.durationDaysOfPcrDetection,
                                                                                params.durationDaysOfPcrDetectionSTD)
        # days until the rapid test will detect. Sampled for this actor.
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
