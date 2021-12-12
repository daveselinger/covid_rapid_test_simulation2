import random
from util import gaussianRandom


class Infection:

    def __init__(self, actor, variant):
        # The actor who is infected
        self.myActor = actor

        # The day the infection started. Date of exposure.
        # TODO: Start doing this off a timestamp instead of a tick
        startDay = None

        # Number of days infected
        self.infectedTime = 0
        
        # Assign variant
        self.variant = variant

        # TODO: Sample these values in the future!

        # Whether this infection will be asymptomatic. Sampled for this actor.
        self.asymptomatic = (random.random() < variant.asymptomaticRate)
        # The days until contagious. Sampled for this actor.
        self.daysToContagious = gaussianRandom(variant.daysToContagious, variant.daysToContagiousSTD)
        self.daysToNotContagious = self.daysToContagious + gaussianRandom(
            variant.daysToRecovery - variant.daysToContagious, variant.daysToRecoverySTD)
        # If symptomatic, the days until symptomatic. Sampled for this actor.
        self.daysToSymptomatic = gaussianRandom(variant.daysToSymptoms, variant.daysToSymptomsSTD)
        self.daysToNotSymptomatic = self.daysToSymptomatic + gaussianRandom(variant.daysToRecovery,
                                                                            variant.daysToRecoverySTD)
        # days_to_pcr detectable. Sampled for this actor.
        self.daysToPcrDetectable = gaussianRandom(variant.daysToPcrDetectable, variant.daysToPcrDetectableSTD)
        self.daysToPcrNotDetectable = self.daysToPcrDetectable + gaussianRandom(variant.durationDaysOfPcrDetection,
                                                                                variant.durationDaysOfPcrDetectionSTD)
        # days until the rapid test will detect. Sampled for this actor.
        # this assures that the sampled daysToAntigenDetectible is after the sampled daysToPcrDetectible
        self.daysToAntigenDetectable = self.daysToPcrDetectable + gaussianRandom(
            variant.daysToAntigenDetectable - variant.daysToPcrDetectable, variant.daysToAntigenDetectableSTD)
        self.daysToAntigenNotDetectable = self.daysToAntigenDetectable + gaussianRandom(
            variant.durationDaysOfAntigenDetection, variant.durationDaysOfAntigenDetectionSTD)

        self.isFatal = (random.random() < params.infectionFatalityRateByAge[self.myActor.ageBracket])

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
