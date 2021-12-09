from dataclasses import dataclass
from enum import Enum,auto
import random
from functools import partial
import scipy.stats

class Actor:
    def __init__(self,simulation) -> None:
        self.simulation=simulation
        self.locations=[]
        self.age=self.simulation.parameters.ageRV()
        self.infected=False

    def infect(self,infect=True):
        self.infected=infect
        if infect:
            self.simulation.infected.append(self)

@dataclass
class ActivityType:
    name:str = ""
    risk:float = 1.0

class Activities(Enum):
    Home=ActivityType("Home",risk=1.0)
    Work=ActivityType("work",risk=0.5)
    School=ActivityType("school",risk=1.5)
    NursingHome=ActivityType("nursinghome",risk=3.5)
    Groceries=ActivityType("groceries",risk=0.5)
    Gym=ActivityType("gym",risk=2.0)


class Location:
    def __init__(self,activity) -> None:
        self.activity=activity
        self.density=1.0
        self.pdf=scipy.stats.norm(0, 10).pdf
        self.actors=[]

    def addActor(self,actor):
        self.actors.append(actor)
        actor.locations.append(self)
    
    # given a list of defaults, create a location for this activity
    def factory(activity,defaults):
        p=defaults[activity]
        location=Location()

@dataclass
class SimulationParameters:
    populationSize:int = 10000
    homeSizeRV=partial(scipy.stats.norm(2,1).rvs)
    workSizeRV=partial(scipy.stats.norm(30,50).rvs)
    ageRV=partial(scipy.stats.norm(1,80).rvs)

class Simulation:
    def __init__(self,parameters) -> None:
        self.parameters=parameters
        self.actors=[]
        self.locations=[]
        self.infected=[]

        # generate actors
        for n in range(self.parameters.populationSize):
            actor=Actor(self)

            self.actors.append(actor)

        # assign everyone to a home
        aiter=iter(self.actors)
        try:
            while True:
                self.locations.append(Location(Activities.Home))
                for n in range(max(1,int(self.parameters.homeSizeRV()))):
                    actor=next(aiter)
                    self.locations[-1].addActor(actor)
        except StopIteration:
            pass

        # assign everyone someplace to go during the day
        adults=list(filter(lambda x: x.age>18,self.actors))
        random.shuffle(adults)
        aiter=iter(adults)
        try:
            while True:
                self.locations.append(Location(Activities.Work))
                for n in range(max(2,int(self.parameters.workSizeRV()))):
                    actor=next(aiter)
                    self.locations[-1].addActor(actor)
                    # todo: periodically add actor to management group
        except StopIteration:
            pass

        for id,actor in enumerate(self.actors):
            actor.id=id
        for id,location in enumerate(self.locations):
            location.id=id

    # test for close contact
    # for each infected person in each location, check for contact biased by
    # gaussian (distance of 2d Uniform distribution)
    def checkContact(self):
        cleanup=[]
        for infected in self.infected:
            if not infected.infected:
                cleanup.append(infected)
                continue
            for loc in infected.locations:
                iidx=loc.actors.index(infected)
                for sidx,susceptible in enumerate(loc.actors):
                    if susceptible.infected:
                        continue
                    p=loc.pdf(abs(iidx-sidx)*loc.density)
                    contact=random.random()<p
                    if(contact):
                        print(f"contact {p:10.5f} {loc.activity}:{loc.id} Actor {infected.id} @{iidx} to {susceptible.id} @{sidx} dist {abs(sidx-iidx)}/{len(loc.actors)}")
        
        # now that we are done iterating on infected, remove actors that are nolonger infected
        for actor in cleanup:
            print(f"remove {actor.id}")
            self.infected.remove(actor)


parameters=SimulationParameters()
simulation=Simulation(parameters)

# for actor in simulation.actors:
#     print(f"Actor {actor.id}")
#     for loc in actor.locations:
#         print(f" {loc.activity.value.name}:{loc.id}")

# for loc in simulation.locations:
#     print(f" {loc.activity.value.name}:{loc.id}={len(loc.actors)}")

#randomly infect some
for actor in random.choices(simulation.actors,k=10):
    actor.infect(True)
    print(f"Infect {actor.id}")

simulation.checkContact()

print(len(simulation.infected))
for infected in simulation.infected:
    infected.infect(False)

simulation.checkContact()

print(len(simulation.infected))
