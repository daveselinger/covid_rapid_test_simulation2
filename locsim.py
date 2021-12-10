from math import sqrt
from dataclasses import dataclass
from enum import Enum,auto
import random
from functools import partial
import scipy.stats

class Actor:
    def __init__(self,simulation) -> None:
        self.simulation=simulation
        self.locationtimes=[]
        self.age=self.simulation.parameters.ageRV()
        self.infected=False
        self.infectionHistory=[]
        self.transmissionHistory=[]

    def infect(self,infected,location=None,when=None):
        if infected is None:
            self.infected=False
        else:
            self.infected=True
            if when is None:
                when = self.simulation.timestamp
            assert location is not None,"infect requires a location"
            self.infectionHistory.append((infected,location,when))
            infected.transmissionHistory.append((self,location,when))
            self.simulation.infected.append(self)



@dataclass
class ActivityType:
    name:str = ""
    risk:float = 1.0

class Activities(Enum):
    External=ActivityType("external",1.0)
    Home=ActivityType("home",1.0)
    Work=ActivityType("work",0.5)
    Admin=ActivityType("admin",0.5)
    School=ActivityType("school",1.5)    
    SchoolAdmin=ActivityType("schooladmin",1.5)
    NursingHome=ActivityType("nursinghome",3.5)
    Groceries=ActivityType("groceries",0.5)
    Gym=ActivityType("gym",2.0)
    Management=ActivityType("management",1.0)

class Location:
    def __init__(self,activity,density=1.0,sigma=10) -> None:
        self.id=-1
        self.activity=activity
        self.density=density
        # todo: What is a reasonable number for sigma and should it be a parameter
        self.pdf=scipy.stats.norm(0, sigma).pdf
        self.actors=[]

    def addActor(self,actor,interval):
        self.actors.append(actor)
        actor.locationtimes.append((self,interval))

@dataclass
class SimulationParameters:
    populationSize:int = 10000
    dayInterval=0.4
    nightInterval=0.5
    homeSizeRV=partial(scipy.stats.norm(2,1).rvs)
    workGroupSizeRV=partial(scipy.stats.norm(15,5).rvs)
    workGroupCntRV=partial(scipy.stats.norm(3,10).rvs)
    schoolClassSizeRV=partial(scipy.stats.norm(25,10).rvs)
    schoolClassCntRV=partial(scipy.stats.norm(25,10).rvs)
    ageRV=partial(scipy.stats.uniform(1,80).rvs)


class Simulation:
    def __init__(self,parameters) -> None:
        self.parameters=parameters
        self.actors=[]
        self.locations=[]
        self.infected=[]
        self.externalActor=Actor(self)
        self.externalLocation=Location(Activities.External)
        self.timestamp=0.0

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
                    self.locations[-1].addActor(actor,self.parameters.nightInterval)
        except StopIteration:
            pass

        # assign jobs to adults
        adults=list(filter(lambda x: x.age>18,self.actors))
        random.shuffle(adults)
        aiter=iter(adults)

        # assign kids to school during the day
        kids=list(filter(lambda x: x.age<=18,self.actors))
        random.shuffle(kids)
        kiter=iter(kids)
        schoolClassCnt=0
        try:
            while True:                 
                if schoolClassCnt <= 0:
                    #create a new school 
                    admin=Location(Activities.SchoolAdmin)
                    self.locations.append(admin)
                    schoolClassCnt=self.parameters.schoolClassCntRV()
                schoolClassCnt-=1
                
                #create a classroom
                self.locations.append(Location(Activities.School,sigma=100))
                # assign a teacher
                teacher=next(aiter)
                self.locations[-1].addActor(teacher,self.parameters.dayInterval)
                admin.addActor(teacher,self.parameters.dayInterval)  
                for n in range(max(2,int(self.parameters.schoolClassSizeRV()))):
                    actor=next(kiter)
                    self.locations[-1].addActor(actor,self.parameters.dayInterval)
        except StopIteration:
            pass

        # assign everyone someplace to go during the day
        # work management team
        workgroupCnt=0
        try:
            while True:
                #create a new work space
                if workgroupCnt <= 0:
                    admin=Location(Activities.Admin)
                    self.locations.append(admin)
                    workgroupCnt=self.parameters.workGroupCntRV()
                workgroupCnt-=1
                self.locations.append(Location(Activities.Work,sigma=10))

                # assign a manager to group and admin team
                manager=next(aiter)
                self.locations[-1].addActor(manager,self.parameters.dayInterval)
                admin.addActor(manager,self.parameters.dayInterval)  
                for n in range(max(2,int(self.parameters.workGroupSizeRV()))):
                    actor=next(aiter)
                    self.locations[-1].addActor(actor,self.parameters.dayInterval)

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
        cleanup=0
        for infected in self.infected:
            if not infected.infected:
                cleanup+=1
                continue
            for loc,iinterval in infected.locationtimes:
                iidx=loc.actors.index(infected)
                for sidx,susceptible in enumerate(loc.actors):
                    # technically you might infect with a new variant, but simplify for now
                    if susceptible.infected:  
                        continue
                    p=loc.pdf(abs(iidx-sidx)*loc.density)
                    contact=random.random()<p
                    if(contact):
                        for tmp,sinterval in susceptible.locationtimes:
                            if tmp == loc:
                                break
                        else:
                            assert "location not in susceptibles locationtime list"

                        #estimate time overlap
                        interval=sqrt(iinterval*sinterval)

                        print(f"contact {p:10.5f} {loc.activity}:{loc.id} Actor {infected.id} @{iidx} to {susceptible.id} @{sidx} dist {abs(sidx-iidx)}/{len(loc.actors)} times {interval} {sinterval}")
                        susceptible.infect(infected,loc)

        # now that we are done iterating on infected, remove actors that are nolonger infected
        if cleanup > 100:
            self.infected=[actor for actor in self.infected if actor.infected]

def selftests(simulation):
    for actor in simulation.actors:
        print(f"Actor {actor.id}")
        for loc,t in actor.locationtimes:
            print(f" {loc.activity.value.name}:{loc.id} {t:0.2f}")

    for loc in simulation.locations:
        print(f" {loc.activity.value.name}:{loc.id}={len(loc.actors)}")

    #randomly infect some actors

    for actor in random.choices(simulation.actors,k=10):
        actor.infect(simulation.externalActor,simulation.externalLocation)
        print(f"Infect {actor.id}")

    for i in range(4):
        simulation.timestamp+=1
        simulation.checkContact()

    def transmissionDump(actor,depth=0):
        if len(actor.transmissionHistory) == 0:
            return
        indent="   "*depth
        for other,oloc,otime in actor.transmissionHistory:
            print(f"{oloc.activity.name:10}, {oloc.id:5d}, {otime:8.4f}, {indent} {other.id:6d}")
            transmissionDump(other,depth+1)

    transmissionDump(simulation.externalActor)

parameters=SimulationParameters()
simulation=Simulation(parameters)
selftests(simulation)