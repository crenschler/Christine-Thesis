#time parameters
start_time = 0
end_time = 100

#how much time passes between each successive calculation
time_step = 1/4 #in years
end_step = int(((end_time-start_time)/time_step))

#number of rabbits when simulation is started
initial_rabbits = 30000

#number of foxes when the simulation is started
initial_foxes = 15

#number of rabbits that must be killed for a fox to be born
rabbits_killed_per_fox_birth = 1000

#chance a rabbit will die when rabbit fox cross paths
chance_a_rabbit_will_die_during_a_meeting = 0.50

#the chance a rabbit and fox will cross paths
chance_of_rabbit_and_fox_meeting = 0.02

#percent of the rabbit poplation that will be born this time step
rabbit_growth_rate = 0.20

#the percent of the fox pop that will die this time step from old age
fox_death_rate = 0.10

#initialization
rabbits_over_time = fill(0.0, end_step+1)
foxes_over_time = fill(0.0, end_step+1)
model_time = fill(0.0, end_step+1)

rabbits = initial_rabbits
foxes = initial_foxes

rabbits_over_time[1] = rabbits
foxes_over_time[1] = foxes

for sim_step = 1:end_step
#get the time from the step
sim_time = start_time + sim_step * time_step
model_time[sim_step] = sim_time

#calculate our rates
rabbit_births = rabbits * rabbit_growth_rate
rabbits_eaten = min(rabbits, chance_a_rabbit_will_die_during_a_meeting * chance_of_rabbit_and_fox_meeting * foxes * rabbits)

fox_births = 1/rabbits_killed_per_fox_birth * rabbits_eaten
fox_deaths = foxes * fox_death_rate

foxes = foxes + fox_births - fox_deaths
rabbits = rabbits + rabbit_births - rabbits_eaten

rabbits_over_time[sim_step+1] = rabbits
foxes_over_time[sim_step+1] = foxes
end

println("Time,Rabbits (Thousands),Foxes")
for i = 1:end_step
        print(model_time[i])
        print(",")
        print(rabbits_over_time[i]/1000)
        print(",")
        println(foxes_over_time[i])
        end
