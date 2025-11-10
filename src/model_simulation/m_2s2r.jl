# In this model, we define a 2s2rladder model where each rung is a two-state model; it is the 2rst 
# version where the kon is 0 or konpositive mimicing a very sharp sigmoidal 
# MOREOVER there is no delay in this version but a "degradation" rate
# Gregory Roth
# First version:            15.08.2023
# Modified:                 15.01.2024
#
#-----------------------------------------------------------------------------------------------
module m_2s2r


using JumpProcesses
using Statistics
using Survival
using StatsBase
using DataFrames


export simulate_bursttraj


function simulateModel(u0,tf, parameters)
       # define the rates and transitions   
   # regime transitions
   rate1(u,p,t) = (u[1]==0)*p.kforward
   function affect1!(integrator)
       integrator.u[1] += 1               
       nothing
   end
   cj1 = ConstantRateJump(rate1, affect1!)

   rate2(u,p,t) = (u[1]==1)*p.kbackward
   function affect2!(integrator)
       integrator.u[1] -= 1               
       nothing
   end
   cj2 = ConstantRateJump(rate2, affect2!)

   # ON transition
   rate3(u,p,t) = (u[2]==0)*(u[1]==0)*p.kon0
   function affect3!(integrator)
       integrator.u[2] += 1               
       nothing
   end
   cj3 = ConstantRateJump(rate3, affect3!)

   rate4(u,p,t) = (u[2]==0)*(u[1]==1)*p.kon1
   function affect4!(integrator)
       integrator.u[2] += 1               
       nothing
   end
   cj4 = ConstantRateJump(rate4, affect4!)


   # OFF transition
   rate5(u,p,t) = (u[2]==1)*p.koff
   function affect5!(integrator)
       integrator.u[2] -= 1               
       nothing
   end
   cj5 = ConstantRateJump(rate5, affect5!)

   # initiation
   rate6(u,p,t) = (u[2]==1)*p.kinitiation
   function affect6!(integrator)
       integrator.u[3] += 1               
       nothing
   end
   cj6 = ConstantRateJump(rate6, affect6!)

   # degradation
   rate7(u,p,t) = u[3]*p.kdegradation
   function affect7!(integrator)
       integrator.u[3] -= 1               
       nothing
   end
   cj7 = ConstantRateJump(rate7, affect7!)
 
   jumpset = JumpSet(cj1,cj2,cj3,cj4,cj5,cj6,cj7)

   tspan = (0,tf)
   p = (kforward = parameters[1], kbackward = parameters[2], kon0 = parameters[3], kon1 = parameters[4],  koff = parameters[5], kinitiation = parameters[6], kdegradation = parameters[7]) 

   # build the model    
   dprob = DiscreteProblem(u0, tspan, p)
   jprob = JumpProblem(dprob, Direct(),  jumpset)

   sol = solve(jprob, SSAStepper())
   return sol
end




function simulate_bursttraj(parameters, tspan, timeres)
    # initial condition and time horizon
    u0 = [1.,0., 0.];
    burntime = 10000.;
    tf = burntime + tspan;
    
    # simulation
    sol = simulateModel(u0,tf,parameters);

    # read trajectories
    times=burntime+1:timeres:tf
    traj = sol(times,idxs=3)[1,:];
    
    return traj
end





end