#!/bin/bash

simname=$1
model=$3


max_sessions=60
session_prefix="sim_m"

N=$2
arr=()
for ((i=1; i<=N; i++)); do
  arr+=("$i")
done


cp /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/config_generic.jl /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/config.jl
sed -i 's/SIMNAMEHERE/'$simname'/g' /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/config.jl
sed -i 's/FITNAMEHERE/'$simname'/g' /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/config.jl
sed -i 's/MODELID/'$model'/g' /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/config.jl

mkdir /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/bestfits/$simname


for i in "${arr[@]}"
do
    simid="${simname}_${i}"

    # Wait until fewer than max_sessions tmux sessions are running
    while [ "$(tmux ls 2>/dev/null | grep -c "^$session_prefix")" -ge "$max_sessions" ]; do
        echo "Waiting for slots to be free..."
        sleep 5  # Wait 5 seconds before checking again
        # List all matching tmux sessions
        for session in $(tmux ls 2>/dev/null | awk -F: -v prefix="$session_prefix" '$1 ~ "^"prefix {print $1}'); do
            # Check if the session is still running something
            if ! tmux list-panes -t "$session" -F "#{pane_active}" | grep -q 1; then
                echo "Killing inactive session: $session"
                tmux kill-session -t "$session"
            fi
        done
    done

    # Start a new tmux session in detached mode
    echo "Launching $i"

    cp /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/fit_expmodel_generic.jl /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/runfit_${i}.jl
    sed -i 's/SIMIDX/'$i'/g' /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/runfit_${i}.jl

    TMUX= tmux new-session -d -s $simid "julia --project /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/fitsim/runfit_${i}.jl"
done


  while [ "$(tmux ls 2>/dev/null | grep -c "^$session_prefix")" -gt 0 ]; do
      echo "Waiting job to finish..."
      echo "$session_prefix"
      sleep 5  # Wait 5 seconds before checking again
      # List all matching tmux sessions
  done

sleep 1  # Simulate long task

# Signal completion
touch /tungstenfs/scratch/ggiorget/Gregory/julia/BurstingDynamicsInference/src/modelassumption_robustness/myscript_done.flag


