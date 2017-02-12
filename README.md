# Final report branch

This is the version of the planner that I used to write the final report -> 80dce39

## Experiments

The following experiments were run for the benchmarks in Release build mode.

### RRT-Connect 
```
./demos/rlPlanDemo/rlPlanDemo ../rl-examples-0.6.2/rlplan/box-2d-gripper_rrtCon.xml
```

### Connect-move-only uncertainty planner  
```
./demos/rlPlanDemo/rlPlanDemo ../rl-examples-0.6.2/rlplan/box-2d-gripper_connect-pcrrt.xml
```

### Full uncertainty planner
```
./demos/rlPlanDemo/rlPlanDemo ../rl-examples-0.6.2/rlplan/box-2d-gripper_pcrrt.xml
```

### POMDP full uncertainty planner
```
./demos/rlPlanDemo/rlPlanDemo ../rl-examples-0.6.2/rlplan/box-2d-gripper-pomdp_pcrrt.xml
```
