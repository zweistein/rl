KIN="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/rlkin"
MDL="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/rlmdl"
PLN="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/rlplan"
SG="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/rlsg"


KIN_="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/reactiveTestFiles/rlkin"
MDL_="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/reactiveTestFiles/rlmdl"
PLN_="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/reactiveTestFiles/rlplan"
SG_="/home/elod/SoMa/ws_reactivePlanning/rl/rl-examples-0.6.2/reactiveTestFiles/rlsg"

#ln -s from to

#rm $SG/box-3d-narrow-V_.xml $SG/narrow-V_.wrl $SG/box-3d-narrow-V_.wrl  $PLN/box-3d-narrow-V_pcrrt_.xml
rm $SG/box-2d-narrow-V_.xml $SG/narrow-V_.wrl $SG/box-2d-narrow-V_.wrl  $PLN/box-2d-narrow-V_pcrrt_.xml #$KIN/box-2d-025025025.xml

#ln -s $PLN_/box-3d-narrow-V_pcrrt.xml $PLN/box-3d-narrow-V_pcrrt_.xml
#ln -s $SG_/box-3d-narrow-V.wrl $SG/box-3d-narrow-V_.wrl
#ln -s $SG_/box-3d-narrow-V.xml $SG/box-3d-narrow-V_.xml
#ln -s $SG_/narrow-V.wrl $SG/narrow-V_.wrl

ln -s $PLN_/box-2d-narrow-V_pcrrt.xml $PLN/box-2d-narrow-V_pcrrt_.xml
ln -s $SG_/box-2d-narrow-V.wrl $SG/box-2d-narrow-V_.wrl
ln -s $SG_/box-2d-narrow-V.xml $SG/box-2d-narrow-V_.xml
ln -s $SG_/cones.wrl $SG/narrow-V_.wrl
#ln -s $SG_/V_45_0.4.wrl $SG/narrow-V_.wrl
#ln -s $SG_/V_45_0.2repaired.wrl $SG/narrow-V_.wrl

#ln -s $KIN_/box-2d-025025025_.xml $KIN/box-2d-025025025.xml

