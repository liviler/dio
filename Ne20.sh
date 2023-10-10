#!/bin/bash
ELE=Ne # ZKr #Ca #Mo # Ti # Cr # Fe # Ni # Zn #Ge # Se # Ti62 # Cr64 # Fe66 # Ni68 # Zn70 # # Ge72 # Se74 # Kr #Sr # Zr # Mg # Si # Si # Mg
A=20
NUC=${ELE}${A}   # the DIC code generates the files with nucleus name as element name

# this is the same as that in compile_script.sh 
#export pathstore=/home/research/RPG/efzhou/data/${NUC} 
# export pathwfs=${pathstore}/wfs
# export GCM_FILES_DIR=${pathwfs}

export pathwork=./ #~/RPG/DIO/${NUC}
pathexe=${pathwork}/exe
mkdir -p ${pathexe}/input
mkdir -p ${pathexe}/output
#--------------------------------------------------------------------
#      generate input files: dic.dat, gcm.dat, betgam.dat
#--------------------------------------------------------------------
# parameter sets that could be used in the beyond mean-field calculations
# Force    =  PC,PC-F1               ! Parameterset of the Lagrangian
# Force    =  PC,PC-PK1              ! Parameterset of the Lagrangian
# Force    =  PC,DD-PC1              ! Parameterset of the Lagrangian
#V0       =  349.500   330.000      ! pairing strength for delta-force
#--------------------------------------------------------------------
cd ${pathexe}/input
echo create input file ....
cat <<EOF > dio.dat 
n0f,n0b  =  8  8                    ! number of oscillator shells
b0       = -2.448                   ! oscillator parameter (fm**-1) of basis
beta0    =  0.00                    ! deformation parameter of basis
betas    =  0.50                    ! deformation beta2 of W.S. potential
bet3s    =  0.00                    ! deformation beta3 of W.S. potential
maxi     =  400                     ! maximal number of iterations
xmix     =  0.50                    ! mixing parameter
inin     =  1                       ! 1 (calc. from beginning); 0 (read saved pot.) 
${ELE} $A                           ! nucleus under consideration
Ide      =  4  4                    ! Pairing control: 1. no  2. Frozen  3. G   4. Delta
Delta    =  0.000000  0.000000      ! Frozen Gaps (neutrons and protons)
Ga       =  0.000000  0.000000      ! Pairing-Constants GG = GA/A
Delta-0  =  2.000000  2.000000      ! Initial values for the Gaps
Force    =  PC-PK1
Vpair    =  349.500   330.000       ! pairing strength for delta force
icstr    =  2                       ! Quadratic constraint (no 0; beta2 1; b2+b3 2)
cspr     =  10.00                   ! Spring constant
cmax     =  1.000                   ! cutoff for dE/db
iwf      =  1                       ! 
c-------------------------------------------------------------------
EOF
#-----------------------------------
# mesh points in deformation q-space
####################################  
#  beta  gamma(deg)
####################################  
#  Number of mesh points in q-space 
#  is specified in gcm.dat file
####################################  
cat <<EOF > b23.dat 
    0.20  0.00   0.00
    0.40  0.00   0.00
    0.60  0.00   0.00
    0.80  0.00   0.00
    1.00  0.00   0.00
    1.20  0.00   0.00
    1.40  0.00   0.00
    0.20  0.20   0.00
    0.40  0.20   0.00
    0.60  0.20   0.00
    0.80  0.20   0.00
    1.00  0.20   0.00
    1.20  0.20   0.00
    1.40  0.20   0.00
    0.20  0.40   0.00
    0.40  0.40   0.00
    0.60  0.40   0.00
    0.80  0.40   0.00
    1.00  0.40   0.00
    1.20  0.40   0.00
    1.40  0.40   0.00
    0.20  0.60   0.00
    0.40  0.60   0.00
    0.60  0.60   0.00
    0.80  0.60   0.00
    1.00  0.60   0.00
    1.20  0.60   0.00
    1.40  0.60   0.00
    0.20  0.80   0.00
    0.40  0.80   0.00
    0.60  0.80   0.00
    0.80  0.80   0.00
    1.00  0.80   0.00
    1.20  0.80   0.00
    1.40  0.80   0.00
    0.20  1.00   0.00
    0.40  1.00   0.00
    0.60  1.00   0.00
    0.80  1.00   0.00
    1.00  1.00   0.00
    1.20  1.00   0.00
    1.40  1.00   0.00
    0.20  1.20   0.00
    0.40  1.20   0.00
    0.60  1.20   0.00
    0.80  1.20   0.00
    1.00  1.20   0.00
    1.20  1.20   0.00
    1.40  1.20   0.00
EOF

cd ../
echo -e "\033[32m run dio...\033[0m"
./dio

echo calculation is finished !

echo Done!
