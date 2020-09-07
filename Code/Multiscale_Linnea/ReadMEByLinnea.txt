ReadMeByLinnea

multiscaleSimulation2

Runs the multiscale simulation where the Boolean constrain the FBA model 
when a critical specific glucose uptake rate is reached the Boolean 
switches from glucose depleeted to glucose rich mode.
-   postTranscriptionalModifications2 implements seperate output for the 
    Boolean PTM.
-   constrainFBA2 implements seperate constraints on the FBA model where 
    generegulations is modified by enzymeUsage and PTM are modified by 
    kcat.

analyseEnzymeUsage

-   This scripts analyzes the enzyme usage from the multiscale 
    simulations and the proteomics data to see understand how the model 
    works and how it should be regulated. 
-   Using functions: mapEnzymeSubsystem
