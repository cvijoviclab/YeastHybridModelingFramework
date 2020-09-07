ReadMe analyseEnzymeUsage

analyzeEnzymeUsage explores the EnzymeUsage for the unregulated model and 
proteomics datasets to find a resonable level of regulation. Regulation in 
the multiscale model is done by calculating the parameter space and the 
allowed range is then altered by a percentage decided by the regulation 
parameter trough constraining the ub/lb depending on up/down regulation. 
Since moast of the regulated proteins work at their extreme bounds, we 
do expect an effect by this approach.

1)  First I look at the simulated data to see how the model behave.
2)  I use two proteomic datasets to compare to the simulated data. The 
    first dataset is for respiration taken at D_1 using chemostats: 
    "Lawless, Craig, et al. "Direct and absolute quantification of over 
    1800 yeast proteins via selected reaction monitoring." Molecular & 
    Cellular Proteomics 15.4 (2016): 1309-1322." Unfortunatly the dataset 
    and the proteins in the model had very low overlap (only about a third 
    of the proteins are also in the proteomics dataset). However, the 
    simulated data is set to a chemostat setting and teh dilutionrate D_1 
    was used to compare with the proteomics data. The second dataset is 
    from a batch experiment harwested at OD=0.6 wich I use to compare with 
    the maximum simulated growthrate for the model at 0.4. 

    To compare the dataset to the simulation I first extract the overlap 
    between the model and the proteomics dataset and normalize the value 
    to percentages (X/sum(protein in dataset overlap)). To vizulize it the 
    simulated data is plotted pairvise to the proteomics data in a scatter 
    plot and coulord based on regulation (red for downregulated proteins 
    and green for upregulated proteins) and the Bhattacharyya distance is 
    calculated. Even though the fit is inherrently bad the reasoning is to 
    use the regulation to get the model to use more realistic values for 
    the regulated enzymes and then adapt the fluxes and other enzymes 
    accordingly. Thats why the regulation parameter is decided on the 
    mean difference between the dataset and the simulated data.

    To evaluate the change the plot is made and the Bhattacharyya distance 
    is calculated for the regulated model as well. To set the regulated 
    model the parameter path2Sim2 in the loadData function is set to the 
    result folder for the correct simulation.