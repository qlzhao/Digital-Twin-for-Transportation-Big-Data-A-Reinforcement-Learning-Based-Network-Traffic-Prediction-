This project is part of the code in the paper (DOI: 10.1109/TITS.2022.3232518) entitled Digital Twin for Transportation Big Data: A Reinforcement Learning-Based Network Traffic Prediction Approach.
The neural network's are mainly trained for the Abilene network, the hyperparameters for training other networks may vary.

A list of all documents is described below:

main_GAN_DQN_TM.m : This file is the main program of our algorithm, in which we build the neural network, train it with the training set, and finally make predictions on the test set data.

function_Train_D.m : This file describes the architecture of the D-neural network and information on the training parameters.

function_Train_G.m : This file describes the architecture of the G-neural network and information on the training parameters.

function_IPFP.m : Use the ipfp.p algorithm for each moment of the incoming matrix.

ipfp.p : The details of the ipfp algorithm recorded therein are used to make corrections to the predicted traffic size.

function_SRE.m : Calculating the spatial relative error of the traffic size matrix.

function_TRE.m : Calculating the temporal relative error of the traffic size matrix.

function_SRE_CDF.m : Calculating the spatial relative error probability density function statistics for the traffic size matrix.

function_TRE.m : Calculating the temporal relative error probability density function statistics for the traffic size matrix.

display_OD_SRE_TRE_CDF_ALL3.m : This file calculates the probability density functions for temporal relative error, spatial relative error, and error, and then plots them.

display_bias_SD_Abilene_GEANT3.m : Calculate the deviation and standard deviation information for the three methods in the three networks.

display_TM_vs_Estimates1.m : Visualisation of our method's predicted traffic size matrix against real traffic data. Example 1 in Abilene. 

display_TM_vs_Estimates2.m : Visualisation of our method's predicted traffic size matrix against real traffic data. Example 2 in GEANT. 

display_TM_vs_Estimates3.m : Visualisation of our method's predicted traffic size matrix against real traffic data. Example 3 in IoT.

function_Pre_Abilene.m : Pre-processing in Abilene network comparison of the improvement rate of our method with PCA and SRMF methods.

function_Pre_GEANT.m : Pre-processing in GEANT network comparison of the improvement rate of our method with PCA and SRMF methods.

function_Pre_IoT.m : Pre-processing in IoT network comparison of the improvement rate of our method with PCA and SRMF methods.

display_improvement2.m : Graphing the improved rate ratio obtained from the function shows.

Comparison of hyperparameters between networks trained for reinforcement learning:

|                           | Abilene | GEANT |  IoT  |
| :-----------------------: | :-----: | :---: | :---: |
|      length_training      |   500   |  500  |  500  |
|        length_data        |  2016   | 3360  | 2016  |
|         precision         |    5    |   3   |   2   |
| K(Number of updates of D) |   100   |   1   |  100  |
|           num_G           |  1000   |   1   |  500  |
|        maxEpoch_G         |    1    |   1   |   1   |
|      miniBatchSize_G      |  15000  |  50   |  500  |
|    InitialLearnRate_G     |  0.001  | 0.001 | 0.001 |
|        maxEpoch_D         |    1    |   1   |   1   |
|      miniBatchSize_D      |  15000  |  50   |  500  |
|    InitialLearnRate_D     |  0.001  | 0.001 | 0.001 |
|     epsilon_threshold     |   0.5   |  0.5  |  0.5  |
|             C             |    5    |   1   |  10   |
|         Iteration         |   100   |   1   |  100  |



