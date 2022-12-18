# WHTS
Weighted hub and time sampling

### Guides
#### Quick start:
1) Running DeepID requires the [python](https://www.python.org/downloads/) (3.7 version or later) runtime environment; 
2) Make sure that extension package including [Numpy](https://numpy.org/), [Pandas](https://pandas.pydata.org/), networkx,sklearn, joblib, and tqdm have installed for current python environment; 
3) ```compare_top_multiply_p_sample_simple_test```  in the script file of ```Compare_nodes.py``` is the main function of this study, which contains the simulation process of Regular-1/7 and WHTS-1/f.
4) The input parameters of this function are nodenum (Nodes of the network), avg_degree (Mean degree of the network), beta (Transmission rate), gamma (Recovery time/day), FN (False-negative rate), s_p (sampling fraction of WHTS-1/f), infect_num (Number of imported cases), respectively.
5) The output is a DataFrame data format, including six columns of the number of testing when alerting (Regular-1/7), alerting time (Regular-1/7), the number of infections when alerting (Regular-1/7), the number of testing when alerting (WHTS-1/f), alerting time (WHTS-1/f), the number of infections when alerting (WHTS-1/f).
For example, we compared WHTS-1/7 with Regular-1/7 across 1000 repeats as follows:
```
outdir = 'D:/yourdir'
rlt = pd.DataFrame(columns=['norm_total','norm_days','norm_I','hub_total','hub_days','hub_I'])
for i in range(1000):
    rlt = rlt.append(compare_top_multiply_p_sample_simple_test(nodenum=50000,avg_degree=4,beta=0.1,gamma=7,FN=0.3,s_p = 1/7))
rlt.to_csv('%s/Rlt_50000_4_0.1_7_0.3_7.csv'%outdir)
```

3) Download the RLM.pdparams, SLM.pdparams, RLM.py and SLM.py, DeepID.py, test_repertoire_level_features.npy, test_sequence_level_features.npy and y_test.npy to the running directory;
4) The command for evaluating RLM on the test_repertoire_level_features is: 
   ```
   python RLM.py input_x_file input_true_Y output_file
   ```
   For example: 
   ```
   python RLM.py test_repertoire_level_features y_test RLM_test_rlt.csv
   ```
5) The command for evaluating SLM on the test_sequence_level_features is: 
   ```
   python SLM.py input_x_file input_true_Y output_file
   ```
   For example: 
   ```
   python SLM.py test_sequence_level_features y_test SLM_rlt_4feat.csv
   ```
6) The command for evaluating DeepID on the test set is (Must run RLM and SLM first):
   ```
   python DeepID.py input_RLM_rlt input_SLM_rlt output_file
   ```
   For example: 
   ```
   python DeepID.py test_rlt.csv rlt_4feat.csv DeepID_rlt.csv
   ```
7) The command for muti-classification on Healthy, HBV, influenza, and COVID-19:
   ```
   python muti-classification.py test_muti-classification y_test_muti-classification output_file
   ```

#### File details:
1) The RLM.pdparams and SLM.pdparams are the model files that have been trained;
2) RLM.py, SLM.py and DeepID.py are the scripts for RLM, SLM and DeepID model;
3) The test_repertoire_level_features.npy and test_sequence_level_features.npy are the 547 repertoire-level features and 160 sequence-level features for test dataset, respectively. The names and order of these features are listed in the Feature names.xlsx. test_repertoire_level_features.npy is a 120\*547 matrix with 120 samples and 547 features; test_sequence_level_features.npy is a 120\*160\*160 matrix, in which the three dimensions are samples, clones and sequence_level_features, respectively;
4) The y_test.npy is the true labels of the test dataset and is only used for accuracy calculation;
5) RLM_test_rlt.csv, SLM_rlt_4feat.csv and DeepID_rlt.csv are output files of the RLM.py, SLM.py and DeepID.py, respectively. The five columns of the CSV files are Probability of 0 (infection), Probability of 1 (healthy), samples are predicted to be class 0, samples are predicted to be class 1 and the true labels.
6) The user also can apply the DeepID to their test dataset by replacing the input data (test_repertoire_level_features.npy, test_sequence_level_features.npy and y_test.npy).
7）For the muti-classification, 0：COVID-19; 1:influenza; 2->HBV; 3:Healthy.
    
### Source: 
[Chen Y, Ye Z, Zhang Y, et al. A Deep Learning Model for Accurate Diagnosis of Infection Using Antibody Repertoires[J]. The Journal of Immunology,20222,08 (12): 2675–2685.](https://doi.org/10.4049/jimmunol.2200063)

If you have any questions or problems, please e-mail chenyuan0510@126.com
