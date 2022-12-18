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
    
 6) The script file of ```Compare_nodes.py``` is an example of parallel acceleration on the linux environment.

    
### Source: 
A more cost-saving or sensitive and robust strategy for regular nucleic acid-testing in COVID-19 pandemic.

If you have any questions or problems, please e-mail chenyuan0510@126.com
