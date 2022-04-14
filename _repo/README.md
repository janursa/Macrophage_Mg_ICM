# How to run calibration
The calibration aims at tuning the parameters of the model called Mg_M.xml to initially produce the same results as Zhao_2021.xml and secondly to match the observations. After changing the target parameters, first sampling should be done from Zhao_2021.xml, which is done by 'Python sample.py'. 
Then, the calibration should be called by 
'python diff_calibration.py'
The tuned values of the free parameters are stored to 'inferred_values.csv'