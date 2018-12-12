# PROTOSPACER PRECISION ANN -  INDIVIDUAL NT PERMUTATION
# Tristan Henser-Brownhill, July, 2018
# For Chakrabarti et al.

# IMPORTS
import numpy as np # numpy v. 1.13.3
import math
import os 
import pandas as pd # pandas v. 0.22.0
from matplotlib import pyplot as plt # matplotlib v. 2.1.1
from matplotlib import ticker as ticker
from sklearn.model_selection import train_test_split # scikit-learn v. 0.19.1
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from scipy.stats.stats import pearsonr, spearmanr, linregress, kendalltau # scipy v. 1.0.0
import mxnet as mx # mxnet v. 1.2.0
from mxnet import nd, autograd, gluon
ctx = mx.cpu()
ctx_list = [mx.cpu()] #set mxnet ndarray context
mx.random.seed(1)

#DATA PREPROCESSING
for i in range(21):

    nt = '{0}{1}'.format('nt',i+1)
    filename = os.path.join('data','permutation','ann','{0}{1}{2}'.format('ANN_permutation_nt_',i+1,'_.csv')) # change this
    with open(filename, "w") as csv_file:
        # print first line
        csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}'.format("BOOTSTRAP,","RMSE,","MAE,","MBE,","LR_slope,","LR_intercept,","LR_r_val,","LR_wald_p_val,","LR_stdev\n"))
        
        for j in range(10):
        
            fil = os.path.join('data','random_train_set80pc.csv')
            fil_test = os.path.join('data','random_test_set20pc.csv')
            # import
            data = pd.read_csv(fil)
            data_test = pd.read_csv(fil_test)
            #data = pd.read_csv(os.path.join('data','random_train_set80pc-1.csv'))
            #data_test = pd.read_csv(os.path.join('data','random_test_set20pc.csv'))
            #data_ob = pd.read_csv(os.path.join('data','ob_all48h.csv'))
            #data_ob = pd.read_csv(os.path.join('data','ob_HCT116_48h.csv'))
            #data_ob = pd.read_csv(os.path.join('data','ob_HEK293_48h.csv'))
            #data_ob = pd.read_csv(os.path.join('data','ob_K562_48h.csv'))

            # prepare the data
            df_train_x = data.iloc[:,2:] #extract pixels; all rows, all except sgRna name (0) and label (1)
            y_train = data.iloc[:,1].values #extract labels
            df_test_x = data_test.iloc[:,2:]
            y_test = data_test.iloc[:,1].values
            #df_ob_x = data_ob.iloc[:,2:]
            #y_ob = data_ob.iloc[:,1].values
            
            # HERE WE PERMUTATE THE NUCLEOTIDES
            np.random.seed(j)
            #nt_choices = ['A','T','C','G']
            np.random.shuffle(df_test_x.iloc[:,i])
            
            # perform label encoding

            # Note label ancoded ['A','C','G','T']
            df_train_x = df_train_x.apply(LabelEncoder().fit_transform)
            df_test_x = df_test_x.apply(LabelEncoder().fit_transform)

            # one hot encoding
            #df_ob_x = df_ob_x.apply(LabelEncoder().fit_transform) 
            enc = OneHotEncoder(sparse = False)
            enc.fit(df_train_x)
            x_train = enc.transform(df_train_x)
            x_test = enc.transform(df_test_x)
            #x_ob = enc.transform(df_ob_x)

            # Convert to iterator data batches for ANN input
            train_data = mx.io.NDArrayIter(x_train, label = y_train, batch_size = 100, shuffle = True)
            test_data = mx.nd.array(x_test).as_in_context(ctx)
            test_labels = mx.nd.array(y_test).as_in_context(ctx)

            # Repackage the train data for checking in one go at the end! :not used by training!
            nonbatch_train_data = mx.nd.array(x_train).as_in_context(ctx)
            nonbatch_train_labels = mx.nd.array(y_train).as_in_context(ctx)
        
            # Build ANN
            net = gluon.nn.Sequential()
            with net.name_scope():
                net.add(gluon.nn.Dense(512, activation = "relu"))
                net.add(gluon.nn.Dense(1, activation = "softrelu"))
            
            # DIRECT LOAD PARAMS
            net.load_params(os.path.join('saved_models','ANN_model.params'), ctx = ctx)
                    
            test_pred = net(test_data)
            test_pred = test_pred[:,0]
            val_mse = ((test_labels - test_pred) ** 2).mean()
            val_mse = val_mse.asscalar()
            val_rmse = math.sqrt(val_mse)
            val_mae = (test_labels - test_pred).mean()
            val_mbe = val_mae.asscalar() # mean bias error give directional bias
            val_mae = abs(val_mae.asscalar())
            
            test_slope, test_intercept, TE_R_VALUE, TE_LIN_P, test_std = linregress(test_labels.asnumpy(), test_pred.asnumpy())
  
            print('nt',i+1,' bootstrap #',j,' completed.')
            if i+1 == 24:
                print('i.e. unaltered sequence')
            csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}{17}'.format(j,",",val_rmse,",",val_mae,",",val_mbe,",",test_slope,",",test_intercept,",",TE_R_VALUE,",",TE_LIN_P,",",test_std,"\n"))
            del net