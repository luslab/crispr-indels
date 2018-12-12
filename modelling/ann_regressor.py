# PROTOSPACER PRECISION ANN
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

# import
data = pd.read_csv(os.path.join('data','random_train_set80pc.csv'))
data_test = pd.read_csv(os.path.join('data','random_test_set20pc.csv'))
#data_ob = pd.read_csv(os.path.join('data','ob_all48h.csv'))
data_ob = pd.read_csv(os.path.join('data','ob_HCT116_48h.csv'))
#data_ob = pd.read_csv(os.path.join('data','ob_HEK293_48h.csv'))
#data_ob = pd.read_csv(os.path.join('data','ob_K562_48h.csv'))

# prepare the data
df_train_x = data.iloc[:,2:] #extract pixels; all rows, all except sgRna name (0) and label (1)
y_train = data.iloc[:,1].values #extract labels
df_test_x = data_test.iloc[:,2:]
y_test = data_test.iloc[:,1].values
df_ob_x = data_ob.iloc[:,2:]
y_ob = data_ob.iloc[:,1].values

# one hot encoding

# Note label encoded order ['A','C','G','T']
df_train_x = df_train_x.apply(LabelEncoder().fit_transform)
df_test_x = df_test_x.apply(LabelEncoder().fit_transform)
df_ob_x = df_ob_x.apply(LabelEncoder().fit_transform)
enc = OneHotEncoder(sparse = False)
enc.fit(df_train_x)
x_train = enc.transform(df_train_x)
x_test = enc.transform(df_test_x)
x_ob = enc.transform(df_ob_x)

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


# DIRECT LOAD (or uncomment to bootstrap validate or train below as necessary)
net.load_params(os.path.join('saved_models','ANN_model.params'), ctx = ctx)


'''
# BOOTSTRAP VALIDATION

for i in range(10):

    # Build ANN
    net = gluon.nn.Sequential()
    with net.name_scope():
        net.add(gluon.nn.Dense(512, activation = "relu"))
        net.add(gluon.nn.Dense(1, activation = "softrelu"))

    # TRAINING    
        
    # Initialize ANN
    net.collect_params().initialize(mx.init.Xavier(magnitude = 2.24), ctx = ctx)

    # Training setup
    L2_loss = gluon.loss.L2Loss() #LOSS
    trainer = gluon.Trainer(net.collect_params(), 'nag', {'learning_rate': 0.001, 'momentum': 0.9})

    # Training loop
    epochs = 1500
    learning_curve_acc = []
    learning_curve_val_acc = []
    learning_curve_loss = []

    boot_train_x, boot_val_x, boot_train_y, boot_val_y = train_test_split(x_train, y_train, test_size = 0.2, random_state = 42, shuffle = True)

    boot_train_data = mx.io.NDArrayIter(boot_train_x, label = boot_train_y, batch_size = 100, shuffle = True)
    boot_val_data = mx.nd.array(boot_val_x).as_in_context(ctx)
    boot_val_labels = mx.nd.array(boot_val_y).as_in_context(ctx)
    boot_train_labels = mx.nd.array(boot_train_y).as_in_context(ctx) # nonbatch
    boot_nonbatch_train_data = mx.nd.array(boot_train_x).as_in_context(ctx) # nonbatch
    
    filename = os.path.join('data','ann_bootstrap','{0}{1}{2}'.format('ANN_bootstrap_validation_fold_',i,'.csv'))
    with open(filename, "w") as csv_file:
        for e in range(epochs):
            # Reset iterator
            boot_train_data.reset()
            # Loop the data iterator
            for batch in boot_train_data:
                # Split the data and labels into slices along the correct axis (batch_axis)
                data = gluon.utils.split_and_load(batch.data[0], ctx_list = ctx_list, batch_axis = 0)
                label = gluon.utils.split_and_load(batch.label[0], ctx_list = ctx_list, batch_axis = 0)
                outputs = []
                with autograd.record():
                    for x, y in zip(data, label):
                        z = net(x)
                        loss = L2_loss(z, y)
                        loss.backward()
                        z = nd.argmax(net(x), axis = 1)
                        outputs.append(z)
                trainer.step(batch.data[0].shape[0])
            avg_loss = loss.mean()
            avg_loss = avg_loss.asscalar()
            
            # Run ANN on train/test to check results at current epoch for learning curves
            
            train_pred = net(boot_nonbatch_train_data)
            train_pred = train_pred[:,0]
            mse = ((boot_train_labels - train_pred) ** 2).mean()
            mse = mse.asscalar()
            rmse = math.sqrt(mse)
            mae = (boot_train_labels - train_pred).mean()
            mbe = mae.asscalar()
            mae = abs(mae.asscalar())
            
            test_pred = net(boot_val_data)
            test_pred = test_pred[:,0]
            val_mse = ((boot_val_labels - test_pred) ** 2).mean()
            val_mse = val_mse.asscalar()
            val_rmse = math.sqrt(val_mse)
            val_mae = (boot_val_labels - test_pred).mean()
            val_mbe = val_mae.asscalar() # mean bias error give directional bias
            val_mae = abs(val_mae.asscalar())
            
            _, _, TR_R_VALUE, TR_LIN_P, _ = linregress(boot_train_labels.asnumpy(), train_pred.asnumpy())
            _, _, TE_R_VALUE, TE_LIN_P, _ = linregress(boot_val_labels.asnumpy(), test_pred.asnumpy())
            
            print('Epoch %d, RMSE = %f. Val RMSE = %f. MAE = %f. Val MAE = %f. Avg L2 loss = %f'% (e+1, rmse, val_rmse, mae, val_mae, avg_loss))
            
            csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}{17}{18}{19}{20}{21}{22}{23}{24}{25}{26}{27}{28}{29}'.format('EPOCH,', e+1,',','RMSE,',rmse,',','VAL_RMSE,',val_rmse,',','MAE,',mae,',','VAL_MAE,',val_mae,',','AVG L2 LOSS,',avg_loss,',','LINREG_R,',TR_R_VALUE,',','LINREG_P,',TR_LIN_P,',','VAL_LINREG_R,',TE_R_VALUE,',','VAL_LINREG_P,',TE_LIN_P,'\n'))
            
            learning_curve_acc.append(rmse)
            learning_curve_val_acc.append(val_rmse)
            learning_curve_loss.append(avg_loss)
    
    del net

# FINISH BOOTSTRAP
'''

'''
# FINAL MODEL TRAINING

# TRAINING    
    
# Initialize ANN
net.collect_params().initialize(mx.init.Xavier(magnitude = 2.24), ctx = ctx)

# Training setup
L2_loss = gluon.loss.L2Loss() #LOSS
trainer = gluon.Trainer(net.collect_params(), 'nag', {'learning_rate': 0.001, 'momentum': 0.9})

# Training loop
epochs = 800
learning_curve_acc = []
learning_curve_val_acc = []
learning_curve_loss = []

with open('model_training_stats.csv', "w") as csv_file:
    for e in range(epochs):
        # Reset iterator
        train_data.reset()
        # Loop the data iterator
        for batch in train_data:
            # Split the data and labels into slices along the correct axis (batch_axis)
            data = gluon.utils.split_and_load(batch.data[0], ctx_list = ctx_list, batch_axis = 0)
            label = gluon.utils.split_and_load(batch.label[0], ctx_list = ctx_list, batch_axis = 0)
            outputs = []
            with autograd.record():
                for x, y in zip(data, label):
                    z = net(x)
                    loss = L2_loss(z, y)
                    loss.backward()
                    z = nd.argmax(net(x), axis = 1)
                    outputs.append(z)
            trainer.step(batch.data[0].shape[0])
        avg_loss = loss.mean()
        avg_loss = avg_loss.asscalar()
        
        # Run ANN on train/test to check results at current epoch for learning curves
        
        train_labels = nonbatch_train_labels
        train_pred = net(nonbatch_train_data)
        train_pred = train_pred[:,0]
        mse = ((train_labels - train_pred) ** 2).mean()
        mse = mse.asscalar()
        rmse = math.sqrt(mse)
        mae = (train_labels - train_pred).mean()
        mbe = mae.asscalar()
        mae = abs(mae.asscalar())
        
        test_pred = net(test_data)
        test_pred = test_pred[:,0]
        val_mse = ((test_labels - test_pred) ** 2).mean()
        val_mse = val_mse.asscalar()
        val_rmse = math.sqrt(val_mse)
        val_mae = (test_labels - test_pred).mean()
        val_mbe = val_mae.asscalar() # mean bias error give directional bias
        val_mae = abs(val_mae.asscalar())
        
        _, _, TR_R_VALUE, TR_LIN_P, _ = linregress(train_labels.asnumpy(), train_pred.asnumpy())
        _, _, TE_R_VALUE, TE_LIN_P, _ = linregress(test_labels.asnumpy(), test_pred.asnumpy())
        
        print('Epoch %d, RMSE = %f. Val RMSE = %f. MAE = %f. Val MAE = %f. Avg L2 loss = %f'% (e+1, rmse, val_rmse, mae, val_mae, avg_loss))
        
        csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}{17}{18}{19}{20}{21}{22}{23}{24}{25}{26}{27}{28}{29}'.format('EPOCH,', e+1,',','RMSE,',rmse,',','VAL_RMSE,',val_rmse,',','MAE,',mae,',','VAL_MAE,',val_mae,',','AVG L2 LOSS,',avg_loss,',','LINREG_R,',TR_R_VALUE,',','LINREG_P,',TR_LIN_P,',','VAL_LINREG_R,',TE_R_VALUE,',','VAL_LINREG_P,',TE_LIN_P,'\n'))
        
        learning_curve_acc.append(rmse)
        learning_curve_val_acc.append(val_rmse)
        learning_curve_loss.append(avg_loss)
 
    
# Save weights        
net.save_params(os.path.join('saved_models','ANN_model.params'))
 
# Print learning curves
fig, ax1 = plt.subplots()
epoch_ax = np.arange(1, epochs+1, 1)
ax1.plot(epoch_ax, learning_curve_acc, 'b')
ax1.plot(epoch_ax, learning_curve_val_acc, 'r')
ax1.set_xlabel('Epoch', fontsize = 14)
ax1.set_ylabel('RMSE', fontsize = 14)
ax1.tick_params('both', labelsize = 12)
max_train_rmse = max(learning_curve_acc)
max_val_rmse = max(learning_curve_val_acc)
if (max_train_rmse > max_val_rmse):
    max_rmse = max_train_rmse
else:
    max_rmse = max_val_rmse
min_train_rmse = min(learning_curve_acc)
min_test_rmse = min(learning_curve_val_acc)
if (min_train_rmse < min_test_rmse):
    min_rmse = min_train_rmse
else:
    min_rmse = min_test_rmse
ax1.set_ylim([0, max_rmse+(max_rmse*0.1)])
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins = 10))
ax1.legend(['Training','Test'])
#ax2 = ax1.twinx()
#ax2.plot(epoch_ax, learning_curve_loss, 'g')
#ax2.set_ylabel('L2 loss', color = 'g')
#ax2.tick_params('y', colors = 'g')
#ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer = True))
#ax2.set_ylim([0, max(learning_curve_loss)])
#ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
fig.tight_layout()
plt.savefig('ANN_reg_train_curves_strap10.png')


# Stats

print('Train preds:')
print(train_labels[:5])
print(train_pred[:5])
print('Test preds:')
print(test_labels[:5])
print(test_pred[:5])

print('final training RMSE:', rmse)
print('final training MAE:', mae)
print('final training MBE:', mbe)
print('final val RMSE:', val_rmse)
print('final val MAE:', val_mae)
print('final val MBE:', val_mbe)


print('train LR:', linregress(train_labels.asnumpy(), train_pred.asnumpy()))
print('val LR:', linregress(test_labels.asnumpy(), test_pred.asnumpy()))
print('train Pearson:', pearsonr(train_labels.asnumpy(), train_pred.asnumpy()))
print('Val Pearson:', pearsonr(test_labels.asnumpy(), test_pred.asnumpy()))
print('train Spearman:', spearmanr(train_labels.asnumpy(), train_pred.asnumpy()))
print('Val Spearman:', spearmanr(test_labels.asnumpy(), test_pred.asnumpy()))
print('train Kendall:', kendalltau(train_labels.asnumpy(), train_pred.asnumpy()))
print('Val Kendall:', kendalltau(test_labels.asnumpy(), test_pred.asnumpy()))


# FINISH MODEL TRAINING
'''

# PREDICTIONS ONLY

test_pred = net(test_data)
test_pred = test_pred[:,0]

slope, intercept, r_value, p_value, std_err = linregress(test_labels.asnumpy(), test_pred.asnumpy())
pear_cor, pear_p = pearsonr(test_labels.asnumpy(), test_pred.asnumpy())
spear_cor, spear_p = spearmanr(test_labels.asnumpy(), test_pred.asnumpy())
ken_tau_b, ken_p = kendalltau(test_labels.asnumpy(), test_pred.asnumpy())
temp_mse = ((test_labels.asnumpy() - test_pred.asnumpy()) ** 2).mean()
temp_rmse = math.sqrt(temp_mse)
with open('prediction_stats.csv', "w") as csv_file:
    csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}{17}{18}{19}{20}{21}{22}{23}{24}'.format('RMSE,',temp_rmse,',R: ,', r_value,',', ' wald p: ,', p_value,',', ' pearson: ,', pear_cor,',', ' pear p: ,', pear_p,',', ' spearman: ,', spear_cor,',', ' spear p: ,', spear_p,',', ' kendall tau b: ,', ken_tau_b,',', ' ken p: ,', ken_p))


# PLOT

plot_text = '{0}{1:.2f}{2}{3}{4:.2E}'.format('R-value: ', r_value,'\n','P-value: ', p_value)

fig2, ax2 = plt.subplots()
ax2.plot(test_labels.asnumpy(), test_pred.asnumpy(),'o') #include for ob: color = 'y'
ax2.plot(test_labels.asnumpy(), intercept + slope*test_labels.asnumpy(), 'r')
ax2.text(0.3,0.85,plot_text,fontsize=18)
ax2.set_title('Unseen test_data\n', fontsize = 18) #repace for ob with: HCT116 (van Overbeek et al) or w/e
ax2.set_ylim([0.0,1.0])
ax2.set_xlim([0.0,1.0])
ax2.tick_params('both', labelsize = 18)
ax2.set_ylabel('Predicted precision', fontsize = 18)
ax2.set_xlabel('Observed precision', fontsize = 18)
ax2.set_aspect('equal', 'box')
fig2.tight_layout()
plt.savefig('ann_predictions_vs_observed.pdf', dpi = 600)
