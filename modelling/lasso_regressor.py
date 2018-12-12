# PROTOSPACER PRECISION LASSO
# Tristan Henser-Brownhill, September, 2018
# For Chakrabarti et al.

# IMPORTS
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.linear_model import LassoCV, Lasso
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, median_absolute_error
from sklearn.externals import joblib
from scipy.stats.stats import pearsonr, spearmanr, linregress, kendalltau # scipy v. 1.0.0
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from time import time
from scipy.stats import randint as sp_randint

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

# Note label order encoded ['A','C','G','T']
df_train_x = df_train_x.apply(LabelEncoder().fit_transform)
df_test_x = df_test_x.apply(LabelEncoder().fit_transform)
df_ob_x = df_ob_x.apply(LabelEncoder().fit_transform)
enc = OneHotEncoder(sparse = False)
enc.fit(df_train_x)
x_train = enc.transform(df_train_x)
x_test = enc.transform(df_test_x)
x_ob = enc.transform(df_ob_x)
'''
#MAIN MODEL SETUP

# parameter setup and model fitting
las = LassoCV(alphas = np.logspace(-5, 2, 30), cv = 10, verbose = 1, n_jobs = 4).fit(x_train, y_train) #fit the pixel values (x) to the labels (y)   

print("Chosen alpha:", las.alpha_)

# save model
joblib.dump(las, os.path.join('saved_models','lasso_model.pkl'))
'''

las = joblib.load(os.path.join('saved_models','lasso_model.pkl'))       
# MAKE PREDICTIONS AND CHECK ACCURACY

# get predictions
train_labels = y_train
train_pred = las.predict(x_train)
test_labels = y_test
test_pred = las.predict(x_test)
ob_labels = y_ob
ob_pred = las.predict(x_ob)

train_slope, train_intercept, train_r_value, train_p_value, train_std = linregress(train_labels, train_pred)
test_slope, test_intercept, test_r_value, test_p_value, test_std = linregress(test_labels, test_pred)
ob_slope, ob_intercept, ob_r_value, ob_p_value, ob_std = linregress(ob_labels, ob_pred)

# print metrics

with open('lasso_stats.csv', "w") as csv_file:

    csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}'.format("SET,","MSE,","MAE,","MAD,","LR_slope,","LR_intercept,","LR_r_val,","LR_wald_p_val,","LR_stdev\n"))
    csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}'.format("TRAIN,",mean_squared_error(train_labels, train_pred),",",mean_absolute_error(train_labels, train_pred),",",median_absolute_error(train_labels, train_pred),",",train_slope,",",train_intercept,",",train_r_value,",",train_p_value,",",train_std,"\n"))
    csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}{16}'.format("TEST,",mean_squared_error(test_labels, test_pred),",",mean_absolute_error(test_labels, test_pred),",",median_absolute_error(test_labels, test_pred),",",test_slope,",",test_intercept,",",test_r_value,",",test_p_value,",",test_std,"\n"))
    csv_file.write('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}'.format("OVERBEEK,",mean_squared_error(ob_labels, ob_pred),",",mean_absolute_error(ob_labels, ob_pred),",",median_absolute_error(ob_labels, ob_pred),",",ob_slope,",",ob_intercept,",",ob_r_value,",",ob_p_value,",",ob_std))
    
print('TRAINING SET RESULTS:')
print('MSE: ', mean_squared_error(train_labels, train_pred))
print('MAE: ', mean_absolute_error(train_labels, train_pred))
print('Median AE: ', median_absolute_error(train_labels, train_pred))
print('LR:', linregress(train_labels, train_pred))
print('Pearson:', pearsonr(train_labels, train_pred))
print('Spearman:', spearmanr(train_labels, train_pred))
print('Kendall:', kendalltau(train_labels, train_pred))

print('\nTEST SET RESULTS:')
print('MSE: ', mean_squared_error(test_labels, test_pred))
print('MAE: ', mean_absolute_error(test_labels, test_pred))
print('Median AE: ', median_absolute_error(test_labels, test_pred))
print('LR:', linregress(test_labels, test_pred))
print('Pearson:', pearsonr(test_labels, test_pred))
print('Spearman:', spearmanr(test_labels, test_pred))
print('Kendall:', kendalltau(test_labels, test_pred))

print('\nOVERBEEK RESULTS:')
print('MSE: ', mean_squared_error(ob_labels, ob_pred))
print('MAE: ', mean_absolute_error(ob_labels, ob_pred))
print('Median AE: ', median_absolute_error(ob_labels, ob_pred))
print('LR:', linregress(ob_labels, ob_pred))
print('Pearson:', pearsonr(ob_labels, ob_pred))
print('Spearman:', spearmanr(ob_labels, ob_pred))
print('Kendall:', kendalltau(ob_labels, ob_pred))


# FEATURE IMPORTANCES

importances = las.coef_
nt_importance = []
allAs = []
allCs = []
allGs = []
allTs = []
As = []
Cs = []
Gs = []
Ts = []
negAs = []
negCs = []
negGs = []
negTs = []
# output |summed| coeffs
for i in range(80):
    if ((i+1)%4 == 0):
        # sum of the abs of coeffs from the individual one-hot encoded features
        nt_importance.append(abs(importances[i-3])+abs(importances[i-2])+abs(importances[i-1])+abs(importances[i]))
        
        allAs.append(importances[i-3])
        if importances[i-3] < 0:
            As.append(0)
            negAs.append(importances[i-3])
        else:
            As.append(importances[i-3])
            negAs.append(0)
        
        allCs.append(importances[i-2])
        if importances[i-2] < 0:
            Cs.append(0)
            negCs.append(importances[i-2])
        else:
            Cs.append(importances[i-2])
            negCs.append(0)
        
        allGs.append(importances[i-1])
        if importances[i-1] < 0:
            Gs.append(0)
            negGs.append(importances[i-1])
        else:
            Gs.append(importances[i-1])
            negGs.append(0)
            
        allTs.append(importances[i])
        if importances[i] < 0:
            Ts.append(0)
            negTs.append(importances[i])
        else:
            Ts.append(importances[i])
            negTs.append(0) 
            
            
        # output table of individual nt letter at nt pos contributions
nt_importance.append(abs(importances[80])+abs(importances[81])+abs(importances[82])+abs(importances[83]))

allAs.append(importances[80])
if importances[80] < 0:
    As.append(0)
    negAs.append(importances[80])
else:
    As.append(importances[80])
    negAs.append(0)

allCs.append(importances[81])
if importances[81] < 0:
    Cs.append(0)
    negCs.append(importances[81])
else:
    Cs.append(importances[81])
    negCs.append(0)

allGs.append(importances[82])
if importances[82] < 0:
    Gs.append(0)
    negGs.append(importances[82])
else:
    Gs.append(importances[82])
    negGs.append(0)
    
allTs.append(importances[83])
if importances[83] < 0:
    Ts.append(0)
    negTs.append(importances[83])
else:
    Ts.append(importances[83])
    negTs.append(0) 

nt_importance.append(abs(importances[84]))
allAs.append(0), allCs.append(0), allGs.append(0), allTs.append(0)
nt_importance.append(abs(importances[85]))
allAs.append(0), allCs.append(0), allGs.append(0), allTs.append(0)

# DRAW PLOTS

#feature coefficients vs alphas
alphas = las.alphas_
coefs = []
for a in alphas:
    check_las = Lasso(alpha=a)
    check_las.fit(x_train, y_train)
    coefs.append(check_las.coef_)

ax = plt.gca()
ax.plot(alphas, coefs)
ax.set_xscale('log')
ax.set_xlim(ax.get_xlim()[::-1])  # reverse axis
plt.xlabel('alpha', fontsize = 18)
plt.ylabel('coefficients',fontsize = 18)
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
plt.title('LASSO coefficients as a function alpha', fontsize = 18)
plt.axis('tight')
plt.tight_layout()
plt.savefig('lasso_coeff_lambda.pdf', dpi = 300)


# Lasso model selection
m_log_alphas = -np.log10(las.alphas_)
plt.figure()
plt.plot(m_log_alphas, las.mse_path_, ':')
plt.plot(m_log_alphas, las.mse_path_.mean(axis=-1), 'k',
         label='Average across the folds', linewidth=2)
plt.axvline(-np.log10(las.alpha_), linestyle='--', color='k',
            label='alpha: CV estimate')
plt.legend()
plt.xlabel('-log(alpha)', fontsize = 18)
plt.ylabel('MSE', fontsize = 18)
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
plt.title('LASSO model selection', fontsize = 18)
plt.axis('tight')
plt.tight_layout()
plt.savefig('lasso_model_selection_cv.pdf', dpi = 300)

# plot the feature importance of LASSO (all)

plt.figure()
ind = np.arange(21)    # the x locations for the groups
width = 0.8       # the width of the bars: can also be len(x) sequence
p1 = plt.bar(ind, As, width, color = 'green')
p2 = plt.bar(ind, Ts, width, bottom = As, color = 'red')
p3 = plt.bar(ind, Cs, width, bottom = np.add(As, Ts), color = 'blue')
p4 = plt.bar(ind, Gs, width, bottom = np.add(Cs, np.add(As, Ts)), color = 'orange')
p5 = plt.bar(ind, negAs, width, color = 'green')
p6 = plt.bar(ind, negTs, width, bottom = negAs, color = 'red')
p7 = plt.bar(ind, negCs, width, bottom = np.add(negAs, negTs), color = 'blue')
p8 = plt.bar(ind, negGs, width, bottom = np.add(negCs, np.add(negAs, negTs)), color = 'orange')
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A', 'T', 'C', 'G'), fontsize  = 14, loc = 2)
plt.title('LASSO nucleotide importances', fontsize = 18)
plt.ylabel('coefficients', fontsize = 18)
plt.ylim([-0.1,0.2])
plt.yticks(fontsize = 14)
plt.xticks(range(23), ['-20', '-19', '-18', '-17', '-16', '-15', '-14', '-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'N', 'G', 'G'], rotation = 60, fontsize = 14)
plt.tight_layout()
plt.savefig('lasso_coeff.pdf', dpi = 300)

# plot the feature importances of LASSO (|summed| coeff)
plt.figure()
plt.title('LASSO nucleotide importances', fontsize = 18)
plt.bar(range(23), nt_importance, width = 0.8, color='gray')
plt.ylabel('absolute coefficients', fontsize = 18)
plt.yticks(fontsize = 12)
plt.xticks(range(23), ['-20', '-19', '-18', '-17', '-16', '-15', '-14', '-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'N', 'G', 'G'], rotation = 60, fontsize = 14)
plt.tight_layout()
plt.savefig('summed_lasso_coeff.pdf', dpi = 300)

##################


# Setting the positions and width for the bars
pos = list(range(len(allAs))) 
width = 0.25 
x_labels = ['-20', '-19', '-18', '-17', '-16', '-15', '-14', '-13', '-12', '-11', '-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', 'N', 'G', 'G']
    
# Plotting the bars
fig, ax = plt.subplots(figsize=(10,5))

# Create a bar with pre_score data,
# in position pos,
plt.bar(pos, 
        #using df['pre_score'] data,
        allAs, 
        # of width
        width, 
        # with color
        color='green') 

# Create a bar with mid_score data,
# in position pos + some width buffer,
plt.bar([p + width for p in pos], 
        #using df['mid_score'] data,
        allTs,
        # of width
        width,
        # with color
        color='red') 

# Create a bar with post_score data,
# in position pos + some width buffer,
plt.bar([p + width*2 for p in pos], 
        #using df['post_score'] data,
        allCs, 
        # of width
        width,
        # with color
        color='blue') 
        
plt.bar([p + width*3 for p in pos], 
        #using df['post_score'] data,
        allGs, 
        # of width
        width,
        # with color
        color='orange') 


# Set the y axis label
ax.set_ylabel('coefficient', fontsize = 20)

# Set the chart's title
ax.set_title('Nucleotide importance (LASSO)', fontsize = 20)

# Set the labels for the x ticks
plt.xticks(range(23), x_labels, rotation = 60, fontsize = 14)

# Set the position of the x ticks
ax.set_xticks([p + 2 * width for p in pos])

# Setting the x-axis and y-axis limits
plt.xlim(min(pos)-width, max(pos)+width*4)

# Adding the legend and showing the plot
plt.legend(['A', 'T', 'C', 'G'], loc='upper right', framealpha = 0, fontsize = 14)

for p in range(22):   
    plt.axvline(p + width*3.5, color = 'black', linewidth = 0.3)

plt.savefig('trail_lasso_coeff.pdf', dpi = 600)


##################

# Plot pred vs observed (test)

plot_text = '{0}{1:.2f}{2}{3}{4:.2E}'.format('R-value: ', test_r_value,'\n','P-value: ', test_p_value)

fig2, ax2 = plt.subplots()
ax2.plot(test_labels, test_pred,'o')
ax2.plot(test_labels, test_intercept + test_slope*test_labels, 'r')
ax2.text(0.3,0.85,plot_text,fontsize=18)
ax2.set_title('Unseen test_data\n', fontsize = 18)
ax2.set_ylim([0.0,1.0])
ax2.set_xlim([0.0,1.0])
ax2.tick_params('both', labelsize = 18)
ax2.set_ylabel('Predicted precision', fontsize = 18)
ax2.set_xlabel('Observed precision', fontsize = 18)
ax2.set_aspect('equal', 'box')
fig2.tight_layout()
plt.savefig('lasso_predictions_vs_labels_test.pdf', dpi = 300)

# Plot pred vs observed (Overbeek)

plot_text = '{0}{1:.2f}{2}{3}{4:.2E}'.format('R-value: ', ob_r_value,'\n','P-value: ', ob_p_value)

fig2, ax2 = plt.subplots()
ax2.plot(ob_labels, ob_pred,'o', color = 'y')
ax2.plot(ob_labels, ob_intercept + ob_slope*ob_labels, 'r')
ax2.text(0.3,0.85,plot_text,fontsize=18)
ax2.set_title('HCT116 (van Overbeek et al)\n', fontsize = 18)
ax2.set_ylim([0.0,1.0])
ax2.set_xlim([0.0,1.0])
ax2.tick_params('both', labelsize = 18)
ax2.set_ylabel('Predicted precision', fontsize = 18)
ax2.set_xlabel('Observed precision', fontsize = 18)
ax2.set_aspect('equal', 'box')
fig2.tight_layout()
plt.savefig('lasso_predictions_vs_labels_overbeek.pdf', dpi = 300)