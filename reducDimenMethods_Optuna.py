# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
## ---------------------------
##
## Script name: all_in_one
##
## Purpose of script:
##
## Author: Diego Salazar
##
## Date Created: 2021-10-04
##
## Copyright (c) Diego Salazar, 2021
## Email: das4019@med.cornell.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

import optuna
import pandas as pd
# import numpy as np
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
# from sklearn.decomposition import KernelPCA
# from sklearn.metrics import mean_squared_error
# from sklearn.model_selection import GridSearchCV
# from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler
# from keras.models import Sequential
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import LeakyReLU
from tensorflow.keras.layers import BatchNormalization

# datos = pd.read_csv('D:/pathFeatureLabelKernels/co-mod_best_results_last_run/all_TCGA_LGG_GBM.csv')
# datos = pd.read_csv('D:/pathFeatureLabelKernels/co-md_records/O_Y_A_Delta_W.csv')
datos = pd.read_csv('D:/pathFeatureLabelKernels/co-md_records/imputation_1.csv')
# datos = datos.loc[datos['delta'] == 1,:].copy()
datos = datos.loc[datos['Y'] != 0,:].copy()
datos.reset_index(inplace=True)

mrna = datos.filter(regex="mrna_"); print(mrna.shape)
# mrna = datos.filter(regex="mirna_"); print(mrna.shape)
# cnv = datos.filter(regex="cnv_"); print(cnv.shape)

# %% Optimze NMF
# def objective(trial):
    
#     n_components = trial.suggest_int("n_components", 50,500)
#     # init = trial.suggest_categorical("init", ['random', 'nndsvd', 'nndsvda', 'nndsvdar'])
#     alpha = trial.suggest_uniform("alpha", 0.0, 1.0)
#     l1_ratio = trial.suggest_uniform("l1_ratio", 0.0, 1.0)
    
#     model = NMF(n_components=n_components, init='random',
#                 random_state=0, tol=1e-7, max_iter=500, 
#                 alpha=alpha, l1_ratio=l1_ratio)
#     model.fit(cnv.T)

#     error = model.reconstruction_err_
    
#     return error

# study = optuna.create_study(direction='minimize')
# study.optimize(objective, n_trials=100)

# ## cnv -> GBM + LGG
# # [I 2022-01-19 23:13:17,206] Trial 7 finished with value: 25.930487330221194 and parameters: 
#     # {'n_components': 468, 'alpha': 0.4031682277410855, 'l1_ratio': 0.12316508198716614}. 
#     # Best is trial 7 with value: 25.930487330221194.

# ## cnv -> OV
# # [I 2022-01-17 10:27:54,571] Trial 1 finished with value: 13.570008889443843 and parameters: 
#     # {'n_components': 298, 'alpha': 0.16606797571125143, 'l1_ratio': 0.37292885046544433}. 
#     # Best is trial 1 with value: 13.570008889443843.
    
# model = NMF(n_components=study.best_trial.params["n_components"], 
#             alpha=study.best_trial.params["alpha"], l1_ratio=study.best_trial.params["l1_ratio"],
#               init='random', random_state=0, tol=1e-7, max_iter=1000)
# model.fit(cnv.T)

# pd.DataFrame(model.components_.T).to_csv("D:/pathFeatureLabelKernels/co-md_records/nmf_W_cnv.csv")

# print("/\ /\ /\ /\ /\ /\    NMF_CNV    /\ /\ /\ /\ /\ /\ ")
# print("|| || || || || ||               || || || || || || ")
# print("                 Input shape: ")
# print("                             "+str(cnv.shape))
# print(" ........................................... ")

# %% Optimze NMF
def objective(trial):
    
    n_components = trial.suggest_int("n_components", 50,500)
    # init = trial.suggest_categorical("init", ['random', 'nndsvd', 'nndsvda', 'nndsvdar'])
    alpha = trial.suggest_uniform("alpha", 0.0, 1.0)
    l1_ratio = trial.suggest_uniform("l1_ratio", 0.0, 1.0)
    
    model = NMF(n_components=n_components, init='random',
                random_state=0, tol=1e-7, max_iter=500, 
                alpha=alpha, l1_ratio=l1_ratio)
    model.fit(mrna.T)

    error = model.reconstruction_err_
    
    return error

study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=100)

## mrna -> GBM + LGG
# [I 2022-01-19 23:13:17,206] Trial 7 finished with value: 25.930487330221194 and parameters: 
    # {'n_components': 468, 'alpha': 0.4031682277410855, 'l1_ratio': 0.12316508198716614}. 
    # Best is trial 7 with value: 25.930487330221194.

## mrna -> OV
# [I 2022-01-17 11:12:53,197] Trial 6 finished with value: 3.2592297824188363 and parameters: 
    # {'n_components': 332, 'alpha': 0.06623971974650289, 'l1_ratio': 0.356629166227025}. 
    # Best is trial 6 with value: 3.2592297824188363.

model = NMF(n_components=study.best_trial.params["n_components"], 
            alpha=study.best_trial.params["alpha"], l1_ratio=study.best_trial.params["l1_ratio"],
              init='random', random_state=0, tol=1e-7, max_iter=1000)
model.fit(mrna.T)

pd.DataFrame(model.components_.T).to_csv("D:/pathFeatureLabelKernels/co-md_records/nmf_W_mrna.csv")

print("/\ /\ /\ /\ /\ /\    NMF_MRNA   /\ /\ /\ /\ /\ /\ ")
print("|| || || || || ||               || || || || || || ")
print("                 Input shape: ")
print("                             "+str(mrna.shape))
print(" ........................................... ")

# %% Optimize PCA

# pca = PCA(n_components=0.99)
# pca.fit(cnv.T)
# pca.components_.T.shape
# pd.DataFrame(pca.components_.T).to_csv("D:/pathFeatureLabelKernels/co-md_records/pca_W_cnv.csv")

# pca.singular_values_
pca = PCA(n_components=0.99)
pca.fit(mrna.T)
pca.components_.T.shape
pd.DataFrame(pca.components_.T).to_csv("D:/pathFeatureLabelKernels/co-md_records/pca_W_mrna.csv")
 
# %% Optimize Kernel-PCA
# from sklearn.metrics import mean_squared_error
# def my_scorer(estimator, X, y=None):
#     X_reduced = estimator.transform(X)
#     X_preimage = estimator.inverse_transform(X_reduced)
#     return -1 * mean_squared_error(X, X_preimage)

# param_grid = [{
#         "gamma": np.linspace(0.03, 1, 10), 
#         "kernel": ["rbf", "sigmoid", "linear", "poly"]
#     }]

# kpca=KernelPCA(fit_inverse_transform=True, n_jobs=-1) 
# grid_search = GridSearchCV(kpca, param_grid, cv=3, scoring=my_scorer)
# grid_search.fit(cnv)

# def objective(trial):
#     n_components = trial.suggest_int("n_components", 10, 135)
#     gamma = trial.suggest_uniform("gamma", 0.0,0.5)
#     degree = trial.suggest_int("degree", 2,3)
#     kernel = trial.suggest_categorical("kernel",  ['linear', 'poly', 'rbf', 'sigmoid'])
    
#     kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma, 
#                      degree=degree, fit_inverse_transform=True, n_jobs=-1)
    
#     kpca.fit(cnv.T)
    
#     X_reduced = kpca.transform(cnv.T)
#     X_preimage = kpca.inverse_transform(X_reduced)
    
#     error = -1 * mean_squared_error(cnv.T, X_preimage)
    
#     return error

# study = optuna.create_study(direction='minimize')
# study.optimize(objective, n_trials=500)

# %% Optimze Autoencoders
## Variational autoencoder (VAE)
### https://blog.keras.io/building-autoencoders-in-keras.html
# train autoencoder for classification with no compression in the bottleneck layer
# from sklearn.datasets import make_classification
# from tensorflow.keras.utils import plot_model
# from matplotlib import pyplot
# define dataset
# X, y = make_classification(n_samples=1000, n_features=100, n_informative=10, n_redundant=90, random_state=1)

# %%
# # number of input columns
# n_inputs = cnv.shape[1]
# # split into train test sets
# y = datos.Y
# X_train, X_test, y_train, y_test = train_test_split(cnv, y, test_size=0.33, random_state=1)
# # scale data
# t = MinMaxScaler()
# t.fit(X_train)
# X_train = t.transform(X_train)
# X_test = t.transform(X_test)

# def objective(trial):
    
#     n_bottleneck = trial.suggest_int("n_bottleneck", 100, 500)
#     enco_deco_level1 = trial.suggest_int("enco_deco_level1", 256,512)
#     enco_deco_level2 = trial.suggest_int("enco_deco_level2", 64,128)
#     num_layers = trial.suggest_int("num_layer", 0, 2)
#     activation = trial.suggest_categorical("activation", ['relu', 'selu', 'tanh', 'softmax'])
    
#     # define encoder
#     visible = Input(shape=(n_inputs,))
#     # encoder level 1
#     e = Dense(enco_deco_level1, activation=activation)(visible)
#     e = BatchNormalization()(e)
#     e = LeakyReLU()(e)
#     e = Dense(units = enco_deco_level2, activation=activation)(e)
#     e = BatchNormalization()(e)
#     e = LeakyReLU()(e)
#     for i in range(num_layers):
#         # encoder level 2
#         e = Dense(enco_deco_level2, activation=activation)(e)
#         e = BatchNormalization()(e)
#         e = LeakyReLU()(e)
#     # bottleneck
#     # n_bottleneck = round(float(n_inputs) / 2.0)
#     bottleneck = Dense(n_bottleneck, activation=activation)(e)
#     # define decoder, level 1
#     d = Dense(enco_deco_level2, activation=activation)(bottleneck)
#     d = BatchNormalization()(d)
#     d = LeakyReLU()(d)
#     for i in range(num_layers):
#     # define decoder, level 1
#         d = Dense(enco_deco_level2, activation=activation)(bottleneck)
#         d = BatchNormalization()(d)
#         d = LeakyReLU()(d)
#     # decoder level 2
#     d = Dense(enco_deco_level1, activation=activation)(d)
#     d = BatchNormalization()(d)
#     d = LeakyReLU()(d)
#     # output layer
#     output = Dense(n_inputs, activation='linear')(d)
#     # define autoencoder model
#     model = Model(inputs=visible, outputs=output)
#     # compile autoencoder model
#     model.compile(optimizer='adam', loss='mse')
#     # plot the autoencoder
#     # plot_model(model, 'autoencoder_compress.png', show_shapes=True)
#     # fit the autoencoder model to reconstruct input
#     history = model.fit(X_train, X_train, epochs=200, batch_size=16, verbose=0, validation_data=(X_test,X_test))
#     # plot loss
#     # pyplot.plot(history.history['loss'], label='train')
#     # pyplot.plot(history.history['val_loss'], label='test')
#     # pyplot.legend()
#     # pyplot.show()
#     # define an encoder model (without the decoder)
#     # encoder = Model(inputs=visible, outputs=bottleneck)
#     # plot_model(encoder, 'encoder_compress.png', show_shapes=True)
#     # save the encoder to file
#     # encoder.save('encoder.h5')
    
#     return history.history['val_loss'][-1]

# study = optuna.create_study(direction='minimize')
# study.optimize(objective, n_trials=100)

# # cnv -> GBM + LGG
# # [I 2022-01-20 09:36:16,734] Trial 8 finished with value: 0.01432463712990284 and parameters: 
#     # {'n_bottleneck': 435, 'enco_deco_level1': 446, 'enco_deco_level2': 70, 'num_layer': 1, 
#     # 'activation': 'tanh'}. Best is trial 8 with value: 0.01432463712990284.

# # cnv -> OV
# # [I 2022-01-17 11:21:15,614] Trial 6 finished with value: 0.02778436988592148 and parameters: 
#     # {'n_bottleneck': 166, 'enco_deco_level1': 417, 'enco_deco_level2': 70, 'num_layer': 0, 
#     # 'activation': 'selu'}. Best is trial 6 with value: 0.02778436988592148.
    
# # define encoder
# n_bottleneck = study.best_trial.params["n_bottleneck"]; enco_deco_level1 = study.best_trial.params["enco_deco_level1"];
# enco_deco_level2 = study.best_trial.params["enco_deco_level2"];
# num_layer = study.best_trial.params["num_layer"]; activation_inside = study.best_trial.params["activation"]
# visible = Input(shape=(n_inputs,))
# # encoder level 1
# e = Dense(n_bottleneck, activation=activation_inside)(visible)
# e = BatchNormalization()(e)
# e = LeakyReLU()(e)
# e = Dense(units = enco_deco_level1, activation=activation_inside)(e)
# e = BatchNormalization()(e)
# e = LeakyReLU()(e)
# for i in range(num_layer):
#     # encoder level 2
#     e = Dense(units = enco_deco_level1, activation=activation_inside)(e)
#     e = BatchNormalization()(e)
#     e = LeakyReLU()(e)
# # bottleneck
# # n_bottleneck = round(float(n_inputs) / 2.0)
# bottleneck = Dense(enco_deco_level2, activation=activation_inside)(e)
# # define decoder, level 1
# d = Dense(enco_deco_level1, activation=activation_inside)(bottleneck)
# d = BatchNormalization()(d)
# d = LeakyReLU()(d)
# for i in range(num_layer):
# # define decoder, level 1
#     d = Dense(units = enco_deco_level1, activation=activation_inside)(bottleneck)
#     d = BatchNormalization()(d)
#     d = LeakyReLU()(d)
# # decoder level 2
# d = Dense(n_bottleneck, activation=activation_inside)(d)
# d = BatchNormalization()(d)
# d = LeakyReLU()(d)
# # output layer
# output = Dense(n_inputs, activation='linear')(d)
# # define autoencoder model
# model = Model(inputs=visible, outputs=output)
# # compile autoencoder model
# model.compile(optimizer='adam', loss='mse')
# # plot the autoencoder
# # plot_model(model, 'autoencoder_compress.png', show_shapes=True)
# # fit the autoencoder model to reconstruct input
# history = model.fit(X_train, X_train, epochs=200, batch_size=16, verbose=0, validation_data=(X_test,X_test))
# # plot loss

# encoder = Model(inputs=visible, outputs=bottleneck)
# cnv_encode = encoder.predict(cnv)
# cnv_encode.shape
# pd.DataFrame(cnv_encode).to_csv("D:/pathFeatureLabelKernels/co-md_records/encoder_W_cnv.csv")

# print("/\ /\ /\ /\ /\ /\ ENCODER_CNV   /\ /\ /\ /\ /\ /\ ")
# print("|| || || || || ||               || || || || || || ")
# print("                 Input shape: ")
# print("                             "+str(cnv.shape))
# print(" ........................................... ")

# %% Optimze Autoencoders
# number of input columns
n_inputs = mrna.shape[1]
# split into train test sets
y = datos.Y
X_train, X_test, y_train, y_test = train_test_split(mrna, y, test_size=0.33, random_state=1)
# scale data
t = MinMaxScaler()
t.fit(X_train)
X_train = t.transform(X_train)
X_test = t.transform(X_test)

def objective(trial):
    
    n_bottleneck = trial.suggest_int("n_bottleneck", 100, 500)
    enco_deco_level1 = trial.suggest_int("enco_deco_level1", 256,512)
    enco_deco_level2 = trial.suggest_int("enco_deco_level2", 64,128)
    num_layers = trial.suggest_int("num_layer", 0, 2)
    activation = trial.suggest_categorical("activation", ['relu', 'selu', 'tanh', 'softmax'])
    
    # define encoder
    visible = Input(shape=(n_inputs,))
    # encoder level 1
    e = Dense(enco_deco_level1, activation=activation)(visible)
    e = BatchNormalization()(e)
    e = LeakyReLU()(e)
    e = Dense(units = enco_deco_level2, activation=activation)(e)
    e = BatchNormalization()(e)
    e = LeakyReLU()(e)
    for i in range(num_layers):
        # encoder level 2
        e = Dense(enco_deco_level2, activation=activation)(e)
        e = BatchNormalization()(e)
        e = LeakyReLU()(e)
    # bottleneck
    # n_bottleneck = round(float(n_inputs) / 2.0)
    bottleneck = Dense(n_bottleneck, activation=activation)(e)
    # define decoder, level 1
    d = Dense(enco_deco_level2, activation=activation)(bottleneck)
    d = BatchNormalization()(d)
    d = LeakyReLU()(d)
    for i in range(num_layers):
    # define decoder, level 1
        d = Dense(enco_deco_level2, activation=activation)(bottleneck)
        d = BatchNormalization()(d)
        d = LeakyReLU()(d)
    # decoder level 2
    d = Dense(enco_deco_level1, activation=activation)(d)
    d = BatchNormalization()(d)
    d = LeakyReLU()(d)
    # output layer
    output = Dense(n_inputs, activation='linear')(d)
    # define autoencoder model
    model = Model(inputs=visible, outputs=output)
    # compile autoencoder model
    model.compile(optimizer='adam', loss='mse')
    # plot the autoencoder
    # plot_model(model, 'autoencoder_compress.png', show_shapes=True)
    # fit the autoencoder model to reconstruct input
    history = model.fit(X_train, X_train, epochs=200, batch_size=16, verbose=0, validation_data=(X_test,X_test))
    # plot loss
    # pyplot.plot(history.history['loss'], label='train')
    # pyplot.plot(history.history['val_loss'], label='test')
    # pyplot.legend()
    # pyplot.show()
    # define an encoder model (without the decoder)
    # encoder = Model(inputs=visible, outputs=bottleneck)
    # plot_model(encoder, 'encoder_compress.png', show_shapes=True)
    # save the encoder to file
    # encoder.save('encoder.h5')
    
    return history.history['val_loss'][-1]

study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=100)


# mrna -> GBM + LGG
# [I 2022-01-20 09:26:12,392] Trial 4 finished with value: 0.020567554980516434 and parameters: 
    # {'n_bottleneck': 350, 'enco_deco_level1': 500, 'enco_deco_level2': 93, 'num_layer': 0, 
    # 'activation': 'tanh'}. Best is trial 2 with value: 0.017486678436398506.

# mrna -> OV
# [I 2022-01-17 11:28:54,427] Trial 6 finished with value: 0.06571829319000244 and parameters: 
    # {'n_bottleneck': 458, 'enco_deco_level1': 293, 'enco_deco_level2': 118, 'num_layer': 0, 
    # 'activation': 'tanh'}. Best is trial 6 with value: 0.06571829319000244.

# define encoder
n_bottleneck = study.best_trial.params["n_bottleneck"]; enco_deco_level1 = study.best_trial.params["enco_deco_level1"];
enco_deco_level2 = study.best_trial.params["enco_deco_level2"];
num_layer = study.best_trial.params["num_layer"]; activation_inside = study.best_trial.params["activation"]
visible = Input(shape=(n_inputs,))
# encoder level 1
e = Dense(n_bottleneck, activation=activation_inside)(visible)
e = BatchNormalization()(e)
e = LeakyReLU()(e)
e = Dense(units = enco_deco_level1, activation=activation_inside)(e)
e = BatchNormalization()(e)
e = LeakyReLU()(e)
for i in range(num_layer):
    # encoder level 2
    e = Dense(units = enco_deco_level1, activation=activation_inside)(e)
    e = BatchNormalization()(e)
    e = LeakyReLU()(e)
# bottleneck
# n_bottleneck = round(float(n_inputs) / 2.0)
bottleneck = Dense(enco_deco_level2, activation=activation_inside)(e)
# define decoder, level 1
d = Dense(enco_deco_level1, activation=activation_inside)(bottleneck)
d = BatchNormalization()(d)
d = LeakyReLU()(d)
for i in range(num_layer):
# define decoder, level 1
    d = Dense(units = enco_deco_level1, activation=activation_inside)(bottleneck)
    d = BatchNormalization()(d)
    d = LeakyReLU()(d)
# decoder level 2
d = Dense(n_bottleneck, activation=activation_inside)(d)
d = BatchNormalization()(d)
d = LeakyReLU()(d)
# output layer
output = Dense(n_inputs, activation='linear')(d)
# define autoencoder model 
model = Model(inputs=visible, outputs=output)
# compile autoencoder model
model.compile(optimizer='adam', loss='mse')
# plot the autoencoder
# plot_model(model, 'autoencoder_compress.png', show_shapes=True)
# fit the autoencoder model to reconstruct input
history = model.fit(X_train, X_train, epochs=200, batch_size=16, verbose=0, validation_data=(X_test,X_test))
# plot loss

encoder = Model(inputs=visible, outputs=bottleneck)
mrna_encode = encoder.predict(mrna)
mrna_encode.shape
pd.DataFrame(mrna_encode).to_csv("D:/pathFeatureLabelKernels/co-md_records/encoder_W_mrna.csv")

print("/\ /\ /\ /\ /\ /\ ENCODER_MRNA  /\ /\ /\ /\ /\ /\ ")
print("|| || || || || ||               || || || || || || ")
print("                 Input shape: ")
print("                             "+str(mrna.shape))
print(" ........................................... ")
