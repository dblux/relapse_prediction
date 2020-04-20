#%% IMPORT PACKAGES

import glob, re
import numpy as np, pandas as pd
from scipy.stats import norm

import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

from sklearn import trees
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC

from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc

#%% IMPORT DATA
rpath_list = sorted(glob.glob("temp/pca_bcm_qpsp/*.tsv"))
print(len(rpath_list))

# Subtype labels
LABEL_RPATH = "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df = pd.read_csv(LABEL_RPATH, sep="\t")
MRD_RPATH = "data/GSE67684/processed/metadata/metadata-label_mrd_subtype.tsv"
mrd_df = pd.read_csv(MRD_RPATH, sep="\t")

#%% IMPORT DATA
features_rpath = rpath_list[5]
raw_data = pd.read_csv(features_rpath, sep="\t")
SUBTYPE = re.search("pca_bcm_qpsp-(.+?).tsv", features_rpath).groups()[0]
print(SUBTYPE)

# Feature selection
features = ["erm1", "l2norm_d0_d8", "d0_normal_proj", "angle_d0_d8"]
X = raw_data.loc[:,features]
y = metadata_df.loc[X.index, "label"]
y_label = y.astype("category")
# y_label = y_label.cat.rename_categories(["remission","relapse"])
X_y = pd.concat([X,y_label], axis=1)
# metadata = metadata_df.loc[X.index,:]
# metadata.groupby(["batch_info", "label"]).size()

SPLIT = 0.3
RAND_STATE = 0
X_train, _, _, _ = train_test_split(X, y,
                                    random_state=RAND_STATE,
                                    test_size=SPLIT)

# Scaling of data
## Scale using positive samples in training dataset in order
## to avoid affected by the proportion of positive and negative samples
scaler = StandardScaler()
X_train_remission = X_train.loc[y==0,]
scaler.fit(X_train_remission)
X_scaled = pd.DataFrame(scaler.transform(X),
                        index=X.index, columns=list(X))

X_scaled_train, X_scaled_test, y_train, y_test = train_test_split(
    X_scaled, y, random_state=RAND_STATE, test_size=SPLIT)

print("Size of training set:", y_train.shape)
print("Size of test set:", y_test.shape)

#%% PCA Visualisaton of features

# pca = PCA(n_components=3)
# X_pca = pca.fit_transform(X_scaled)
# col = np.where(y == 0, "steelblue", "red")
# X_pca_remission = X_pca[y==0,:]
# X_pca_relapse = X_pca[y==1,:]

# fig1, ax1 = plt.subplots(1,1)
# ax1.scatter(X_pca[:,0], X_pca[:,1], c=col)
# plt.show()

# # PCA: 3D
# fig2, ax2 = plt.subplots(1, 1)
# ax2 = Axes3D(fig2, elev=-150, azim=110)
# ax2.scatter(X_pca_remission[:,0], X_pca_remission[:,1],
#             X_pca_remission[:,2], c="blue")
# ax2.scatter(X_pca_relapse[:,0], X_pca_relapse[:,1],
#             X_pca_relapse[:,2], c="red")
# ax2.set_xlabel("PC1")
# ax2.set_ylabel("PC2")
# ax2.set_zlabel("PC3")
# plt.show()

#%% STRIP PLOT: Unnormalised

# fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(9,7))
# sns.stripplot(data=X_y, x="erm1", y="label", ax=ax1)
# sns.stripplot(data=X_y, x="l2norm_d0_d8", y="label", ax=ax2)
# sns.stripplot(data=X_y, x="d0_normal_proj", y="label", ax=ax3)
# sns.stripplot(data=X_y, x="angle_d0_d8", y="label", ax=ax4)
# plt.tight_layout()
# plt.show()

# RAW_WPATH = "dump/raw-{}.pdf".format(SUBTYPE)
# fig.savefig(RAW_WPATH)

#%% STRIP PLOT: Scaled
# X_stack = X_scaled.stack().reset_index().rename(
#     columns={"level_0": "pid",
#              "level_1": "features",
#              0: "values"}
# ).set_index("pid")

# X_stack_y = pd.concat([X_stack, y[X_stack.index]], axis=1)
# print(X_stack_y)

# plt.figure(figsize=(8,4))
# ax = sns.stripplot(x="values", y="features", data=X_stack_y, hue="label")
# ax.legend_.remove()
# plt.tight_layout()

# plt.savefig("dump/scaled-{}.pdf".format(SUBTYPE))

#%% Likelihood calculation
## Estimate parameters using the remission samples from training set
X_remission_train = X_scaled_train.loc[y_train == 0,:]
mu_remission_train = X_remission_train.mean(axis=0)
sd_remission_train = X_remission_train.std(axis=0)

#%% Calculate probability of each feature

likelihood = X_scaled.apply(norm.pdf, axis=1,
                            loc=mu_remission_train,
                            scale=sd_remission_train)

likelihood_df = pd.DataFrame(likelihood.values.tolist(),
                             index=X.index, columns=list(X))

likelihood_y = pd.concat([likelihood_df,y], axis=1)
print(likelihood_y)

## Assuming conditional independence
product_likelihood = likelihood_df.mean(axis=1)
p = 1-product_likelihood

product_likelihood_y = pd.concat(
    [product_likelihood, p, y], axis=1
).rename(columns={0: "joint", 1: "p"})
product_likelihood_y.label = product_likelihood_y.label.astype("category")
print(product_likelihood_y)

## Evaluation
likelihood_train, likelihood_test = train_test_split(product_likelihood_y,
                                                     test_size=SPLIT,
                                                     random_state=RAND_STATE)

# X_scaled_y = pd.concat([X_scaled,y], axis=1)
# print(X_scaled_y)
s
#%% STRIP PLOT: 1-joint

fig, (ax1, ax2) = plt.subplots(2, 1)
sns.stripplot(data=likelihood_train, x="p", y="label", ax=ax1)
sns.stripplot(data=likelihood_test, x="p", y="label", ax=ax2)
ax1.set_title("Training set")
ax2.set_title("Test set")
plt.tight_layout()
plt.show()

PROBA_WPATH = "dump/predict_proba-{}.pdf".format(SUBTYPE)
# fig.savefig(PROBA_WPATH)

#%% STRIP PLOT: Feature likelihoods

likelihood_stack = likelihood_df.stack().reset_index().rename(
    columns={"level_0": "pid",
             "level_1": "likelihood",
             0: "values"}
    ).set_index("pid")

likelihood_stack_y = pd.concat([likelihood_stack,
                                y[likelihood_stack.index]], axis=1)
print(likelihood_stack_y)

plt.figure(figsize=(8,4))
ax = sns.stripplot(data=likelihood_stack_y,
                   x="values", y="likelihood", hue="label")
ax.legend_.remove()
plt.tight_layout()

LIKELIHOOD_WPATH = "dump/feat_likelihood-{}.pdf".format(SUBTYPE)
# plt.savefig(LIKELIHOOD_WPATH)

#%% ROC curve from proba
fpr, tpr, _ = roc_curve(likelihood_test.iloc[:,2],
                        likelihood_test.iloc[:,1],
                        drop_intermediate=False)
AUC = auc(fpr, tpr)

print(product_likelihood_y.sort_values(by=["joint"]))

fig, ax = plt.subplots(1, 1, figsize=(6,6))
lab = "ROC (AUC: {:.3f})".format(AUC)
ax.plot(fpr, tpr, label=lab)
ax.set_xlabel("FPR")
ax.set_ylabel("TPR")
ax.legend(loc="lower right")
plt.tight_layout()
plt.show()

ROC_WPATH = "dump/roc-{}.pdf".format(SUBTYPE)s
fig.savefig(ROC_WPATH)

#%% PCA Visualisaton of likelihood features
pca = PCA(n_components=3)
X_pca = pca.fit_transform(likelihood_df)
col = np.where(y == 0, "steelblue", "red")
X_pca_remission = X_pca[y==0,:]
X_pca_relapse = X_pca[y==1,:]

fig1, ax1 = plt.subplots(1,1)
ax1.scatter(X_pca[:,0], X_pca[:,1], c=col)
plt.show()

# PCA: 3D
fig2, ax2 = plt.subplots(1, 1)
ax2 = Axes3D(fig2, elev=-170, azim=80)
ax2.scatter(X_pca_remission[:,0], X_pca_remission[:,1],
            X_pca_remission[:,2], c="blue")
ax2.scatter(X_pca_relapse[:,0], X_pca_relapse[:,1],
            X_pca_relapse[:,2], c="red")
ax2.set_xlabel("PC1")
ax2.set_ylabel("PC2")
ax2.set_zlabel("PC3")
plt.show()

#%% LIKELIHOOD RATIO

gnb_clf = GaussianNB()
# Fit to X1
gnb_clf.fit(X_train, y_train)
# Fitted parameters for each class
mean_arr = gnb_clf.theta_
var_arr = gnb_clf.sigma_
sd_arr = np.sqrt(var_arr)

# Calculate for each feature
def calcLR(col):
    col_idx = features.index(col.name)
    ## P(x|y=relapse)
    likelihood_1 = norm.pdf(col, mean_arr[1,col_idx],
                            sd_arr[1,col_idx])
    ## P(x|y=remission)
    likelihood_0 = norm.pdf(col, mean_arr[0,col_idx], # Remission
                            sd_arr[0,col_idx])
    ## For calculation stability (replace 0 with 1e-12)
    likelihood_0[likelihood_0 == 0] = 1e-12
    likelihood_1[likelihood_1 == 0] = 1e-12

    return likelihood_1/likelihood_0

lr = X.apply(calcLR, axis=0)

lr_train, lr_test, y_train, y_test = train_test_split(lr, y,
                                                      random_state=0,
                                                      test_size=SPLIT)

#%% LOGISTIC REGRESSION
logreg = LogisticRegression(random_state=0, solver="lbfgs")
logreg.fit(lr_train, y_train)
yhat_train = logreg.predict(lr_train)
yhat_test = logreg.predict(lr_test)
proba_train = logreg.predict_proba(lr_train)
proba_test = logreg.predict_proba(lr_test)

print(confusion_matrix(y_train, yhat_train))
print(confusion_matrix(y_test, yhat_test))

print(logreg.score(lr_train, y_train))
print(logreg.score(lr_test, y_test))

#%% Plot prediction probabilities
mrd = np.log10(mrd_df.d33_mrd)

## Training set
pid_train = [pid[:4] for pid in y_train.index]
mrd_train = mrd.loc[pid_train]
col_train = np.where(y_train == 0, "steelblue", "red")

## Test set
pid_test = [pid[:4] for pid in y_test.index]
mrd_test = mrd.loc[pid_test]
col_test = np.where(y_test == 0, "steelblue", "red")

## Plot
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,4))
fig.suptitle("{}: Likelihood ratios".format(SUBTYPE))

ax1.scatter(mrd_train, proba_train[:,1], c=col_train)
ax1.set_xlabel("log10(D33_MRD)")
ax1.set_ylabel("p")
ax1.set_title("Training set")

ax2.scatter(mrd_test, proba_test[:,1], c=col_test)
ax2.set_xlabel("log10(D33_MRD)")
ax2.set_ylabel("p")
ax2.set_title("Test set")

plt.tight_layout()
# fig.savefig("dump/{}_lr-proba.pdf".format(SUBTYPE))

#%% SVM
svm = LinearSVC(random_state=0, tol=1e-5)
svm.fit(lr_train, y_train)
yhat_train = svm.predict(lr_train)
yhat_test = svm.predict(lr_test)
proba_train = svm.predict_proba(lr_train)
proba_test = svm.predict_proba(lr_test)

print(confusion_matrix(y_train, yhat_train))
print(confusion_matrix(y_test, yhat_test))

print(svm.score(lr_train, y_train))
print(svm.score(lr_test, y_test))

#%%
## CHANGE MRD33 values from 0 to NA

### Likelihood ratios from GNB ###
## OPTION 1
subtypes = ["TEL-AML1", "T-ALL", "BCR-ABL",
            "E2A-PBX1", "Hyperdiploid", "Others"]
SUBTYPE = subtypes[0]
# OPTION 2
# SUBTYPE = "qpsp"

# OPTION 1
RPATH1 = "dump/{}-lr.tsv".format(SUBTYPE)
RPATH2 = "dump/qpsp_bcmD0_PCA/features-{}.tsv".format(SUBTYPE)
# OPTION 2
# RPATH1 = "dump/global_{}-lr.tsv".format(SUBTYPE)
# RPATH2 = "dump/features-{}.tsv".format(SUBTYPE)
# ## OPTION 3
# RPATH1 = "dump/{}-lr.tsv".format(SUBTYPE)
# RPATH2 = "dump/features2-{}.tsv".format(SUBTYPE)
print(RPATH1, RPATH2)

lr = pd.read_csv(RPATH1, sep="\t")
X = lr.iloc[:,0:4]
y = lr.iloc[:,4]

print(X.shape[0])
print("Features: ", list(X))

# Create new instance of logistic regression classifier
logreg = LogisticRegression(random_state=0, solver="lbfgs")
logreg.fit(X,y)
logreg.score(X,y)
proba = logreg.predict_proba(X)

# Attributes of the classifier
beta = logreg.coef_
beta_0 = logreg.intercept_

print()
with np.printoptions(precision=3, suppress=True):
    beta_all = np.concatenate([beta_0, np.squeeze(beta)])
    print("Coefficients =", beta_all)
print()

# # np.sum sums ALONG the axis
# # Before logistic function
# X_beta = np.sum(X*beta, axis = 1) + beta_0
# # Prediction tallies with coefficients
# pred_class = 1/(1+np.exp(-X_beta))

### Original features ###
# Import features
features = pd.read_csv(RPATH2, sep="\t")
### Set MRD33 values of 0 to NaN
features.loc[features["d33_mrd"] == 0, "d33_mrd"] = float("nan")
# OPTION 2: No MRD33
X_selected = features.iloc[:,[0,6,8,10]]

logreg1 = LogisticRegression(random_state=0, solver="lbfgs")
logreg1.fit(X_selected,y)
logreg1.score(X_selected,y)
proba1 = logreg1.predict_proba(X_selected)

# Attributes of the classifier
beta1 = logreg1.coef_
beta1_0 = logreg1.intercept_

# # np.sum sums ALONG the axis
# # Before logistic function
# X_beta1 = np.sum(X_selected*beta1, axis = 1) + beta1_0
# # Prediction tallies with coefficients
# pred_class1 = 1/(1+np.exp(-X_beta1))

# EDA
mrd = np.log10(features.d33_mrd)
lab_col = np.where(y == 0, "steelblue", "red")
fig1, ax1 = plt.subplots(1,1)
ax1.scatter(mrd, proba[:,1], c=lab_col)
ax1.set_xlabel("log10(D33_MRD)")
ax1.set_ylabel("p")
ax1.set_title("{}: LogReg - Likelihood ratios (Original)".format(SUBTYPE))
# fig1.savefig("dump/{}_lr-proba_original.pdf".format(SUBTYPE))

fig2, ax2 = plt.subplots(1,1)
ax2.scatter(mrd, proba1[:,1], c=lab_col)
ax2.set_xlabel("log10(D33_MRD)")
ax2.set_ylabel("p")
ax2.set_title("{}: LogReg - Features (Original)".format(SUBTYPE))
# fig2.savefig("dump/{}_features-proba_original.pdf".format(SUBTYPE))

# # Logits
# plt.scatter(mrd, X_beta, c = y+1)
# plt.scatter(mrd, X_beta1, c = y+1)

#%%
# Plot features
X_y = pd.concat([X,y], axis=1)
fig, ax = plt.subplots(1,1)
sns.stripplot(x="d0_normal_proj", y="label", data=X_y, ax=ax)
plt.show()

#%%

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(5,8))
fig.subplots_adjust(hspace=0.5)
fig.suptitle("{} ({}; {})".format(SUBTYPE, score1, score2))

sns.stripplot(x="erm1", y="label", data=data1, ax=ax1)
sns.stripplot(x="l2norm_d0_d8", y="label", data=data1, ax=ax2)
sns.stripplot(x="angle_d0_d8", y="label", data=data1, ax=ax3)
sns.stripplot(x="d0_normal_proj", y="label", data=data1, ax=ax4)

#ax1.axvline(class_mean1[0,0]); ax1.axvline(class_mean1[1,0], color="orange")
#ax2.axvline(class_mean1[0,1]); ax2.axvline(class_mean1[1,1], color="orange")
#ax3.axvline(class_mean1[0,2]); ax3.axvline(class_mean1[1,2], color="orange")
#ax4.axvline(class_mean2[0,2]); ax4.axvline(class_mean2[1,2], color="orange")

ax1.axvline(0.8, color="r")
ax2.axvline(1, color="r"); ax2.axvline(1.96, color="g")
ax3.axvline(49, color="g")
ax4.axvline(2.6, color="r")


fig.savefig("dump/strplt_threshold-{}.pdf".format(SUBTYPE))


#%%
### INCLUDE DUMMY VARIABLES FOR SUBTYPE IN LOGISTIC REGRESSION
## CHANGE MRD33 values from 0 to NA

### Likelihood ratios from GNB ###
# subtypes = ["TEL-AML1", "T-ALL", "BCR-ABL",
#             "E2A-PBX1", "Hyperdiploid", "Others"]
# SUBTYPE = subtypes[5]
### OPTION 2
SUBTYPE = "sampled_Others"

# ARG1 = 2

### OPTION 1
# RPATH1 = "dump/{}-lr.tsv".format(SUBTYPE)
# RPATH2 = "dump/qpsp_bcmD0_PCA/features-{}.tsv".format(SUBTYPE)
### OPTION 2
# RPATH1 = "dump/global_{}-lr.tsv".format(SUBTYPE)
# RPATH2 = "dump/features-{}.tsv".format(SUBTYPE)
### OPTION 3
RPATH1 = "dump/{}-lr.tsv".format(SUBTYPE)
RPATH2 = "dump/features2-{}.tsv".format(SUBTYPE)
RPATH3 = "data/GSE67684/processed/metadata-label_mrd_subtype.tsv"
print(RPATH1, RPATH2)

lr = pd.read_csv(RPATH1, sep="\t")
features = pd.read_csv(RPATH2, sep="\t")
metadata = pd.read_csv(RPATH3, sep="\t")
# Only select patients that are present in data
subtype_info = metadata.subtype.to_frame().loc[lr.index,:]

onehot_encoder = OneHotEncoder(drop="first", sparse=False)
onehot_encoder.fit(subtype_info)
enc_categories = onehot_encoder.categories_[0] # List of diff factors
dummy_X = onehot_encoder.transform(subtype_info)
# Check that m-1 dummy variables are created
# dummy_X.shape[1]
# len(enc_categories)

dummy_df = pd.DataFrame(dummy_X, index=lr.index, columns=enc_categories[1:])
X = pd.concat([lr.iloc[:,0:4], dummy_df], axis=1)
y = lr.y

print(X.shape[0])
print("Features: ", list(X))

# Create new instance of logistic regression classifier
logreg = LogisticRegression(random_state=0, solver="lbfgs", max_iter=1000)
logreg.fit(X,y)
logreg.score(X,y)

logreg.predict(X)
proba = logreg.predict_proba(X)

# Attributes of the classifier
beta = logreg.coef_
beta_0 = logreg.intercept_

print()
with np.printoptions(precision=3, suppress=True):
    beta_all = np.concatenate([beta_0, np.squeeze(beta)])
    print("Coefficients =", beta_all)
print()

# # np.sum sums ALONG the axis
# # Before logistic function
# X_beta = np.sum(X*beta, axis = 1) + beta_0
# # Prediction tallies with coefficients
# pred_class = 1/(1+np.exp(-X_beta))

### Original features ###
# Import features
### Set MRD33 values of 0 to NaN
features.loc[features.d33_mrd == 0, "d33_mrd"] = float("nan")

### Classifier using features
X1 = features.iloc[:,[0,6,8,10]]
X_selected = pd.concat([X_selected, dummy_df], axis=1)


logreg1 = LogisticRegression(random_state=0, solver="lbfgs", max_iter=1000)
logreg1.fit(X_selected,y)
logreg1.score(X_selected,y)
proba1 = logreg1.predict_proba(X_selected)

# Attributes of the classifier
beta1 = logreg1.coef_
beta1_0 = logreg1.intercept_

# # np.sum sums ALONG the axis
# # Before logistic function
# X_beta1 = np.sum(X_selected*beta1, axis = 1) + beta1_0
# # Prediction tallies with coefficients
# pred_class1 = 1/(1+np.exp(-X_beta1))

# EDA
mrd = np.log10(features.d33_mrd)
lab_col = np.where(y == 0, "steelblue", "red")
fig1, ax1 = plt.subplots(1,1)
ax1.scatter(mrd, proba[:,1], c=lab_col)
ax1.set_xlabel("log10(D33_MRD)")
ax1.set_ylabel("p")
ax1.set_title("{}: LogReg - Likelihood ratios".format(SUBTYPE))
fig1.savefig("dump/{}_lr_modmat-proba.pdf".format(SUBTYPE))

fig2, ax2 = plt.subplots(1,1)
ax2.scatter(mrd, proba1[:,1], c=lab_col)
ax2.set_xlabel("log10(D33_MRD)")
ax2.set_ylabel("p")
ax2.set_title("{}: LogReg - Features".format(SUBTYPE))
fig2.savefig("dump/{}_features_modmat-proba.pdf".format(SUBTYPE))

#%%
##### MRD33 classification results #####
### Global GSS ###
global_features = pd.read_csv("dump/features-cs_quantile.tsv", sep="\t")
### Set MRD33 values of 0 to NaN
global_features.loc[global_features.d33_mrd == 0, "d33_mrd"] = float("nan")
y_all = global_features["label"]
n = global_features.shape[0]

mrd33 = global_features.iloc[:,14:16]

# LOW RISK
low_risk = mrd33[mrd33.iloc[:,0] <= 0.0001]
false_low = low_risk.index[low_risk.label == 1]
print(false_low.shape, low_risk.shape[0])
# Relapse classfied as low: 9/100

# HIGH RISK
high_risk = mrd33[mrd33.iloc[:,0] >= 0.01]
false_high = high_risk.index[high_risk.label == 0]
print(false_high.shape, high_risk.shape[0])
# Remission classified as high: 8/25

# MEDIUM RISK
medium_risk = mrd33[(0.0001 < mrd33.d33_mrd) & (mrd33.d33_mrd < 0.01)]
plt.scatter(np.arange(medium_risk.shape[0]), medium_risk.d33_mrd,
            c = medium_risk.label)
plt.ylim(0,0.01)

# Examine misclassifications
csquantile_features.loc[false_low,["erm1","l2norm_d0_d8",
                                   "d0_normal_proj","angle_d0_d8"]]
csquantile_features.loc[false_high,["erm1","l2norm_d0_d8",
                                    "d0_normal_proj","angle_d0_d8"]]

# Exploring the feasibility of building a likelihood ratio in GLOBAL GSS

plt.hist(csquantile_features.d0_normal_proj[
    csquantile_features.label == 0],bins=30)
plt.hist(csquantile_features.d0_normal_proj[
    csquantile_features.label == 1], bins=30)

csquantile_features = pd.read_csv("dump/features-cs_quantile.tsv", sep="\t")

# In the global GSS can ERM1 be used to salvage MRD33?
# Is there any correlation between ERM1 and MRD33?
fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.scatter(np.log10(csquantile_features.d33_mrd),
            csquantile_features.erm1,
            c=csquantile_features.label)
ax.set_xlabel("log10(D33_MRD)")
ax.set_ylabel("ERM1")
ax.set_title("CS-Quantile - Global GSS")
fig.savefig("dump/csquantile-erm1_logmrd.pdf")

# No correlation. Also there is no independent threshold

qpsp_features = pd.read_csv("dump/features-qpsp.tsv", sep="\t")

# In the global GSS can ERM1 be used to salvage MRD33?
# Is there any correlation between ERM1 and MRD33?
fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.scatter(np.log10(qpsp_features.d33_mrd),
            qpsp_features.erm1,
            c=qpsp_features.label)
ax.set_xlabel("log10(D33_MRD)")
ax.set_ylabel("ERM1")
ax.set_title("QPSP - Global GSS")
fig.savefig("dump/qpsp-erm1_logmrd.pdf")

### Plot pairwise in R for global GSS
### There is some merit in trying out whether SVM can separate global GSS

#%% DECISION TREE - IMPORT DATA

# FEATURES_RPATH = "dump/features-qpsp.tsv"
FEATURES_RPATH = "dump/features-quantile.tsv"
# FEATURES_RPATH = "dump/features-cs_quantile.tsv"

data = pd.read_csv(FEATURES_RPATH, sep="\t")
# Feature selection
X = data.iloc[:,np.r_[0,8,10]]
# OPTION 2
# X = data.iloc[:,np.r_[0,8:11]]
list(X)

#%% DECISION TREE

# for i in range(len(rpath_list)):
#     i = 0
#     print(rpath_list[i])
#     data = pd.read_csv(rpath_list[i], sep="\t")

data = pd.read_csv(FEATURES_RPATH, sep="\t")

# Feature selection
X = data.iloc[:,np.r_[0,8:11]]
print(list(X))
print("X.shape =", X.shape)
# Create labels
y = data.iloc[:,15]

# Create and train model
clf = tree.DecisionTreeClassifier(random_state=0, max_depth=3)
clf = clf.fit(X, y)
# Results of training set
print(clf.score(X,y))

# Decision rules
tree_rules = tree.export_text(clf, feature_names=list(X))
print(tree_rules)

#%%

subtype = subtype_labels.loc[X.index.values, "subtype"]
# Concatenate subtype and truth labels
X1 = pd.concat([X, subtype, y],1).sort_values(["erm1","l2norm_d0_d8"])
list_dfs = list(X1.groupby("subtype"))
X5 = list_dfs[7][1]
print(list_dfs[7][0])
X2 = X5[X5["erm1"] < 74]
# X3 = X2[X2["l2norm_d0_d8"] > 0.54]
# X4 = X3[X3["angle_d0_d8"] > 49]
# X2.sort_values("l2norm_d0_d8")
# X1[X1["label"] == 1]

X2.sort_values("l2norm_d0_d8")


# In[295]:


# Construct decision tree
X2 = X1[X1["erm1"] < 0.8]
X3 = X2[X2["l2norm_d0_d8"] < 1.96]
X4 = X3[X3["angle_d0_d8"] > 49]
X4

# Random stratified split
# X_train, X_test, y_train, y_test = train_test_split(X, y,
#                                                     random_state=1,
#                                                     test_size=0.2)


# In[292]:

get_ipython().run_line_magic('matplotlib', 'notebook')
# Visualise rules of decision tree
decision_plot = tree.plot_tree(clf, feature_names=list(X),
               class_names=("remission", "relapse"), filled=True)

# In[467]:

i = 4
print(rpath_list[i])
data = pd.read_csv(rpath_list[i], sep="\t")

# Feature selection
X = data.iloc[:,np.r_[0,8:11]]
print("X.shape =", X.shape)
print(list(X))
# Create labels
y = data.iloc[:,15]

# X1 = pd.concat([X,y],1).sort_values(["erm1","l2norm_d0_d8"])
X1 = pd.concat([X,y],1).sort_values(["erm1"])
X1


# In[473]:
# Construct decision tree
X2 = X1[X1["erm1"] < 0.8]
X3 = X2[X2["l2norm_d0_d8"] < 1.96]
X4 = X3[X3["angle_d0_d8"] > 49]

# In[430]:

y_predict = data["d33_mrd"] > 0.0001
y_predict = y_predict.astype(int)

sum(y_predict != y)
pd.concat([y_predict, data.iloc[:,14:16]], 1)
