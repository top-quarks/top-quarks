#Fit a logistic regression model for canditate triples, and output the parameters to a file
#Overwrites the previous parameters
#We fit one model for everything

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, auc, roc_curve
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from bisect import bisect
import numpy as np
import sys, os

from sklearn.externals import joblib

data = np.loadtxt('data3')

y = data[:,0]
x = data[:,1:]
poly = PolynomialFeatures(interaction_only=False, include_bias=False, degree=3)

#d1 = 1/x.std(axis=0)
#x *= d1
#print(x.mean(axis=0))
#print(x.std(axis=0))

x = poly.fit_transform(x)
print(poly.get_feature_names())
d2 = 1/x.std(axis=0)
x *= d2

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.5, random_state=43)

#x_train = x
#y_train = y
#x_test = x
#y_test = y

fp, tp = sum(y_test == 0), sum(y_test == 1)
print(fp, tp)

run = 1
if run:
    logreg = LogisticRegression(penalty="l2", C=1e2, class_weight = {0:1,1:100})#, solver="lbfgs", n_jobs=8)#, solver='saga')#penalty="l2", dual=False, tol=0.0001, C=1.0, fit_intercept=True, intercept_scaling=1, class_weight=None, random_state=None, solver=’liblinear’, max_iter=100, multi_class="ovr", verbose=0, warm_start=False, n_jobs=8)
    result = logreg.fit(x_train, y_train)
    joblib.dump(logreg, 'pair_logreg4/logreg_all3.pkl')
else:
    logreg = joblib.load('pair_logreg4/logreg_all3.pkl')

#print(poly.get_feature_names())
pred = logreg.predict_proba(x_test)[:,1]

fpr, tpr, thresholds = roc_curve(y_test, pred)
print("AUC:", auc(fpr, tpr))

ti = bisect(fpr, 1e7/2000/fp)
#ti = bisect(tpr, 0.99)
#for i in range(len(thresholds)-1, 0, -1):
#    if fpr[i]*fp < tpr[i]*tp*10 and tpr[i] != 1:
#        ti = i
#        break
thres = thresholds[ti]
print(confusion_matrix(y_test, pred>=thres))

f = open("pair_logreg4/model_all3", "w")
print(thres)
thres = -np.log(1/thres-1)
f.write(" ".join(str(i) for i in list(logreg.intercept_)+list(logreg.coef_[0,:]*d2)+[float(thres)]))
f.close()
