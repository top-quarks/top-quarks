#Fit a logistic regression model for canditate pairs, and output the parameters to a file
#Overwrites the previous parameters
#We fit one model for each of the 100 most important pairs of layers

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, auc, roc_curve
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from bisect import bisect
import numpy as np
import sys, os

from sklearn.externals import joblib

index = 0
if len(sys.argv) >= 2: index = int(sys.argv[1])

data = np.loadtxt('trained/tmp/data%d'%index)#_%d'%(ma, mb))

y = data[:,0]
x = data[:,1:]
poly = PolynomialFeatures(interaction_only=True, include_bias=False, degree=2)

#d1 = 1/x.std(axis=0)
#x *= d1
#print(x.mean(axis=0))
#print(x.std(axis=0))

x = poly.fit_transform(x)
d2 = 1/x.std(axis=0)
x *= d2

#x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.33, random_state=43)

x_train = x
y_train = y
x_test = x
y_test = y

modeli = 3
fp, tp = sum(y_test == 0), sum(y_test == 1)
print(fp, tp)
if (tp < 30):
    os.system('rm trained/pair_logreg%d/model%d'%(modeli, index))
    exit(0)
run = 1
if run:
    logreg = LogisticRegression(penalty="l2", C=1e2, class_weight = {0:1,1:100})#, solver="lbfgs", n_jobs=8)#, solver='saga')#penalty="l2", dual=False, tol=0.0001, C=1.0, fit_intercept=True, intercept_scaling=1, class_weight=None, random_state=None, solver=’liblinear’, max_iter=100, multi_class="ovr", verbose=0, warm_start=False, n_jobs=8)
    result = logreg.fit(x_train, y_train)
    joblib.dump(logreg, 'trained/pair_logreg%d/logreg%d.pkl'%(modeli, index))
else:
    logreg = joblib.load('trained/pair_logreg%d/logreg%d.pkl'%(modeli, index))

#print(poly.get_feature_names())
pred = logreg.predict_proba(x_test)[:,1]

fpr, tpr, thresholds = roc_curve(y_test, pred)
print("AUC:", auc(fpr, tpr))

ti = bisect(tpr, 0.99)
for i in range(len(thresholds)-1, 0, -1):
    if fpr[i]*fp < tpr[i]*tp*10 and tpr[i] != 1:
        ti = i
        break
thres = thresholds[ti+1]
print(confusion_matrix(y_test, pred>=thres))

f = open("trained/pair_logreg%d/model%d"%(modeli, index), "w")
print(thres)
thres = -np.log(1/thres-1)
f.write(" ".join(str(i) for i in list(logreg.intercept_)+list(logreg.coef_[0,:]*d2)+[float(thres)]))
f.close()
