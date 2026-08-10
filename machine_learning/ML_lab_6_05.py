#Logistic Regression with Manual Normalization
import numpy as np
import pandas as pd

from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score


def normalize(data):

    min_val=np.min(data,axis=0)
    max_val=np.max(data,axis=0)

    return (data-min_val)/(max_val-min_val)


def sonar():



    df=pd.read_csv("sonar.csv",header=None)

    X=df.iloc[:,:60].values
    y=df.iloc[:,60].values

    X=normalize(X)

    kf=KFold(n_splits=10,shuffle=True,random_state=42)

    scores=[]

    for train,test in kf.split(X):

        X_train,X_test=X[train],X[test]
        y_train,y_test=y[train],y[test]

        model=LogisticRegression(max_iter=1000)

        model.fit(X_train,y_train)

        pred=model.predict(X_test)

        scores.append(accuracy_score(y_test,pred))

    print("Average Accuracy =",sum(scores)/len(scores))


def mymain():
    sonar()


if __name__=="__main__":
    mymain()