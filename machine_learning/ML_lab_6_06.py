#Breast Cancer Dataset Encoding & Logistic Regression
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder,LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report,accuracy_score


def breast_cancer():



    df=pd.read_csv("breast-cancer.csv",header=None)

    X=df.iloc[:,:-1]
    y=df.iloc[:,-1]

    encoder=OneHotEncoder(sparse_output=False)

    X=encoder.fit_transform(X)

    y=LabelEncoder().fit_transform(y)

    X_train,X_test,y_train,y_test=train_test_split(
        X,y,test_size=0.2,random_state=42)

    model=LogisticRegression(max_iter=1000)

    model.fit(X_train,y_train)

    pred=model.predict(X_test)

    print(classification_report(y_test,pred))

    print("Accuracy =",accuracy_score(y_test,pred))


def mymain():
    breast_cancer()


if __name__=="__main__":
    mymain()