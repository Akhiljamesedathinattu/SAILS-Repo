from sklearn.linear_model import LinearRegression
from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score





def load_data():
    [X,y]=fetch_california_housing(return_X_y=True)

    return(X,y)






















def mymain():
    [X,y]=load_data()
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30, random_state=999)
    print('----------TRAINING----------')
    print("N=%d"% (len(X)))
    #training  LINEAR REGRESSION
    model=LinearRegression()
    #train the model
    model.fit(X_train, y_train)
    #prediction on a test set
    y_pred = model.predict(X_test)
    #perfomance measure   r2 score  this will give us a number if it close to 1 the model is good
    r2=r2_score(y_test, y_pred)
    print("r2 score is %0.2f (closer to 1 is good)"% r2)
    print("done")





if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()