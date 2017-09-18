import numpy as np
import matplotlib.pyplot as plt
import argparse
import logistic_regression as log_reg
import softmax as soft_reg
from h1_util import MlModel, load_train_data, load_test_data

        
class LogisticClassifierTwoClass(MlModel):

    def __init__(self):
        self.w = None
        self.name = 'Logistic_Regression_Two_Class'
        
    def train(self, X_, y, reg=1e-4, lr=0.1, epochs=5, batch_size=16):
        """ Train a model using mini_batch_grad_descent and assign the parameters to self.w 

        Args:
          X_: np.array shape (n,d) float - Each row is a data point
          y: np.array shape (n,)  int - Labels 
          w: np.array shape (d,)  float - Initial parameter vector
          reg: scalar - regularization parameter
          epochs: scalar - number of epochs to run
          batch_size: scalar - number of elements per mini_batch
        """

        X = np.c_[np.ones(X_.shape[0]), X_]  # Add one for bias to the first columns
        ### YOUR CODE HERE
        self.w = log_reg.mini_batch_grad_descent(X, y, self.w, reg, lr, batch_size, epochs)
        ### END CODE
        
    def visualize_model(self, ax):
        """ visualize the model by plotting the weight vector as an image """
        ax.imshow(self.w[1:].reshape(28, -1, order='F').T, cmap='bone')
        
    def predict(self, X_):
        """ Predict class for each data point in X with this model
        For all data points if the output probability is < 0.5 return class 0, else return class 1
        Args:
         X: np.array  shape (n,d)

        Returns:
         predictions: np.array shape (n,) dtype int64        
        """
        X = np.c_[np.ones(X_.shape[0]), X_]  # Add one for bias to the first columns
        predictions = np.zeros(X.shape[0])
        ### YOUR CODE HERE
        predictions = np.round(log_reg.logistic(X.dot(self.w)))
        ### END CODE
        assert predictions.shape == (X.shape[0],)
        return predictions.astype('int64')

    def probability(self, X_):
        """ Return the probability for class 1, for point x that is sigmoid(w^intercal x) 
        for each point x in the input. Needed for all-vs-one
        
        Args:
         X_: np.array  shape (n,d)
        Returns:
         probs: np.array shape (n,) 
        """
        X = np.c_[np.ones(X_.shape[0]), X_]  # Add one for bias to the first columns
        probs = np.zeros(X.shape[0])
        ### YOUR CODE HERE
        probs = log_reg.logistic(X.dot(self.w))
        ### END CODE
        assert probs.shape == (X.shape[0],)
        return probs
    
class LogisticClassifier(MlModel):
    """ Logistic Regression model for more than two classes using one-vs-all """

    def __init__(self):
        self.models = None
        self.classes = None
        self.name = 'Logistic_Regression'

        
    def visualize_model(self, ax):
        """ Plots the model weights on the given axes """
        tr = np.c_[[model.w[1:] for model in self.models]].T
        tr2 = tr.reshape(28,28, 10, order='F')
        tr3 = np.transpose(tr2, axes=[1, 0, 2])
        ax.imshow(tr3.reshape(28, -1, order='F'), cmap='bone')
        
    def train(self, X, y, reg=1e-4, lr=0.1, epochs=5, batch_size=16):
        """
        Fill self.models with a one-vs-all model one for each class y must be in [0, num_classes-1]
        Each model should be a LogisticClassifierTwoClass
        For each class c you should make a label vector that is one for y=c and 0 otherwise and train a one-vs-all classifer with for that class.       
        That way the higher the proability output by the classifier for class c the more likely that the classifier thinks the data point is in class c.
        
        Args:
          X_: np.array shape (n,d) float - Each row is a data point
          y: np.array shape (n,)  int - Labels 
          w: np.array shape (d,)  float - Initial parameter vector
          reg: scalar - regularization parameter
          epochs: scalar - number of epochs to run
          batch_size: scalar - number of elements per mini_batch
        """
        self.classes = np.sort(np.unique(y))
        self.models = [None]*10
        ### YOUR CODE HERE 5-10 lines
        for i in range(10):
            yi = np.array(y == i).astype(int)
            model_i = LogisticClassifierTwoClass()
            model_i.train(X, yi, reg, lr, epochs, batch_size)
            self.models[i] = model_i
        ### END CODE
        assert len(self.models) == len(self.classes)
        
    def predict(self, X):
        """ Predict class for each data point in X with this model
                
        np.c_ and np.argmax may be useful
        Args:
         X: numpy array, shape (n,d) each row is a data point
        
        Returns:
         predictions: numpy array shape (n,) int, prediction on each input point 
        """
        pred = np.zeros(X.shape[0])
        
        ### YOUR CODE HERE 1-3 lines
        pred = np.array([x.probability(X) for x in self.models]).argmax(axis=1)
        ### END CODE
        assert pred.shape == (X.shape[0],)
        return pred


class SoftmaxClassifier(MlModel):
    """ Softmax Model Classifier trained using mini-batch gradient descent """
    
    def __init__(self):
        self.w = None
        self.name = 'Softmax'

    def train(self, X_, y, reg=1e-4, lr=0.1, epochs=5, batch_size=16):
        """ Train a softmax model using mini_batch_grad_descent and assign the parameters to self.w 
       
        Args:
            X_: np.array shape (n,d) dtype float - Features 
            y: np.array shape (n,) dtype int - Labels 
        """
        # set up classes and make y into a matrix 1 in k encoded
        self.classes = np.sort(np.unique(y))
        num_classes = self.classes.size
        y_as_matrix = np.zeros((y.size, num_classes))
        y_as_matrix[np.arange(y.shape[0]), y] = 1
        X = np.c_[np.ones(X_.shape[0]), X_] # add bias variable 1
        ### YOUR CODE HERE        
        ### END CODE
            
    def predict(self, X_):
        """ Write the prediction algorithm for softmax classification
        np.dot (or @) and np.argmax may be useful
        
        Args:
          X_: numpy array shape (n,d) dtype float - Each row is a data point

        Returns:
          pred: numpy array shape (n,) with prediction on each input point.
        """
        X = np.c_[np.ones(X_.shape[0]), X_] # add bias variable 1
        pred = np.zeros(X.shape[0])
        ### YOUR CODE HERE
        ### END CODE
        return pred
        
    def probability(self, X_):
        """ Write the algorithm for computing the probability distribution over classes for each input point 

        Args:
          X_: np.array shape (n,d) dtype float - Features 
        """
        X = np.c_[np.ones(X_.shape[0]), X_] # add bias variable 1
        prob = np.zeros(X.shape[0], self.num_classes)
        ### YOUR CODE HERE
        ### END CODE
        return prob

    def visualize_model(self, ax):
        """ Visualize the model by plotting the weight matrix learned in training

        Args: 
          ax: matplotlib axes to draw image on
        """
        rs = self.w[1:,:].reshape(28, 28, 10, order='F')
        rs2 = np.transpose(rs, axes=[1,0,2])
        ax.imshow(rs2.reshape(28, -1, order='F'), cmap='bone')
        

def model_accuracy(model, X, y):
    """ Compute the accuracy of mode defined by w on data X, Y

    np.mean may be useful here
    Args:
    model: Object supporting predict
    X: np.array shape (n,d) dtype float - Features 
    y: np.array shape (n,) dtype int - Labels 

    Returns: 
    acc: scalar  Percentage of correct predictions (=y) for model on X 
    """
    acc = None
    ### YOUR CODE HERE 1-2 lines   
    ypredict = model.predict(X)
    acc = (y == ypredict).mean()
    ### END CODE
    return acc


    
def run_validation(model, X, y, params, val_size=0.2, **kwargs):
    """ Compute the best regularization parameter using validation and return a model trained with that parameter set on the full data

    Args:
        X: np.array shape (n,d) dtype float - Features 
        y: np.array shape (n,) dtype int - Labels 
        model: An MlModel (support train, predict)
        params: list of regularization parameters to try, usually something like [0.3**i for i in range(10)]
        val_size: scalar, fraction to use for validation
    Returns:
       model: input MlModel trained with best regularization params on all training data (X,y)
       acc: numpy array, the validation scores for params input
    """
    acc =  np.zeros(len(params))
    n = y.size
    val_size = int(n*0.2)
    val_train = X[0:val_size,:]
    val_target = y[0:val_size]
    train = X[val_size+1:,:]
    target = y[val_size+1:]
    acc = np.zeros(len(params))
    ### YOUR CODE HERE 5-10 lines
    ### END CODE
    return model, acc
        
