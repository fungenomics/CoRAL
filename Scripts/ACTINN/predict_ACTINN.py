import numpy as np
import pandas as pd
import sys
import math
import collections
import tensorflow as tf
tf.compat.v1.disable_eager_execution()
import argparse
import timeit
run_time = timeit.default_timer()
from tensorflow.python.framework import ops
import pickle
import os
import random

#Set seed
random.seed(123456) 

def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label



def get_parser(parser=None):
    if parser == None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-ts", "--test_set", type=str, help="Training set file path.")
    parser.add_argument("-mp", "--model_path", type=str, help="Model file path.")
    parser.add_argument("-pp", "--pred_path", type=str, help="Output prediction file path.")
    parser.add_argument("-thr", "--threshold", type=float, help="Threshold for rejection file path.")
    return parser


# Get common genes, normalize  and scale the sets
def scale_sets(sets):
    # input -- set to be scaled
    # output -- scaled set
    sets = sets.to_numpy()
    sets = np.divide(sets, np.sum(sets, axis=0, keepdims=True)) * 10000
    sets = np.log2(sets+1)
    return sets



# Turn labels into matrix
def one_hot_matrix(labels, C):
    # input -- labels (true labels of the sets), C (# types)
    # output -- one hot matrix with shape (# types, # samples)
    C = tf.constant(C, name = "C")
    one_hot_matrix = tf.one_hot(labels, C, axis = 0)
    sess = tf.compat.v1.Session()
    one_hot = sess.run(one_hot_matrix)
    sess.close()
    return one_hot



# Make types to labels dictionary
def type_to_label_dict(types):
    # input -- types
    # output -- type_to_label dictionary
    type_to_label_dict = {}
    all_type = list(set(types))
    for i in range(len(all_type)):
        type_to_label_dict[all_type[i]] = i
    return type_to_label_dict



# Convert types to labels
def convert_type_to_label(types, type_to_label_dict):
    # input -- list of types, and type_to_label dictionary
    # output -- list of labels
    types = list(types)
    labels = list()
    for type in types:
        labels.append(type_to_label_dict[type])
    return labels



# Function to create placeholders
def create_placeholders(n_x, n_y):
    X = tf.compat.v1.placeholder(tf.float32, shape = (n_x, None))
    Y = tf.compat.v1.placeholder(tf.float32, shape = (n_y, None))
    return X, Y



# Initialize parameters
def initialize_parameters(nf, ln1, ln2, ln3, nt):
    # input -- nf (# of features), ln1 (# nodes in layer1), ln2 (# nodes in layer2), nt (# types)
    # output -- a dictionary of tensors containing W1, b1, W2, b2, W3, b3
    tf.compat.v1.set_random_seed(3) # set seed to make the results consistant
    W1 = tf.compat.v1.get_variable("W1", [ln1, nf], initializer = tf.compat.v1.keras.initializers.VarianceScaling(scale=1.0, mode="fan_avg", distribution="uniform", seed = 3))
    b1 = tf.compat.v1.get_variable("b1", [ln1, 1], initializer = tf.compat.v1.zeros_initializer())
    W2 = tf.compat.v1.get_variable("W2", [ln2, ln1], initializer = tf.compat.v1.keras.initializers.VarianceScaling(scale=1.0, mode="fan_avg", distribution="uniform", seed = 3))
    b2 = tf.compat.v1.get_variable("b2", [ln2, 1], initializer = tf.compat.v1.zeros_initializer())
    W3 = tf.compat.v1.get_variable("W3", [ln3, ln2], initializer = tf.compat.v1.keras.initializers.VarianceScaling(scale=1.0, mode="fan_avg", distribution="uniform", seed = 3))
    b3 = tf.compat.v1.get_variable("b3", [ln3, 1], initializer = tf.compat.v1.zeros_initializer())
    W4 = tf.compat.v1.get_variable("W4", [nt, ln3], initializer = tf.compat.v1.keras.initializers.VarianceScaling(scale=1.0, mode="fan_avg", distribution="uniform", seed = 3))
    b4 = tf.compat.v1.get_variable("b4", [nt, 1], initializer = tf.compat.v1.zeros_initializer())
    parameters = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
    return parameters



# Forward propagation function
def forward_propagation(X, parameters):
    # function model: LINEAR -> RELU -> LINEAR -> RELU -> LINEAR -> SOFTMAX
    # input -- dataset with shape (# features, # sample), parameters "W1", "b1", "W2", "b2", "W3", "b3"
    # output -- the output of the last linear unit
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    # forward calculations
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4



# Compute cost
def compute_cost(Z4, Y, parameters, lambd=0.01):
    # input -- Z3 (output of forward propagation with shape (# types, # samples)), Y (true labels, same shape as Z3)
    # output -- tensor of teh cost function
    logits = tf.transpose(a=Z4)
    labels = tf.transpose(a=Y)
    cost = tf.reduce_mean(input_tensor=tf.nn.softmax_cross_entropy_with_logits(logits = logits, labels = labels)) + \
    (tf.nn.l2_loss(parameters["W1"]) + tf.nn.l2_loss(parameters["W2"]) + tf.nn.l2_loss(parameters["W3"]) + tf.nn.l2_loss(parameters["W4"])) * lambd
    return cost



# Get the mini batches
def random_mini_batches(X, Y, mini_batch_size=32, seed=1):
    # input -- X (training set), Y (true labels)
    # output -- mini batches
    ns = X.shape[1]
    mini_batches = []
    np.random.seed(seed)
    # shuffle (X, Y)
    permutation = list(np.random.permutation(ns))
    shuffled_X = X[:, permutation]
    shuffled_Y = Y[:, permutation]
    # partition (shuffled_X, shuffled_Y), minus the end case.
    num_complete_minibatches = int(math.floor(ns/mini_batch_size)) # number of mini batches of size mini_batch_size in your partitionning
    for k in range(0, num_complete_minibatches):
        mini_batch_X = shuffled_X[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch_Y = shuffled_Y[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    # handling the end case (last mini-batch < mini_batch_size)
    if ns % mini_batch_size != 0:
        mini_batch_X = shuffled_X[:, num_complete_minibatches * mini_batch_size : ns]
        mini_batch_Y = shuffled_Y[:, num_complete_minibatches * mini_batch_size : ns]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    return mini_batches



# Forward propagation for prediction
def forward_propagation_for_predict(X, parameters):
    # input -- X (dataset used to make prediction), papameters after training
    # output -- the output of the last linear unit
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4



# Predict function
# def predict(X, parameters):
#     # input -- X (dataset used to make prediction), papameters after training
#     # output -- prediction
#     W1 = tf.convert_to_tensor(value=parameters["W1"])
#     b1 = tf.convert_to_tensor(value=parameters["b1"])
#     W2 = tf.convert_to_tensor(value=parameters["W2"])
#     b2 = tf.convert_to_tensor(value=parameters["b2"])
#     W3 = tf.convert_to_tensor(value=parameters["W3"])
#     b3 = tf.convert_to_tensor(value=parameters["b3"])
#     W4 = tf.convert_to_tensor(value=parameters["W4"])
#     b4 = tf.convert_to_tensor(value=parameters["b4"])
#     params = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
#     x = tf.compat.v1.placeholder("float")
#     z4 = forward_propagation_for_predict(x, params)
#     p = tf.argmax(input=z4)
#     sess = tf.compat.v1.Session()
#     prediction = sess.run(p, feed_dict = {x: X})
#     return prediction

def predict_probability(X, parameters,n_cell_Types):
    # input -- X (dataset used to make prediction), papameters after training
    # output -- prediction
    W1 = tf.convert_to_tensor(value=parameters["W1"])
    b1 = tf.convert_to_tensor(value=parameters["b1"])
    W2 = tf.convert_to_tensor(value=parameters["W2"])
    b2 = tf.convert_to_tensor(value=parameters["b2"])
    W3 = tf.convert_to_tensor(value=parameters["W3"])
    b3 = tf.convert_to_tensor(value=parameters["b3"])
    W4 = tf.convert_to_tensor(value=parameters["W4"])
    b4 = tf.convert_to_tensor(value=parameters["b4"])
    params = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
    x = tf.compat.v1.placeholder("float")
    z4 = forward_propagation_for_predict(x, params)
    # z4 = tf.convert_to_tensor(z4, dtype=tf.float32)
    #We need to do this step because otherwise it doesn't work since the 
    z4 = tf.reshape(z4, [n_cell_Types, -1])
    # z4 = tf.reshape(z4, [-1,n_cell_Types])
    p = tf.nn.softmax(z4, axis=0)
    sess = tf.compat.v1.Session()
    prediction = sess.run(p, feed_dict = {x: X})
    return prediction
# Build the model
def model(X_train, Y_train, X_test, starting_learning_rate = 0.0001, num_epochs = 1500, minibatch_size = 128, print_cost = True):
    # input -- X_train (training set), Y_train(training labels), X_test (test set), Y_test (test labels),
    # output -- trained parameters
    ops.reset_default_graph() # to be able to rerun the model without overwriting tf variables
    tf.compat.v1.set_random_seed(3)
    seed = 3
    (nf, ns) = X_train.shape
    nt = Y_train.shape[0]
    costs = []
    # create placeholders of shape (nf, nt)
    X, Y = create_placeholders(nf, nt)
    # initialize parameters
    parameters = initialize_parameters(nf=nf, ln1=100, ln2=50, ln3=25, nt=nt)
    # forward propagation: build the forward propagation in the tensorflow graph
    Z4 = forward_propagation(X, parameters)
    # cost function: add cost function to tensorflow graph
    cost = compute_cost(Z4, Y, parameters, 0.005)
    # Use learning rate decay
    global_step = tf.Variable(0, trainable=False)
    learning_rate = tf.compat.v1.train.exponential_decay(starting_learning_rate, global_step, 1000, 0.95, staircase=True)
    # backpropagation: define the tensorflow optimizer, AdamOptimizer is used.
    optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=learning_rate)
    trainer = optimizer.minimize(cost, global_step=global_step)
    # initialize all the variables
    init = tf.compat.v1.global_variables_initializer()
    # start the session to compute the tensorflow graph
    with tf.compat.v1.Session() as sess:
        # run the initialization
        sess.run(init)
        # do the training loop
        for epoch in range(num_epochs):
            epoch_cost = 0.
            num_minibatches = int(ns / minibatch_size)
            seed = seed + 1
            minibatches = random_mini_batches(X_train, Y_train, minibatch_size, seed)
            for minibatch in minibatches:
                # select a minibatch
                (minibatch_X, minibatch_Y) = minibatch
                # run the session to execute the "optimizer" and the "cost", the feedict contains a minibatch for (X,Y).
                _ , minibatch_cost = sess.run([trainer, cost], feed_dict={X: minibatch_X, Y: minibatch_Y})
                epoch_cost += minibatch_cost / num_minibatches
            # print the cost every epoch
            if print_cost == True and (epoch+1) % 5 == 0:
                print ("Cost after epoch %i: %f" % (epoch+1, epoch_cost))
                costs.append(epoch_cost)
        parameters = sess.run(parameters)
        print ("Parameters have been trained!")
        # calculate the correct predictions
        correct_prediction = tf.equal(tf.argmax(input=Z4), tf.argmax(input=Y))
        # calculate accuracy on the test set
        accuracy = tf.reduce_mean(input_tensor=tf.cast(correct_prediction, "float"))
        print ("Train Accuracy:", accuracy.eval({X: X_train, Y: Y_train}))
        return parameters



if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # load model 
    print('@ LOAD MODEL')
    model = pickle.load(open(args.model_path, 'rb'))
    
    dict_path = os.path.dirname(os.path.abspath(args.model_path))
    label_to_type_dict = pickle.load(open(dict_path + '/label_to_type_dict.pkl', 'rb'))

    print('@ DONE')

    # Read query 
    query = pd.read_csv(args.test_set, index_col=0)

    # transpose query
    query = query.transpose()
    
    # save barcodes
    barcode = list(query.columns)
    
    print("Dimension of the matrix after removing all-zero rows:", query.shape)
    
    # noramlize query 
    query = scale_sets(query)
    
    df = pd.DataFrame(predict_probability(query, model,n_cell_Types = len(label_to_type_dict)))
    df.index = [label_to_type_dict[x] for x in range(df.shape[0])]
    df.columns = barcode
    df = df.T
    # test_predict.to_csv("predicted_probabilities.txt", sep="\t")
    # Create a new column 'max_column' with the column name containing the maximum value for each row
    df['pred_label'], df['proba_label'] = zip(*df.apply(get_max_column_and_value, axis=1))
    df['pred_label_reject'] = df.apply(lambda row: 'Unknown' if row['proba_label'] < args.threshold else row['pred_label'], axis=1)

    print('@ WRITTING PREDICTIONS')
    pred_df = pd.DataFrame({'cell': df.index, "ACTINN": df.pred_label_reject})
    pred_df.to_csv(args.pred_path, index = False)
    print('@ DONE')
    
    # Save the prob matrix
    print('@ WRITTING PROB MATRIX ')
    out_other_path = os.path.dirname(args.pred_path)
    df.to_csv(out_other_path + '/ACTINN_pred_score.csv',
              index=True) #True because we want to conserve the rownames (cells)
    print('@ DONE ')

print("Run time:", timeit.default_timer() - run_time)
