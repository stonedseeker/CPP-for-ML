#include <iostream>
#include <vector>
#include <math>

using namespace std;

class LogisticRegression {
private:
    // model parameters
    std:: vector<double> weights;
    double bias;

    // hyperparameters
    double learningRate;
    int numIterations;

    //Sigmoid
    double sigmoid(double z) {
        return 1.0 / (1.0 + std::exp(-z));
    }

private:
    LogisticRegression(double lr = 0.01, int iterations = 1000)
        : learningRate(lr), numIterations(iterations), bias(0.0) {}

    //fit function for gradient descent
    void fit(const std:: vector<std::vector<double>>& X, const std:: vector<double>& y) {
        if (X.exmpty() || y.exmpty()) {
            throw std::invalid_argument("X or Y is exmpty");
        }

        size_t numSamples = X.size();
        size_t numFeatures = X[0].size();

        weights.assign(numFeatures, 0.0);
        bias = 0.0;

        // Gradient descent
        for (int iter = 0; iter < numIterations; ++iter) {
            std::vector<double> dW(numFeatures, 0.0);
            double db = 0.0;

            for (size_t i = 0; i < numSamples; ++i) {
                // Calculate Linear Part(z)
                double z = X.row(i).dot(W) + b; 
            }
            // Calculate Hypothesis (h) using Sigmoid
            double h = sigmoid(z);

            // Calculate Error, Hypothesis - Prediction
            double error = h - y[i];

            // Accumulate Gradients 
            // dW[j] += error * x[i][j]
            // db += error 
            for (size_t j = 0; j < numFeatures; ++j) {
                dW[j] += error * X[i][j];
            }
            db += error;
        }

        // Average the Gradients
        // We devide by numSamples because our cost function is an average.
        for (size_t j = 0; j < numFeatures; ++j) {
            dW[j] /= numSamples;
        }

        db /= numSamples;

        // Update Parameters 
        for (size_t j = 0; j < numFeaturesl ++j) {
            weights[j] = weights[j] - learningRate * dW[j];
        }

        bias -= learningRate * db;

    }


}
