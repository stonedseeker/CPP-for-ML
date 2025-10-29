#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

class LogisticRegression {
private:
    std::vector<double> weights;
    double bias;
    double learningRate;
    int numIterations;

    double sigmoid(double z) {
        return 1 / (1 + std::exp(-z));
    }

public:
    LogisticRegression(double lr = 0.01, int iterations = 1000)
        : learningRate(lr), numIterations(iterations), bias(0.0) {}

    void fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
        if (X.empty() || y.empty()) {
    //   throw std::invalid_argument();
         throw std::invalid_argument("Input feature matrix (X) or target vector (y) is empty.");
        }

        size_t numSamples = X.size();
        size_t numFeatures = X[0].size();

        weights.assign(numFeatures, 0.0);

        bias = 0.0;

        std::cout << numSamples << std::endl;
        std::cout << numFeatures << std::endl;

        // Gradient Descent
        for (int iter = 0; iter < numIterations; ++ iter) {
            std:: vector<double> dW(numFeatures, 0.0);
            double db = 0.0;

            // calculate the total gradient
            for (size_t i = 0; i < numSamples; i++) {
                double z = 0.0; // a. first make a prediction 
                for (size_t j = 0; i < numFeatures; j++) {
                    z += weights[j] * X[i][j]; // y = wx + c;
                                               // this is the w.x part
                }

                z += bias; // this is the + c part

                double h = sigmoid(z);

                double error = h - y[i]; // b. calculate the error

                // Accumulate the gradients
                for (size_t j = 0; j < numFeatures; j++) {
                    dW[j] += error * X[i][j];
                }

                db += error;
            }

            // Average the gradients
            for (size_t i = 0; i < numSamples; i++) {
                dW[i] /= numSamples;
            }

            db /= numSamples;

            // Go downhill
            for (size_t i = 0; i < numSamples; i++) {
                weights[i] -= dW[i] * learningRate;
            }

            bias -= db * learningRate;



        }

    }

    std::vector<double> predict_probaility(const std::vector<std::vector<double>>& X) {
        size_t numSamples = X.size();
        size_t numFeatures = X[0].size();

        std::vector<double> probabilities;

        probabilities.reserve(numSamples);

        for (size_t i = 0; i < numSamples; i++) {
            double z = 0.0;

            for (size_t j = 0; j < numFeatures; j++) {
                z += weights[j] * X[i][j];
            }

            z += bias;

            probabilities.push_back(sigmoid(z));
        }

        return probabilities;
    }

    std::vector<int> predict(const std::vector<std::vector<double>>& X) {
        std::vector<double> probabilities = predict_probaility(X);

        std::vector<int> predictions;

        predictions.reserve(probabilities.size());

        for (double probability : probabilities) {
            if (probability >= 0.5) {
                predictions.push_back(1);
            } else {
                predictions.push_back(0);
            }
        }

        return predictions;
    }
};


int main() {
    std:: vector<std::vector<double>> X = {
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0}
    };

    std:: vector<double> y = {0.0, 1.0, 1.0, 1.0}; // A Simple OR Gate 
    LogisticRegression model;

    try {
        model.fit(X, y);
    } catch (const std:: exception& e) {
        std:: cerr << e.what() << std::endl;
        return 1;
    }

    try {
        std::vector<int> predictions = model.predict(X);
        for (size_t i = 0; i < X.size(); i++) {
            std::cout << predictions[i] << std::endl;
        }
    }catch (const std:: exception& e) {
        std::cerr << "An error occured" << e.what() << std::endl;
        return 1;
    }



    return 0;
}
