#include <iostream>
#include <vector>
#include <cmath>

using namespace std; 

double predict(double slope, double weight, double bias) {
    return slope * weight + bias;
}

double compute_loss(const vector<double>& Weights, const vector<double>& Predictions, double slope, double bias) {
    double total_error = 0.0;
    int n = Weights.size();

    int i = 0;
    for (; i < n; i++) {
        double prediction = predict(slope, Weights[i], bias);
        total_error += std :: pow(Predictions[i] - prediction, 2);
    }

    return total_error / n; // Mean Squared Erorr
}

void gradient_descent(const vector<double>& Weights, const vector<double>& Predictions, double& slope, double& bias, double& lr) {
    double slope_gradient = 0.0;
    double bias_gradient = 0.0;

    int n = Weights.size();

    int i = 0;
    for (; i < n; i++) {
        double prediction = predict(slope, Weights[i], bias);
        double error = Predictions[i] - prediction; 
        slope_gradient += (-2.0 / n) * Weights[i] * error;
        bias_gradient += (-2.0 / n) * error;
    }

    slope -= lr * slope_gradient;
    bias -= lr * bias_gradient;
}

int main() {

    cout << "Linear Regression from Scratch" << endl;

    // Examplery Data 
    vector<double> Weights = {1,2,3,4,5,6,7,8,9,10};
    vector<double> Predictions = {2,4,6,8,10,12,14,16,18,20};

    double slope = 0.0;
    double bias = 0.0;
    double lr = 0.01;
    int epochs = 1000;

    for (int i = 0; i < epochs; i++) {
        gradient_descent(Weights, Predictions, slope, bias, lr);

        if (i % 100 == 0) {
            double loss = compute_loss(Weights, Predictions, slope, bias);
            cout << "Epoch " << i << " | Loss: " << loss << " | m: " << slope << " | b: " << bias << endl;
        }
    }

    cout << "\nFinal model: Prediction = " << slope << "Weights + " << bias << endl;

    // Make a Prediction
    double test_weight = 6;
    cout << "Prediction for x = " << test_weight << ": " << predict(slope, test_weight, bias) << endl;
    return 0;

}
