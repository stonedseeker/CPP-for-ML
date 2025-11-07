#include <iostream>
#include <vector>
#include <stdexcept>
#include <ctime>
#include <cmath>
#include <random>
#include <limits>

class KMeans {
private:
    int k;
    int maxIterations;

    std:: vector<std::vector<double>> centroids;

    double euclideanDistanceSquared(const std::vector<double>& point1, const std::vector<double>& point2) {

            double distance = 0.0;
            for (size_t i = 0; i < point1.size(); ++i) {
                distance += std::pow(point1[i] - point2[i], 2);
            }

            return distance;
    }


public:
    KMeans(int numClusters, int maxIter) :
        k(numClusters), maxIterations(maxIter) {}
    

    void fit (const std::vector<std::vector<double>>& X) {
        size_t numSamples = X.size();
        size_t numFeatures = X[0].size();

        centroids.clear()

        std::mt19937 gen(std::time(0));
        std::uniform_int_distribution<> dis(0, numSamples - 1);

        for (int i = 0; i < k; i++) {
            int randIndex = dis(gen);
            centroids.push_back(X[randIndex]);
            
            for (size_t j = 0; j < numFeatures; j++) {
            std::cout << centroids[i][j] << (j == numFeatures - 1 ? "" : ", ");
            }
            std::cout << ")" << std::endl;
        }
    }

    for (int iter = 0; iter < maxIterations; ++iter) {
        std::vector<int> assignments(numSamples, 0);

        for (size_t i = 0; i < numSamples; ++i) {
            double minDistance = std::numeric_limits<double>::max();
            int closestCentroidIndex = 0;

            for (int j = 0; j < k; ++j) {
                double dist = euclideanDistanceSquared(X[i], centroids[j]);

                if (dist < minDistance) {
                    minDistance = dist;
                    closestCentroidIndex = j;
                }
            }
            assignments[i] = closestCentroidIndex;
        }
    }

    std::vector<int> predict(const std::vector<std::vector<double>>& X) {
        std::vector<int> assignments;
        return assignments;
    }
};

int main() {
    std::vector<std::vector<double>> X = {
        {1.0, 2.0},
        {1.5, 1.8},
        {5.0, 8.0},
        {8.0, 8.0},
        {1.0, 0.6},
        {9.0, 11.0}
    };

    KMeans model(2, 100);

    model.fit(X);

    model.predict(X);

    return 0;

}
