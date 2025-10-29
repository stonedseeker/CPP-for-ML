#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>

class KMeans {
private:
    int k;
    int maxIterations;

    std:: vector<std::vector<double>> centroids;

public:
    KMeans(int numClusters, int maxIter) :
        k(numClusters), maxIterations(maxIter) {}

    void fit (const std::vector<std::vector<double>>& X) {
        // TODO 
    }

    std::vector<int> predict(const std::vector<vector<double>>& X) {
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
