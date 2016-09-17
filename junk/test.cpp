#include <iostream>

int main() {
    int sum = 0;
    int k = 10000;
    #pragma omp parallel for reduction(+:sum)
    for(int i=1; i<=k; ++i) sum += i;
    std::cout << sum << std::endl;
    std::cout << k * (k+1) / 2 << std::endl;
    return 0;
}
