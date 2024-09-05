#include "/opt/homebrew/opt/libomp/include/omp.h"
#include <iostream>
#include <cmath>
#include "src/mathomp.cpp"


// Declaration of the isPrimellu function
bool isPrimellu(unsigned long long* pn);

int main() {
    unsigned long long number;
    std::cout << "Enter a number to test for primality: ";
    std::cin >> number;
    
    bool result = isPrimellu(&number);
    
    if (result) {
        std::cout << number << " is a prime number." << std::endl;
    } else {
        std::cout << number << " is not a prime number." << std::endl;
    }

    return 0;
}


