/*c++ test primatlity tests*/
#include <iostream>
#include <limits>
#include <random>
#include <numeric>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <vector>
#include <iterator>
#include <chrono>

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-start \
    ).count() << " ms " << std::endl; 
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

// Utility function to do modular exponentiation. 
// It returns (x^y) % p 
unsigned long long power(unsigned long long x, unsigned int y, long long p)
{
	unsigned long long res = 1;      // Initialize result 
	x = x % p;  // Update x if it is more than or 
				// equal to p 
	while (y > 0)
	{
		// If y is odd, multiply x with result 
		if (y & 1)
			res = (res * x) % p;

		// y must be even now 
		y = y >> 1; // y = y/2 
		x = (x * x) % p;
	}
	return res;
}

// This function is called for all k trials. It returns 
// false if n is composite and returns false if n is 
// probably prime. 
// d is an odd number such that  d*2<sup>r</sup> = n-1 
// for some r >= 1 
bool miillerTest(unsigned long long& d, unsigned long long n)
{
	// Pick a random number in [2..n-2] 
	// Corner cases make sure that n > 4 
	unsigned long long a = 2 + rand() % (n - 4);

	// Compute a^d % n 
	unsigned long long x = power(a, d, n);

	if (x == 1 || x == n - 1)
		return true;

	// Keep squaring x while one of the following doesn't 
	// happen 
	// (i)   d does not reach n-1 
	// (ii)  (x^2) % n is not 1 
	// (iii) (x^2) % n is not n-1 
	while (d != n - 1)
	{
		x = (x * x) % n;
		d *= 2;

		if (x == 1)      return false;
		if (x == n - 1)    return true;
	}

	// Return composite 
	return false;
}

// It returns false if n is composite and returns true if n 
// is probably prime.  k is an input parameter that determines 
// accuracy level. Higher value of k indicates more accuracy. 
bool millerisPrime(unsigned long long n, int k)
{
	// Corner cases 
	if (n <= 1 || n == 4)  return false;
	if (n <= 3) return true;

	// Find r such that n = 2^d * r + 1 for some r >= 1 
	unsigned long long d = n - 1;
	while (d % 2 == 0)
		d /= 2;

	// Iterate given nber of 'k' times 
	for (int i = 0; i < k; i++)
		if (!miillerTest(d, n))
			return false;

	return true;
}

/*modernc++17 based found*/
bool c17Prime(unsigned long long const n) {
	if (n <= 3) return n > 1;
	else if (n % 2 == 0 || n % 3 == 0) return false;
	else
	{
		for (long long i = 4; i * i <= n; i += 6)
		{
			if (n % i == 0 || n % (i + 2) == 0)
			{
				return false;
			}
		}
		return true;
	}
}

/*basic validates if integer n is prime number
returns true or false*/
bool bisprime(unsigned long long n) {
	bool k = true;
	if (n != 2) {
		for (int i = 2; i < (int)std::sqrt(n) + 1; i++) {
			if (n % i == 0) {
				k = false;
				break;
			}
		}
	}
	return k;
}


void getResult(bool x) {
	if (x) {
		std::cout << "test successful - number is prime" << std::endl;
	} else {
		std::cout << "test result: NOT prime OR unclear Miller-Rabbin" << std::endl;

	}
	return;
}


int main() {

	int k = 8;
	INIT_TIMER
	unsigned long long number;
	std::cout << "Enter a long long int for primtality testing: ";
	std::cin >> number;

	std::cout << std::endl;
	std::cout << "Beginning primality testing for n = " << number << " ---" << std::endl;

	START_TIMER
	getResult(bisprime(number));
	STOP_TIMER("for basic isprime ")
	std::cout << std::endl;

	START_TIMER
	getResult(millerisPrime(number, k));
	STOP_TIMER( "for Miller-Rabbin isPrime ")
	std::cout << std::endl;

	START_TIMER
	getResult(c17Prime(number));
	STOP_TIMER (" for cpp modern isPrime " )
	std::cout << std::endl;


}