/*c++ test primatlity tests
Only 3 working so far
basic, miller-rabbin, c++ 17 other*/
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
    std::chrono::duration_cast<std::chrono::microseconds>( \
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

/*Modulo used for Solovay-Strassen*/
unsigned long long modulo(unsigned long long base, unsigned long long exponent, unsigned long long mod) { 
    unsigned long long x = 1; 
    unsigned long long y = base; 
    while (exponent > 0) 
    { 
        if (exponent % 2 == 1) 
            x = (x * y) % mod; 
  
        y = (y * y) % mod; 
        exponent = exponent / 2; 
    } 
  
    return x % mod; 
} 

/* To calculate Jacobian symbol of a given number */
unsigned long long calculateJacobian(unsigned long long a, unsigned long long n){ 
    if (!a) 
        return 0;// (0/n) = 0 

    int ans = 1; 
	
    if (a < 0)  { 
        a = -a; // (a/n) = (-a/n)*(-1/n) 
        if (n % 4 == 3) 
            ans = -ans; // (-1/n) = -1 if n = 3 (mod 4) 
    }
  
    if (a == 1) 
        return ans;// (1/n) = 1 
  
    while (a) { 
        if (a < 0) { 
            a = -a;// (a/n) = (-a/n)*(-1/n) 
            if (n % 4 == 3) 
                ans = -ans;// (-1/n) = -1 if n = 3 (mod 4) 
        } 
  
        while (a % 2 == 0){ 
            a = a / 2; 
            if (n % 8 == 3 || n % 8 == 5) 
                ans = -ans; 
  
        } 
  
        std::swap(a, n); 
  
        if (a % 4 == 3 && n % 4 == 3) 
            ans = -ans; 
        a = a % n; 
  
        if (a > n / 2) 
            a = a - n; 
  
    } 
  
    if (n == 1) 
        return ans; 
  
    return 0; 
} 

/*Solovay-Strassen*/
bool solovoyStrassen(unsigned long long p, int iterations) { 
    if (p < 2) 
        return false; 
    if (p != 2 && p % 2 == 0) 
        return false; 
  
    for (int i = 0; i < iterations; i++){ 
        // Generate a random number a 
        unsigned long long a = rand() % (p - 1) + 1; 
        unsigned long long jacobian = (p + calculateJacobian(a, p)) % p; 
        unsigned long long mod = modulo(a, (p - 1) / 2, p); 
  
        if (!jacobian || mod != jacobian) 
            return false; 
    } 
    return true; 
} 

// function to calculate the coefficients 
// of (x - 1)^n - (x^n - 1) with the help 
// of Pascal's triangle . 
void coef(unsigned long long n,  std::vector< long long>& c) 
{ 
    c[0]=1; 
    for (long long i = 0; i < n; c[0] = -c[0], i++) { 
        c[1 + i] = 1; 
  
        for ( long long  j = i; j > 0; j--) 
            c[j] = c[j - 1] - c[j]; 
    } 
} 
  
// function to check whether 
// the number is prime or not 
bool AKSPrime(unsigned long long n,  std::vector< long long>& c) 
{ 
    // Calculating all the coefficients by 
    // the function coef and storing all 
    // the coefficients in c array . 
    coef(n, c); 
  
    // subtracting c[n] and adding c[0] by 1 
    // as ( x - 1 )^n - ( x^n - 1), here we 
    // are subtracting c[n] by 1 and adding 
    // 1 in expression. 
    c[0]++, c[(long long)n]--; 
  
    // checking all the coefficients whether 
    // they are divisible by n or not. 
    // if n is not prime, then loop breaks 
    // and (i > 0). 
    long long i = n; 
    while (i-- && c[i] % n == 0) 
        ; 
  
    // Return true if all coefficients are 
    // divisible by n. 
    return i < 0; 
} 

void getResult(bool x) {
	if (x) {
		std::cout << "test successful - number IS PRIME!!" << std::endl;
	} else {
		std::cout << "test result: NOT prime OR if Miller|SolovoyStrassen maybe unclear?" << std::endl;

	}
	return;
}


/*hack using tgamma function (n-1)!^2 mod n'equals 1 if n is prime and 0 if it is not
*/
bool tgammaPrime(unsigned long long n){
    //'
    unsigned long long nt = std::pow(std::tgamma(n), 2);
	//std::cout << "nt= " << nt << std::endl;
    nt %= n;
	//std::cout << "nt= " << nt << std::endl;
    if (nt == 1.0)
        return true;
    return false;

}


int main() {

	/*precision values*/
	int k = 13;
    std::vector< long long> c;
	
	INIT_TIMER
	unsigned long long number;
	std::cout << "Enter a long long int for primtality testing: ";
	std::cin >> number;

	std::cout << std::endl;
	std::cout << "Beginning primality testing for n = " << number << " ---" << std::endl;
	std::cout << "ms results represent microseconds = " << number << " ---" << std::endl;

	START_TIMER
	getResult(bisprime(number));
	STOP_TIMER("for basic isprime ")
	std::cout << std::endl;

	/*START_TIMER
	getResult(millerisPrime(number, k));
	STOP_TIMER( "for Miller-Rabbin isPrime ")
	std::cout << std::endl;
	*/

	START_TIMER
	getResult(c17Prime(number));
	STOP_TIMER (" for cpp modern isPrime " )
	std::cout << std::endl;

	START_TIMER
	getResult(tgammaPrime(number));
	STOP_TIMER (" for cpp tgamma isPrime " )
	std::cout << std::endl;

}