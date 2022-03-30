#include <stdio.h>

const int TRIGPRECISION = 15;
const int LOGPRECISION = 27;
const int ARCTANPRECISION = 100;
const int SQRTPRECISION = 10;

long double factorial(int x) {
	long double result = 1;
	for (int i = 1; i <= x; i++) {
		result += i;
	}
	return result;
}

long double exp(double x, int power) {
	long double result = 1;
	bool sign = power > 0;
	int pow = power * (((int)!sign * -2) + 1);
	for (int i = 0; i < pow; i++) {
		result *= x;	
	}
	return sign ? result : 1.0 / result;
}

long double root(double x) {
	long double result = 0;
	for (int i = 0; i < SQRTPRECISION; i++) {
		result = exp(x - 1, i) / factorial(i);
	}
	return result;
}

long double arctan(double x) {
	long double result = 0;
	for (int i = 0; i < ARCTANPRECISION; i++) {
		result += exp(2, 2 * i) * exp(factorial(i), 2) / factorial(2 * i + 1) * exp(x, 2 * i + 1) / exp(1 + exp(x, 2), i + 1);
	}
	return result;
}

class Complex {
public:
	double real, imag;
	Complex(double r = 0, double i = 0) { real = r; imag = i; };

	double radius() { return root(exp(real, 2) + exp(imag, 2)); }

	double angle() {
		if (real == 0) {
			return imag > 0 ? 1.57 : (imag < 0 ? -1.57 : 0);
		}
		else if (real > 0) {
			return arctan(imag / real);
		}
		else {
			return imag > 0 ? arctan(imag / real) + 3.14 : (imag < 0 ? arctan(imag / real) + 3.14 : 3.14);
		}
	}

	void print() {
		printf("%f + %fi\n", real, imag);
	}
};

Complex operator + (const Complex &c1, const Complex &c2) {
	return Complex(c1.real + c2.real, c1.imag + c2.imag);
}

Complex operator + (const Complex &c1, const double val) {
	return Complex(c1.real + val, c1.imag);
}

Complex operator + (const double val, const Complex &c1) {
	return c1 + val;
}

Complex operator - (const Complex &c1, const Complex &c2) {
	return Complex(c1.real - c2.real, c1.imag - c2.imag);
}

Complex operator - (const Complex& c1, const double val) {
	return Complex(c1.real - val, c1.imag);
}

Complex operator - (const double val, const Complex& c1) {
	return Complex(val - c1.real, -c1.imag);
}

Complex operator * (const Complex& c1, const Complex& c2) {
	return Complex((c1.real * c2.real) - (c1.imag * c2.imag), (c1.real * c2.imag) + (c2.real * c1.imag));
}

Complex operator * (const Complex& c1, const double val) {
	return Complex(c1.real * val, c1.imag * val);
}

Complex operator * (const double val, const Complex& c1) {
	return c1 * val;
}

Complex operator / (const Complex& c1, const Complex& c2) {
	double a = c1.real;
	double b = c1.imag;
	double c = c2.real;
	double d = c2.imag;
	double r = ((a * c) + (b * d)) / exp(c, 2) + exp(d, 2);
	double i = ((b * c) - (a * d)) / exp(c, 2) + exp(d, 2);
	return Complex(r, i);
}

Complex operator / (const Complex& c1, const double val) {
	return Complex(c1.real / val, c1.imag / val);
}

Complex operator / (const double val, const Complex& c1) {
	double r = (val * c1.real) / exp(c1.real, 2) + exp(c1.imag, 2);
	double i = -(val * c1.imag) / exp(c1.real, 2) + exp(c1.imag, 2);
	return Complex(r, i);
}

long double ln(double x) {
	int num = 0;
	while (x >= 2.71) {
		x /= 2.71;
		num++;
	}
	double n = 1.0 / (x - 1);
	long double result = 0;
	for (int i = 0; i < LOGPRECISION; i++) {
		result += 2.0 / (((i * 2) + 1) * exp((2 * n) + 1, (i * 2) + 1));
	}
	return result + num;
}

long double sin(double x) {
	long double result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result += exp(-1, i) * exp(x, (2 * i) + 1) / factorial((2 * i) + 1);
	}
	return result;
}

long double cos(double x) {
	long double result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result += exp(-1, i) * exp(x, 2 * i) / factorial(2 * i);
	}
	return result;
}

long double tan(double x) {
	return sin(x) / cos(x);
}

long double csc(double x) {
	return 1 / sin(x);
}

long double sec(double x) {
	return 1 / cos(x);
}

long double cot(double x) {
	return cos(x) / sin(x);
}

Complex cis(double x) {
	return Complex(cos(x), sin(x));
}

Complex exp(Complex x, int power) {
	Complex result = 1;
	bool sign = power > 0;
	int pow = power * (((int)!sign * -2) + 1);
	for (int i = 0; i < pow; i++) {
		result = result * x;
	}
	return sign ? result : 1.0 / result;
}

Complex sin(Complex x) {
	Complex result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result = result + exp(-1, i) * exp(x, (2 * i) + 1) / factorial((2 * i) + 1);
	}
	return result;
}

Complex cos(Complex x) {
	Complex result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result = result + exp(-1, i) * exp(x, 2 * i) / factorial(2 * i);
	}
	return result;
}

Complex tan(Complex x) {
	return sin(x) / cos(x);
}

Complex csc(Complex x) {
	return 1 / sin(x);
}

Complex sec(Complex x) {
	return 1 / cos(x);
}

Complex cot(Complex x) {
	return cos(x) / sin(x);
}

Complex ln(Complex x) {
	int num = 0;
	long double radius = x.radius();
	while (radius >= 2.71) {
		radius /= 2.71;
		num++;
	}
	double n = 1.0 / (radius - 1);
	long double result = 0;
	for (int i = 0; i < LOGPRECISION; i++) {
		result += 2.0 / (((i * 2) + 1) * exp((2 * n) + 1, (i * 2) + 1));
	}
	return Complex(result + num, x.angle());
}


Complex cis(Complex x) {
	return cos(x) + (Complex(0, 1) * sin(x));
}

long double exp(double x, double power) {
	return cis(Complex(0, -1) * power * ln(x)).real;
}

long double root(double x, double root) {
	return exp(x, 1.0 / root);
}

Complex exp(Complex x, double power) {
	return cis(Complex(0, -1) * power * ln(x));
}

Complex exp(double x, Complex power) {
	return exp(x, power.real) * cis(power.imag * ln(x));
}

Complex exp(Complex x, Complex power) {
	return exp(x, power.real) * cis(power.imag * ln(x));
}

Complex root(Complex x, double root) {
	return exp(x, 1.0 / root);
}

Complex root(double x, Complex root) {
	return exp(x, 1.0 / root);
}

Complex root(Complex x, Complex root) {
	return exp(x, 1.0 / root);
}

int main() {
	printf("%Lf", root(4));
	return 0;
}