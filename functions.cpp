#include <stdio.h>

const int TRIGPRECISION = 15;
const int LOGPRECISION = 27;
const int ARCTANPRECISION = 100;
const int SQRTPRECISION = 100;
const int EPRECISION = 26;

long double factorial(int x) {
	long double result = 1;
	for (int i = 1; i <= x; i++) {
		result *= i;
	}
	return result;
}

long double exp(long double x, int power) {
	long double result = 1;
	bool sign = power > 0;
	int pow = power * (((int)!sign * -2) + 1);
	for (int i = 0; i < pow; i++) {
		result *= x;	
	}
	return sign ? result : 1.0 / result;
}

long double root(long double x) {
  double upper = x;
  double lower = 0;
  double middle;
  for (int i = 0; i < SQRTPRECISION; i++) {
    middle = (upper + lower) / 2;
    if (exp(middle, 2) > x) {
      upper = middle;
    }
    else if (exp(middle, 2) < x) {
      lower = middle;
    }
    else if (exp(middle, 2) == x) { break; }
  }
  return middle;
}

long double arctan(long double x) {
	long double result = 0;
	for (int i = 0; i < ARCTANPRECISION; i++) {
	  result += exp(2, 2 * i) * 
      exp(factorial(i), 2) / 
      factorial(2 * i + 1) * 
      exp(x, 2 * i + 1) / 
      exp(1 + exp(x, 2), i + 1);
	}
	return result;
}

long double pi(){
  return arctan(1) * 4;
}

long double e(){
  long double result = 0;
  for (int i = 0; i < EPRECISION; i++){
    result += 1/factorial(i);
  }
  return result;
}

const long double PI = pi();
const long double E = e();

class Complex {
public:
	long double real, imag;
	Complex(long double r = 0, long double i = 0) { real = r; imag = i; }

  //long double get_value() { return real; }

	long double radius() { return root(exp(real, 2) + exp(imag, 2)); }

	long double angle() {
		if (real == 0) {
			return imag > 0 ? PI/2 : (imag < 0 ? -PI/2 : 0);
		}
		else if (real > 0) {
			return arctan(imag / real);
		}
		else {
			return imag > 0 ? arctan(imag / real) + PI : (imag < 0 ? arctan(imag / real) + PI : PI);
		}
	}

	void print() {
		printf("%Lf + %Lfi\n", real, imag);
	}
};

Complex operator + (const Complex &c1, const Complex &c2) {
	return Complex(c1.real + c2.real, c1.imag + c2.imag);
}

Complex operator + (const Complex &c1, const long double val) {
	return Complex(c1.real + val, c1.imag);
}

Complex operator + (const long double val, const Complex &c1) {
	return c1 + val;
}

Complex operator - (const Complex &c1, const Complex &c2) {
	return Complex(c1.real - c2.real, c1.imag - c2.imag);
}

Complex operator - (const Complex& c1, const long double val) {
	return Complex(c1.real - val, c1.imag);
}

Complex operator - (const long double val, const Complex& c1) {
	return Complex(val - c1.real, -c1.imag);
}

Complex operator * (const Complex& c1, const Complex& c2) {
	return Complex((c1.real * c2.real) - (c1.imag * c2.imag), (c1.real * c2.imag) + (c2.real * c1.imag));
}

Complex operator * (const Complex& c1, const long double val) {
	return Complex(c1.real * val, c1.imag * val);
}

Complex operator * (const long double val, const Complex& c1) {
	return c1 * val;
}

Complex operator / (const Complex& c1, const Complex& c2) {
	long double a = c1.real;
	long double b = c1.imag;
	long double c = c2.real;
	long double d = c2.imag;
	long double r = ((a * c) + (b * d)) / (exp(c, 2) + exp(d, 2));
	long double i = ((b * c) - (a * d)) / (exp(c, 2) + exp(d, 2));
	return Complex(r, i);
}

Complex operator / (const Complex& c1, const long double val) {
	return Complex(c1.real / val, c1.imag / val);
}

Complex operator / (const long double val, const Complex& c1) {
	long double r = (val * c1.real) / (exp(c1.real, 2) + exp(c1.imag, 2));
	long double i = -(val * c1.imag) / (exp(c1.real, 2) + exp(c1.imag, 2));
	return Complex(r, i);
}

long double ln(long double x) {
	int num = 0;
	while (x >= E) {
		x /= E;
		num++;
	}
	long double n = 1.0 / (x - 1);
	long double result = 0;
	for (int i = 0; i < LOGPRECISION; i++) {
		result += 2.0 / (((i * 2) + 1) * exp((2 * n) + 1, (i * 2) + 1));
	}
	return result + num;
}

Complex negln(long double x){
  return Complex(ln(-x), PI);
}

long double log(long double x, long double base){
  return ln(x)/ln(base);
}

Complex neglog(long double x, long double base){
  bool signX = x >= 0;
  bool signBase = base >= 0;
  if (signX){
    if (signBase){
      return ln(x)/ln(base);
    }
    else{
      return ln(x)/negln(base);
    }
  }
  else{
    if (signBase){
      return negln(x)/ln(x);
    }
    else{
      return negln(x)/negln(base);
    }
  }
}

long double sin(long double x) {
	long double result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result += exp(-1, i) * exp(x, (2 * i) + 1) / factorial((2 * i) + 1);
	}
	return result;
}

long double cos(long double x) {
	long double result = 0;
	for (int i = 0; i < TRIGPRECISION; i++) {
		result += exp(-1, i) * exp(x, 2 * i) / factorial(2 * i);
	}
	return result;
}

long double tan(long double x) {
	return sin(x) / cos(x);
}

long double csc(long double x) {
	return 1 / sin(x);
}

long double sec(long double x) {
	return 1 / cos(x);
}

long double cot(long double x) {
	return cos(x) / sin(x);
}

Complex cis(long double x) {
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
	while (radius >= E) {
		radius /= E;
		num++;
	}
	long double n = 1.0 / (radius - 1);
	long double result = 0;
	for (int i = 0; i < LOGPRECISION; i++) {
		result += 2.0 / (((i * 2) + 1) * exp((2 * n) + 1, (i * 2) + 1));
	}
	return Complex(result + num, x.angle());
}

Complex log(Complex x, long double base){
  return ln(x)/ln(base);
}

Complex log(Complex x, Complex base){
  return ln(x)/ln(base);
}

Complex cis(Complex x) {
	return cos(x) + (Complex(0, 1) * sin(x));
}

long double exp(long double x, long double power) {
	if (x < 0) {
    return cis(Complex(0, -1) * power * negln(x)).real;
  }
  else {
    return cis(Complex(0, -1) * power * ln(x)).real;
  }
}

long double root(long double x, long double root) {
	return exp(x, 1.0 / root);
}

Complex exp(Complex x, long double power) {
	return cis(Complex(0, -1) * power * ln(x));
}

Complex exp(long double x, Complex power) {
	if (x < 0) {
    return exp(x, power.real) * cis(power.imag * negln(x));
  }
  else {
    return exp(x, power.real) * cis(power.imag * ln(x)); 
  }
}

Complex exp(Complex x, Complex power) {
	return exp(x, power.real) * cis(power.imag * ln(x));
}

Complex root(Complex x, long double root) {
	return exp(x, 1.0 / root);
}

Complex root(long double x, Complex root) {
	return exp(x, 1.0 / root);
}

Complex root(Complex x, Complex root) {
	return exp(x, 1.0 / root);
}

Complex arcsin(long double x){
  return Complex(0, -1) * ln(root(1-exp(x, 2))) + (Complex(0, 1) * x);
}

Complex arcsin(Complex x){
  return Complex(0, -1) * ln(root(1-exp(x, 2), 2)) + (Complex(0, 1) * x);
}

Complex arccos(long double x){
  return PI/2 + (Complex(0, 1) * ln(root(1-exp(x, 2))) + (Complex(0, 1) * x));
}

Complex arccos(Complex x){
  return PI/2 + (Complex(0, 1) * ln(root(1-exp(x, 2), 2)) + (Complex(0, 1) * x));
}

int main() {
  printf("%Lf", exp(2.2, 2.3));
	return 0;
}

// It will look like water