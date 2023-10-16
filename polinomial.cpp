#pragma once
#include <iostream>
#include<Windows.h>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <math.h>
#define M_PI (3.141592653589793)    
#define M_2PI (2.*M_PI)

using namespace std;
template <typename T>


class Polinomial {
	T* _coeff;
	int _degree;

	static constexpr float ACCURACY = 0.01f;

public:
	Polinomial();
	Polinomial(const T* coeff, int degree);
	Polinomial(const complex<T> coeff, int degree);
	Polinomial(const int max_degree);
	Polinomial(const Polinomial<T>& other);
	Polinomial(const Polinomial<complex<T>>& other);
	~Polinomial();

	void set_coeff(const T coeff, int degree);
	void set_coeff(double re, double im, int degree);
	void set_degree(const int degree);
	T get_coeff(int degree) const;
	int get_degree() const;

	Polinomial<T>& operator+(const Polinomial<T>& other) const;
	Polinomial<complex<T>>& operator+(const Polinomial<complex<T>>& other) const;
	Polinomial<complex<T>>& operator+=(const Polinomial<complex<T>>& other);
	Polinomial<T>& operator+=(const Polinomial<T>& other);
	Polinomial<T>& operator-(const Polinomial<T>& other) const;
	Polinomial<complex<T>>& operator-(const Polinomial<complex<T>>& other) const;
	Polinomial<T>& operator-=(const Polinomial<T>& other);
	Polinomial<complex<T>>& operator-=(const Polinomial<complex<T>>& other);
	Polinomial<T>& operator*(double value);
	Polinomial<complex<T>>& multiplus(double value);
	T& operator[](const int degree);

	bool operator==(const Polinomial<T>& other) const;
	bool operator==(const Polinomial<complex<T>>& other) const;
	bool operator!=(const Polinomial<T>& other)const;
	bool operator!=(const Polinomial<complex<T>>& other) const;

	bool shrink_to_fit();
	bool expand(const int degree);
	T calculate(const T x) const;
	complex<T> calculate(const complex<T> x) const;

};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Polinomial<T> poly);


template<typename T>
void Polinomial<T>::set_coeff(const T coeff, int degree) {
	if (degree < 0 || degree > _degree) {
		throw std::runtime_error("Degree must be positive and less, that degree polimomial");
	}
	if (degree == 0) {
		_coeff[degree] = 1;
	}
	_coeff[degree] = coeff;
}

template<typename T>
void Polinomial<T>::set_coeff( double re, double im, int degree)
{    
	if (degree < 0 || degree > _degree) {
		throw std::runtime_error("Degree must be positive and less, that degree polimomial");
	}
	if (degree == 0) {
		_coeff[degree].real(1);
		_coeff[degree].imag(1);
	}
	_coeff[degree].real(re);
	_coeff[degree].imag(im);
}

template<typename T>
void Polinomial<T>::set_degree(const int degree) {
	if (degree < 0 || degree > _degree) {
		throw std::runtime_error("Degree must be positive and less, that degree polimomial");
	}
	this->_degree = degree;
}

template <typename T>
int Polinomial<T>::get_degree() const {
	return _degree;
}

template <typename T>
T Polinomial<T>::get_coeff(int deg) const {
	if (deg <= _degree) {
		return _coeff[deg];
	}
	return T();
}

//три конструктора
template<typename T>
Polinomial<T>::Polinomial() : _degree(0) {
	_coeff = new T[1];
	_coeff[0] = 0;
}

template<typename T>
Polinomial<T>::Polinomial(const T* coeff, int degree) : _degree(degree) {
	_coeff = new T[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = coeff[i];
	};
}

template<typename T>
 Polinomial<T>::Polinomial(const complex<T> coeff, int degree) : _degree(degree) {
	_coeff = new complex<T>[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i]->real(coeff[i].real());
		_coeff[i]->imag(coeff[i].imag());
	}
}

template<typename T>
Polinomial<T>::Polinomial(const Polinomial<T>& other) {
	_degree = other.get_degree();
	_coeff = new T[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = other._coeff[i];
	}
}

template<typename T>
Polinomial<T>::Polinomial(const Polinomial<complex<T>>& other) {
	_degree = other.get_degree();
	_coeff = new complex<T>[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i]->real(_coeff[i].real());
		_coeff[i]->imag(_coeff[i].imag());
	};
}

template<typename T>
Polinomial<T>::Polinomial(const int max_degree) : _degree(max_degree) {
	_coeff = new T[_degree + 1];
	_coeff[0] = 1;
	for (int i = 1; i <= _degree; ++i) {
		_coeff[i]=0;
	}
}

template<typename T>
Polinomial<T>::~Polinomial() {
	delete[]  _coeff;
	_coeff = nullptr;
	_degree = 0;
}

template<typename T>
bool Polinomial<T>::expand(const int degree) {
	if (degree <= _degree)
		throw std::runtime_error("Degree must be positive and less, that degree polimomial");
	T* new_coeff = new T[degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		new_coeff[i] = _coeff[i];
		
	}
	_degree = degree;
	delete[] _coeff;
	_coeff = new_coeff;
	return true;
}

template <typename T>
Polinomial<T>& Polinomial<T>::operator+=(const Polinomial<T>& other) {
	if (_degree >= other._degree) {

		for (int i = 0; i <= other._degree; i++)
		{
			_coeff[i] += other._coeff[i];
		}
		

	}
	else {
		int r = _degree;
		for (int i = 0; i <= _degree; ++i) {
			_coeff[i] += other._coeff[i];
		}
		this->expand(other._degree); 
		for (int i = 0; i <= _degree; ++i) {
			if (i > r) {
				_coeff[i]= other._coeff[i];
			}
		}


	}return *this;
	
}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>::operator+=(const Polinomial<complex<T>>& other) {
	if (_degree >= other._degree) {
		for (int i = 0; i <= other._degree; i++) {
			_coeff[i].real() += other._coeff[i].real();
			_coeff[i].imag() += other._coeff[i].imag();
		}
	}
	else {
		int r = _degree;
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			_coeff[i].real() += other._coeff[i].real();
			_coeff[i].imag() += other._coeff[i].imag();
		}
		this->expand(other._degree);
		for (int i = 0; i <= _degree; ++i) {
			if (i > r) {
				_coeff[i].real() = other._coeff[i].real();
				_coeff[i].imag() = other._coeff[i].imag();
			}
		}
	}
	return *this;

}

template <typename T>
Polinomial<T>& Polinomial<T>:: operator+(const Polinomial<T>& other) const {
	Polinomial<T>* result = new Polinomial<T>(*this);
	*result += other;
	return *result;
}


template <typename T>
Polinomial<complex<T>>& Polinomial<T>:: operator+(const Polinomial<complex<T>>& other) const {
	Polinomial<complex<T>>* result = new Polinomial<complex<T>>(*this);
	*result += other;
	return *result;
}

template <typename T>
Polinomial<T>& Polinomial<T>:: operator-=(const Polinomial<T>& other) {
	if (_degree >= other._degree) {
		for (int i = 0; i <= other._degree; i++)
		{
			_coeff[i] -= other._coeff[i];
		}
	}
	else {
		int r = _degree;
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			other._coeff[i]-=_coeff[i];
		}
		this->expand(other._degree);
		for (int i = 0; i <= _degree; ++i) {
			if (i > r) {
				_coeff[i] = other._coeff[i];
			}
		}
	}
	this->shrink_to_fit();
	return *this;

}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>:: operator-=(const Polinomial<complex<T>>& other) {
	if (_degree <= other._degree) {
		for (int i = 0; i <= other._degree; i++) {
			_coeff[i].real() -= other._coeff[i].real();
			_coeff[i].imag() -= other._coeff[i].imag();
		}
	}
	else {
		int r = _degree;
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			other._coeff[i].real() -= _coeff[i].real();
			other._coeff[i].imag()-=_coeff[i].imag();
		}
		this->expand(other._degree);
		for (int i = 0; i <= _degree; ++i) {
			if (i > r) {
				other._coeff[i].real() = _coeff[i].real();
				other._coeff[i].imag() = _coeff[i].imag();
			}
		}
	}
	this->shrink_to_fit();
	return *this;

}

template <typename T>
Polinomial<T>& Polinomial<T>:: operator-(const Polinomial<T>& other) const {
	Polinomial<T>* result = new Polinomial<T>(*this);
	*result -= other;
	return *result;
}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>:: operator-(const Polinomial<complex<T>>& other) const {
	Polinomial<complex<T>>* result = new Polinomial<complex<T>>(*this);
	*result -= other;
	return *result;
}

template <typename T>
Polinomial<T>& Polinomial<T>::operator*(double value) {
	Polinomial<T>* result = new Polinomial<T>(*this);
	for (int i = 0; i <= _degree; ++i) {
		result->_coeff[i] = value * _coeff[i];
	}
	return *this;
}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>::multiplus(double value) {
	Polinomial<complex <T>>* result = new Polinomial <complex<T>>(*this);
	for (int i = 0; i <= _degree; ++i) {
		result[i].real()=_coeff[i].real()* value;
		result[i].imag() = _coeff[i].imag() * value;
		
	}
	return *result;
}

template <typename T>
T& Polinomial<T>::operator[](const int degree) {
	if (degree < 0 || degree > _degree) {
		throw std::runtime_error("Degree must be positive and less, that degree polimomial");
	}
	return _coeff[degree];
}


template<typename T>
T Polinomial<T>::calculate(const T x) const {
	T result = 0;
	for (int i = 0; i <= _degree; ++i) {
		result += _coeff[i] * pow(x, i);
	}
	return result;
}

template<typename T>
complex<T> Polinomial<T>::calculate(const complex<T> x) const {
	complex<T> result = 0;
	for (int i = 0; i <= _degree; ++i) {
		result.real() += _coeff[i].real() * pow(x.real(), i);
		result.imag() += _coeff[i].imag() * pow(x.real(), i);
	}
	return result;
}

template<typename T>
bool Polinomial<T>::shrink_to_fit() {
	int count = 0;
	T zero = { 0 };
	for (int i = _degree; i > 0; --i) {
		if (_coeff[i] == zero)
			++count;
		else break;
	}
	if (count == 0)
		return false;

	T* new_coeff = new T[_degree - count + 1];
	_degree -= count;
	for (int i = 0; i <= _degree; ++i) {
		new_coeff[i] = _coeff[i];
	}
	delete[] _coeff;
	_coeff = new_coeff;

	return true;
}



template<typename T>
bool Polinomial<T>::operator==(const Polinomial<T>& other) const {

	if (_degree != other._degree)
		return false;

	for (int i = 0; i <=_degree; ++i) {
		if (fabs(_coeff[i] - other._coeff[i]) > ACCURACY)
			return false;
	}
	return true;
}

template<typename T>
bool Polinomial<T>::operator==(const Polinomial<complex<T>>& other) const {
	
	if (_degree != other._degree)
		return false;

	for (int i = 0; i <= this._degree; ++i) {
		if ((fabs(_coeff[i].real() - other._coeff[i].real()) > ACCURACY) && (fabs(_coeff[i].imag() - other._coeff[i].imag())> ACCURACY))
			return false;
	}
	return true;
}

template<typename T>
bool Polinomial<T>::operator!=(const Polinomial<T>& other) const {
	if (*this == other)
		return false;
	return true;

}

template<typename T>
bool Polinomial<T>::operator!=(const Polinomial<complex<T>>& other) const {
	if (*this == other)
		return false;
	return true;

}



template<typename T>
std::ostream& operator<<(std::ostream& os, const Polinomial<complex<T>>& poly) {
	for (int i = poly.get_degree(); i >= 1; --i) {
		os <<"(" << poly.get_coeff(i).real() << "+" << poly.get_coeff(i).imag() <<"i"<<")x^" << i << "+";
		if ((poly.get_coeff(i).real() == 0) && (poly.get_coeff(i).imag() == 0))
			continue;
	}
	os << "(" << poly.get_coeff(0).real() << "+" << poly.get_coeff(0).imag() << "i" << ")" << endl;



	return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Polinomial<T> poly) {

	for (int i = poly.get_degree(); i >= 1; --i) {
		os << poly.get_coeff(i) << "x^" << i << "+";
		if (poly.get_coeff(i) == 0)
			continue;
	}
	os << poly.get_coeff(0) << endl;
		
	

	return os;
}

template<typename T>
T find_zero( T a, T b, T c, T d) {
	T p, q;
	std::complex< double > w2(-0.5, sqrt(3)/2);
	std::complex< double > w3(-0.5, -sqrt(3) / 2);
	p = ( - (b * b / (3 * a * a)) + (c / a));
	
	q = (((2 * b * b* b) / (27 * a * a * a)) - ((b * c) / (3 * a * a)) + (d / a));
	double Q;
	Q = ((p*p*p/27.) +(q*q/4.));
	T A, B;
	double x1, x2, x3, x2real, x3real, x2im, x3im, y1, y2, y3;
	A = pow(((-q / 2.) + sqrt(Q)), 1 / 3.);
	B = pow(((-q / 2.) - sqrt(Q)), 1 / 3.);
	cout << "Дискриминант кубического уравнения Q :" << Q<<endl;
	if (Q > 0) {
		y1 = A + B;
		x1 = y1 - (b / 3. * a);

		cout << "Уравнение имеет один действительный корень: X1=" << x1<<endl;
	}
	if (Q == 0) {
		y1 = ( - 2. * pow((q / 2.), 1 / 3.));
		y2= pow((q / 2.), 1 / 3.);
		y3 = pow((q / 2.), 1 / 3.);
		if (b == 0.) {
			cout << "Уравнение имеет три действительных корня:" << endl;
			cout<< "X1=" << y1 << endl;
			cout << "X2=" << y2 << endl;
			cout << "X3=" << y3 << endl;
		}
		else {
			x1 = (y1 - (b / 3. * a));
			x2 = (y2 - (b / 3. * a));
			x3 = (y3 - (b / 3. * a));
			cout << "Уравнение имеет три действительных корня:" << endl;
			cout<< "X1=" <<x1 << endl;
			cout << "X2=" << x2 << endl;
			cout << "X3=" << x3 << endl;
		}
	}
	if (Q < 0) {
		double i, j, k, l, m, n, f;
		i = pow((((q * q) / 4) - Q), 0.5);
		j = pow(i, 0.33333333333333333333333333333333); 
		k = acos((q / (2 * i)) * -1);
		l = j * -1;
		m = cos(k / 3);
		n = sqrt(3) * sin(k / 3);
		f = (b / (3 * a)) * -1;
		x1 = (2 * j) * m - (b / (3 * a));
		cout << "x1 = " << x1 << endl;
		x2 = l * (m + n) + f;
		cout << "x2 = " << x2 << endl;
		x3 = l * (m - n) + f;
		cout << "x3 = " << x3 << endl;
		
	};
	return 0;
}


int main() {

	SetConsoleOutputCP(1251);

	Polinomial<double> l = Polinomial<double>(3);
	l.set_coeff(3, 2);
	l.set_coeff(4, 1);
	l.set_coeff(5, 3);
	std::cout << l;

	Polinomial<double> t = Polinomial<double>(4);
	t.set_coeff(2, 2);
	t.set_coeff(6, 1);
	t.set_coeff(7, 3);
	t.set_coeff(7, 4);
	std::cout << t;

	std::cout << "Сложение многочленов:" << endl;
	Polinomial<double> c = t + l;
	std::cout << (l+=t);
	
	std::cout << "Вычитание многочленов:" << endl;
	Polinomial<double> w = Polinomial<double>(4);
	w.set_coeff(1, 2);
	w.set_coeff(2, 1);
	w.set_coeff(1, 3);
	w.set_coeff(7, 4);
	std::cout << t;
	Polinomial<double> d= l-w;
	std::cout << d;

	Polinomial<double> k = Polinomial<double>(3);
	k.set_coeff(3, 2);
	k.set_coeff(4, 1);
	k.set_coeff(5, 3);
	std::cout << "Вычисление значения при фиксированном значении х:" << endl;
	std::cout << k.calculate(4)<<endl;
	
	Polinomial<double> a = Polinomial<double>(3);
	a.set_coeff(3, 2);
	a.set_coeff(4, 1);
	a.set_coeff(0, 3);
	std::cout << "Проверка метода shrink_to_fit:";
	std::cout << a;
	std::cout << a.shrink_to_fit()<<endl;
	std::cout << a;

	Polinomial<double> b = Polinomial<double>(3);
	b.set_coeff(3, 2);
	b.set_coeff(8, 1);
	b.set_coeff(3, 3);
	
	std::cout << "Проверка метода expand:"<<endl;
	std::cout << b << endl;
	std::cout << b.expand(5)<<endl;
	b.set_coeff(3, 4);
	b.set_coeff(4, 5);
	std::cout << b;

	std::cout << "Проверка оператора[]:" << endl;
	std::cout << b[2]<<endl;

	Polinomial<double> f = Polinomial<double>(3);
	f.set_coeff(3, 2);
	f.set_coeff(4, 1);
	f.set_coeff(5, 3);

	std::cout << "Проверка операторов == и !=:" << endl;
	std::cout << (f==l)<<endl;
	std::cout << (f != l) << endl;

	Polinomial <std::complex<double>> v = Polinomial <std::complex<double>>(3);
	v.set_coeff(4, 2, 2);
	v.set_coeff(1, 1, 1);
	v.set_coeff(2, 5, 3);
	std::cout << v;

	Polinomial <std::complex<double>> s = Polinomial <std::complex<double>>(3);
	s.set_coeff( 3, 2, 2);
	s.set_coeff(3, 4, 1);
	s.set_coeff(3, 4, 3);
	std::cout << s;
	std::cout << "Сложение и вычитание комплексных многочленов:"<<endl;
	Polinomial <std::complex<double>> m = s + v;
	std::cout << m<<endl;
	m = s - v;
	std::cout << m << endl;

	Polinomial<double> r = Polinomial<double>(3);
	r.set_coeff(1, 3);
	r.set_coeff(0, 2);
	r.set_coeff(-12, 1);
	r.set_coeff(16, 0);
	std::cout << r << endl;

	std::cout << find_zero(1, 0, -12, 16)<<endl;

	Polinomial<double> two = Polinomial<double>(3);
	two.set_coeff(1, 3);
	two.set_coeff(3, 2);
	two.set_coeff(-3, 1);
	two.set_coeff(-14, 0);
	std::cout << two << endl;

	std::cout << find_zero(1, 3, -3, -14)<<endl;
	Polinomial<double> tree = Polinomial<double>(3);
	tree.set_coeff(1, 3);
	tree.set_coeff(0, 2);
	tree.set_coeff(-19, 1);
	tree.set_coeff(30, 0);
	std::cout << tree << endl;
	std::cout << find_zero(1, 0, -19, 30)<<endl;


	return 0;
}