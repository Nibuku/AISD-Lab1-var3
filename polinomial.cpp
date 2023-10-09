#pragma once
#include <iostream>
#include<Windows.h>
#include <stdexcept>
#include <complex>
#include <cmath>

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
	complex<T> set_coeff(int degree);
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
	if (degree == 1) {
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
	if (degree == 1) {
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
		_coeff[i]->imag(coeff[i].real());
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
		_coeff[i]->imag(_coeff[i].real());
	};
}

template<typename T>
Polinomial<T>::Polinomial(const int max_degree) : _degree(max_degree) {
	_coeff = new T[_degree + 1];
	_coeff[0] = 1;
	for (int i = 1; i <= _degree; ++i) {
		_coeff[i]=1;
	}
}

template<typename T>
Polinomial<T>::~Polinomial() {
	delete[]  _coeff;
	_coeff = nullptr;
	_degree = 0;
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
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			_coeff[i] += other._coeff[i];
		}
	}
	return *this;

}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>::operator+=(const Polinomial<complex<T>>& other) {
	if (_degree <= other._degree) {
		for (int i = 0; i <= other._degree; i++) {
			_coeff[i].real() += other._coeff[i].real();
			_coeff[i].imag() += other._coeff[i].imag();
		}
	}
	else {
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{ _coeff[i].real() += other._coeff[i].real();
		  _coeff[i].imag() += other._coeff[i].imag();
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
	if (_degree <= other._degree) {
		for (int i = 0; i <= other._degree; i++)
		{
			_coeff[i] -= other._coeff[i];
		}
	}
	else {
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			_coeff[i] -= other._coeff[i];
		}
	}
	return *this;

}

template <typename T>
Polinomial<complex<T>>& Polinomial<T>:: operator-=(const Polinomial<complex<T>>& other) {
	if (_degree <= other._degree) {
		for (int i = 0; i <= other._degree; i++) {
			_coeff[i].real() += other._coeff[i].real();
			_coeff[i].imag() += other._coeff[i].imag();
		}
	}
	else {
		_degree = other._degree;
		for (int i = 0; i <= _degree; i++)
		{
			_coeff[i].real() += other._coeff[i].real();
			_coeff[i].imag() += other._coeff[i].imag();
		}
	}
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
		result.imag() += _coeff[i].imag() * pow(x.imag(), i);
	}
	return result;
}

template<typename T>
bool Polinomial<T>::shrink_to_fit() {
	int count = 0;
	for (int i = _degree; i > 0; --i) {
		if (_coeff[i] == 0)
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
bool Polinomial<T>::expand(const int degree) {
	if (degree <= _degree)
		return false;
	T* new_coeff = new T[degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		new_coeff[i] = _coeff[i];
		
	}
	_degree = degree;
	delete[] _coeff;
	_coeff = new_coeff;
	return true;
}

template<typename T>
bool Polinomial<T>::operator==(const Polinomial<T>& other) const {
	Polinomial<T> this_copy(*this);
	Polinomial<T> other_copy(other);

	this_copy.shrink_to_fit();
	other_copy.shrink_to_fit();

	if (this_copy._degree != other_copy._degree)
		return false;

	for (int i = 0; i <= this_copy._degree; ++i) {
		if ((this_copy._coeff[i] - other_copy._coeff[i]) > ACCURACY)
			return false;
	}
	return true;
}

template<typename T>
bool Polinomial<T>::operator==(const Polinomial<complex<T>>& other) const {
	Polinomial<complex<T>> this_copy(*this);
	Polinomial<complex<T>> other_copy(other);


	if (this_copy._degree != other_copy._degree)
		return false;

	for (int i = 0; i <= this_copy._degree; ++i) {
		if (this_copy._coeff[i].real() - other_copy._coeff[i].real() > ACCURACY)
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
	os << "(" << poly.get_coeff(0).real() << "+" << poly.get_coeff(0).imag() << "i" << ")x^" << endl;



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
T find_zero(const Polinomial<T> poly) {
	if (poly.get_degree() == 3)
	{
		throw std::runtime_error("Uncorrect degree");
	}
	T c = poly.get_coeff(0);
	T d = poly.get_coeff(1);
	T b = poly.get_coeff(2);
	T a = poly.get_coeff(3);

	T discr = b * b - 3 * a * c;
	T sqrt_discr = std::sqrt(std::abs(discr));
	T cuberoot;
	if (discr < 0)
		return false;
	if (discr > 0) {
		cuberoot = (-b + sqrt_discr) / (3 * a);
	}
	else {
		cuberoot = (-b) / (3 * a);
	}

	std::complex<T> first_root = cuberoot;
	std::complex<T> omega = std::complex<double>(-0.5, std::sqrt(3.0) / 2.0);
	std::complex<T> second_root = omega * cuberoot;
	std::complex<T> third_root = omega * omega * cuberoot;

	std::complex<T> first_value = a * first_root * first_root * first_root+ b * first_root * first_root+ c * first_root+ d;
	std::complex<T> second_value = a * second_root * second_root * second_root+ b * second_root * second_root+ c * second_root+ d;
	std::complex<T> third_value = a * third_root * third_root * third_root+ b * third_root * third_root+ c * third_root+ d;
	cout << first_value << endl;
	cout << second_value << endl;
	cout << third_value << endl;
}




int main() {

	SetConsoleOutputCP(1251);

	Polinomial<double> l = Polinomial<double>(3);
	l.set_coeff(3, 2);
	l.set_coeff(4, 1);
	l.set_coeff(5, 3);
	std::cout << l;

	Polinomial<double> t = Polinomial<double>(3);
	t.set_coeff(2, 2);
	t.set_coeff(6, 1);
	t.set_coeff(7, 3);
	std::cout << t;

	std::cout << "Сложение многочленов:" << endl;
	Polinomial<double> c = t + l;
	std::cout << c;
	
	std::cout << "Вычитание многочленов:" << endl;
	Polinomial<double> d= t -l;
	std::cout << d;

	std::cout << "Вычисление значения при фиксированном значении х:" << endl;
	std::cout << l.calculate(4)<<endl;
	
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
	std::cout << v;

	Polinomial <std::complex<double>> s = Polinomial <std::complex<double>>(3);
	s.set_coeff( 3, 2, 2);
	s.set_coeff(3, 4, 1);
	std::cout << s;
	std::cout << "Сложение комплексных и вычитание многочленов:"<<endl;
	Polinomial <std::complex<double>> m = s + v;
	std::cout << m<<endl;
	m = s - v;
	std::cout << m << endl;

	cout << find_zero(t) <<endl;

	return 0;
}