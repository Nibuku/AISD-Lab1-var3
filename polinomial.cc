#pragma once
#include <iostream>
#include <stdexcept>
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
	Polinomial(const int max_degree);
	Polinomial(const Polinomial<T>& other);
	~Polinomial();

	void set_coeff(const T coeff, int degree);
	void set_degree(const int degree);
	T get_coeff(int degree) const;
	int get_degree() const;

	Polinomial<T>& operator+(const Polinomial& other) const;
	Polinomial<T>& operator+=(const Polinomial<T>& other);
	Polinomial<T>& operator-(const Polinomial& other) const;
	Polinomial<T>& operator-=(const Polinomial& other);
	Polinomial<T>& operator*(double value) const;
	T& operator[](const int degree);

	bool operator==( const Polinomial& other) const;
	bool operator!=(const Polinomial& other)const;

	bool shrink_to_fit();
	bool expand(const int degree);
	T calculate(const T x) const;
	

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
	}
}

template<typename T>
Polinomial<T>::Polinomial(const Polinomial<T>& other) {
	_degree = other.get_degree;
	_coeff = new T[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = other._coeff[i];
	}
}

template<typename T>
Polinomial<T>::Polinomial(const int max_degree) : _degree(max_degree) {
	_coeff = new T[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = 1;
	}
}


template<typename T>
Polinomial<T>::~Polinomial() {
	for (int i = 0; i < _degree; ++i) {
		delete _coeff[i];
	}
	_coeff = nullptr;
	_degree = 0;
}

template <typename T>
Polinomial<T>& Polinomial<T>::operator+=(const Polinomial<T>& other) {
	if (_degree >= other._degree) {
		for (int i = 0; i <= other._degree; i++)
			_coeff[i] += other._coeff[i];
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
Polinomial<T>& Polinomial<T>:: operator+(const Polinomial& other) const {
	Polinomial<T>* result = new Polinomial<T>(*this);
	result += other;
	return result;
}

template <typename T>
Polinomial<T>& Polinomial<T>:: operator-=(const Polinomial& other) {
	if (_degree >= other._degree) {
		for (int i = 0; i <= other._degree; i++)
			_coeff[i] -= other._coeff[i];
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
Polinomial<T>& Polinomial<T>:: operator-(const Polinomial& other) const {
	Polinomial<T>* result = new Polinomial<T>(*this);
	result -= other;
	return result;
}

template <typename T>
Polinomial<T>& Polinomial<T>::operator*(double value) const {
	Polinomial<T>* result = new Polinomial<T>(*this);
	for (int i = 0; i <= _degree; ++i) {
		result->_coeff[i] = value * _coeff[i];
	}
	return *this;
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
		if ((this_copy._coeff[i] - other_copy._coeff[i]) >= ACCURACY)
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
std::ostream& operator<<(std::ostream& os, const Polinomial<T> poly) {
	
	for (int i = poly._degree; i >= 0; --i) {
		if (poly._coeff[i] == 0)
			continue;
		if (poly._degree == 0)
			os << poly[i]<<endl;
		os << poly[i] << "x^" << i << "+";
	}

	return os;
}

