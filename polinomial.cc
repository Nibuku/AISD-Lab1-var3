#include "polinomial.h"
#include <stdexcept>
#include <string>
#include <math.h> 

using namespace std;
//гетеры и сетеры
template<typename T>
void Polinomial<T>::set_coeff(const T coeff){
	this->_coeff = coeff;
}

template<typename T>
void Polinomial<T>::set_degree(const int degree) {
	this->_degree = degree;
}

template <typename T> 
T Polinomial<T>::get_coeff() const {
	return _coeffs;
}

template <typename T>
int Polinomial<T>::get_degree() const {
	return _degree;
}

//три конструктора
template<typename T>
Polinomial<T>::Polinomial() : _degree(0) {
	_coeff = new T[1];
	_coeff[0] = 0;
}
template<typename T>

Polinomial<T>::Polinomial(const T* coeff, int degree): _degree(degree) {
	_coeff = new T[_degree + 1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = coeff[i];
	}
}

template<typename T>
Polinomial<T>::Polinomial(const Polinomial<T>& other) {
	_degree = other.get_degree;
	_coeff = new T[_degree+1];
	for (int i = 0; i <= _degree; ++i) {
		_coeff[i] = other._coeff[i];
	}
}

