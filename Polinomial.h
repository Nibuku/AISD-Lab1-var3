#pragma once
#include <iostream>

using namespace std;
template <typename T>

class Polinomial {
	T* _coeff;
	int _degree;

	void shrink_to_fit();

public:
	Polinomial();
	Polinomial(const T* coeff, int degree);
	Polynomial(const Polynomial<T>& other);
	~Polynomial();

	void set_coeff(const T coeff);
	void set_degree(const int degree);
	T get_coeff() const;
	T get_degree() const;

	Polynomial operator+(const Polinomial& other) const;
	Polynomial operator-(const Polinomial& other) const;
	Polynomial operator*(double value) const;

	bool operator==(Polynomial& other) const;
	bool operator!=(Polynomial& other)const;

	friend std::ostream& operator <<(std::ostream& os, const Polynomial poly);

};

