#pragma once
#include <iostream>

using namespace std;
template <typename T>

class Polinomial {
	T* _coeff;
	int _degree;


public:
	Polinomial();
	Polinomial(const T* coeff, int degree);
	Polinomial(const Polinomial<T>& other);
	~Polinomial();

	void set_coeff(const T coeff);
	void set_degree(const int degree);
	T get_coeff() const;
	int get_degree() const;

	Polinomial operator+(const Polinomial& other) const;
	Polinomial operator-(const Polinomial& other) const;
	Polinomial operator*(double value) const;

	bool operator==(Polinomial& other) const;
	bool operator!=(Polinomial& other)const;

	void shrink_to_fit();

	friend std::ostream& operator <<(std::ostream& os, const Polinomial poly);

	};
