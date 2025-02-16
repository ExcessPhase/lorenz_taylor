#include <iostream>
#include <cmath>
#include <vector>


namespace peter
{
template<typename L, typename R>
struct subtraction;
/// a zero value
struct zero
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<std::vector<double> >&, const std::vector<double>&) const
	{	return 0.0;
	}
	friend std::ostream &operator<<(std::ostream &_rS, const zero&)
	{	return _rS << "0";
	}
};
/// a unary negative sign
template<typename T>
struct negate
{	struct derivative
	{	typedef typename negate<
			typename T::derivative::type
		>::type type;
	};
	double operator()(const std::vector<std::vector<double> >&_rI, const std::vector<double>&_rP) const
	{	return -T()(_rI, _rP);
	}
	typedef negate type;
	friend std::ostream &operator<<(std::ostream &_rS, const negate<T>&)
	{	return _rS << "(-" << T() << ")";
	}
};
/// -0 == 0
template<>
struct negate<zero>
{	typedef zero type;
};
enum enumParameter
{	eSigma,
	eRho,
	eBeta
};
/// a parameter value
template<enumParameter VAR>
struct parameter
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<std::vector<double> >&, const std::vector<double>&_rP) const
	{	return _rP.at(VAR);
	}
	friend std::ostream &operator<<(std::ostream &_rS, const parameter<VAR>&)
	{	switch (VAR)
		{	default:
				return _rS << "parameter<" << VAR << ">";
			case eSigma:
				return _rS << "sigma";
			case eRho:
				return _rS << "rho";
			case eBeta:
				return _rS << "beta";
		}
	}
};
template<typename L, typename R>
struct multiplication;
/// binary +
template<typename L, typename R>
struct addition
{	struct derivative
	{	typedef typename addition<
			typename L::derivative::type,
			typename R::derivative::type
		>::type type;
	};
	double operator()(const std::vector<std::vector<double> >& _rI, const std::vector<double>&_rP) const
	{	return L()(_rI, _rP) + R()(_rI, _rP);
	}
	typedef addition type;
	friend std::ostream &operator<<(std::ostream &_rS, const addition<L, R>&)
	{	return _rS << "(" << L() << "+" << R() << ")";
	}
};
/// L + 0 == L
template<typename L>
struct addition<L, zero>
{	typedef L type;
};
/// 0 + L == L
template<typename L>
struct addition<zero, L>
{	typedef L type;
};
template<typename L, typename R>
struct addition<L, negate<R> >
{	typedef typename subtraction<L, R>::type type;
};
/// binary minus
template<typename L, typename R>
struct subtraction
{	struct derivative
	{	typedef typename subtraction<
			typename L::derivative::type,
			typename R::derivative::type
		>::type type;
	};
	double operator()(const std::vector<std::vector<double> >& _rI, const std::vector<double>&_rP) const
	{	return L()(_rI, _rP) - R()(_rI, _rP);
	}
	typedef subtraction type;
	friend std::ostream &operator<<(std::ostream &_rS, const subtraction<L, R>&)
	{	return _rS << "(" << L() << "-" << R() << ")";
	}
};
/// L - 0 == L
template<typename L>
struct subtraction<L, zero>
{	typedef L type;
};
/// 0 - L == -L
template<typename L>
struct subtraction<zero, L>
{	typedef typename negate<L>::type type;
};
/// multiplication
template<typename L, typename R>
struct multiplication
{	struct derivative
	{	typedef typename addition<
			typename multiplication<
				typename L::derivative::type,
				R
			>::type,
			typename multiplication<
				L,
				typename R::derivative::type
			>::type
		>::type type;
	};
	double operator()(const std::vector<std::vector<double> >& _rI, const std::vector<double>&_rP) const
	{	return L()(_rI, _rP) * R()(_rI, _rP);
	}
	typedef multiplication type;
	friend std::ostream &operator<<(std::ostream &_rS, const multiplication<L, R>&)
	{	return _rS << "(" << L() << "*" << R() << ")";
	}
};
/// L*0 == 0
template<typename L>
struct multiplication<L, zero>
{	typedef zero type;
};
/// 0*L == 0
template<typename L>
struct multiplication<zero, L>
{	typedef zero type;
};
template<typename L, typename R>
struct multiplication<L, negate<R> >
{	typedef typename negate<
		typename multiplication<L, R>::type
	>::type type;
};
/// 0*L == 0
template<typename L, typename R>
struct multiplication<negate<L>, R>
{	typedef typename negate<
		typename multiplication<L, R>::type
	>::type type;
};
/// the parameters
typedef parameter<eSigma> sigma;
typedef parameter<eRho> rho;
typedef parameter<eBeta> beta;
#if 0
dx/dt = sigma * (y - x)

dy/dt = x * (rho - z) - y

dz/dt = x * y - beta * z
#endif
enum enumIndependent
{	eX,
	eY,
	eZ
};
/// an independent variable VAR with derivative order
template<enumIndependent VAR, unsigned char ORDER>
struct X
{	struct derivative
	{	typedef X<VAR, ORDER+1> type;
	};
	double operator()(const std::vector<std::vector<double> >&_rI, const std::vector<double>&) const
	{	return _rI.at(VAR).at(ORDER);
	}
	friend std::ostream &operator<<(std::ostream &_rS, const X<VAR, ORDER>&)
	{	switch (VAR)
		{	default:
				return _rS << "X<" << VAR << ", " << int(ORDER) << ">";
			case eX:
				return _rS << "X<" << int(ORDER) << ">";
			case eY:
				return _rS << "Y<" << int(ORDER) << ">";
			case eZ:
				return _rS << "Z<" << int(ORDER) << ">";
		}
	}
};
/// independent variable x
/// operator() is inherited
struct x:X<eX, 0>
{	struct derivative
	{	typedef typename multiplication<
			sigma,
			typename subtraction<X<eY, 0>, X<eX, 0>>::type
		>::type type;
	};
};
/// independent variable y
/// operator() is inherited
struct y:X<eY, 0>
{	struct derivative
	{	typedef typename subtraction<
			typename multiplication<
				X<eX, 0>,
				typename subtraction<rho, X<eZ, 0>>::type
			>::type,
			X<eY, 0>
		>::type type;
	};
};
/// independent variable z
/// operator() is inherited
struct z:X<eZ, 0>
{	struct derivative
	{	typedef typename subtraction<
			typename multiplication<
				X<eX, 0>,
				X<eY, 0>
			>::type,
			typename multiplication<
				beta,
				X<eZ, 0>
			>::type
		>::type type;
	};
};
/// recursive implementation of factorial
std::size_t factorial(const std::size_t _i)
{	if (_i)
		return _i*factorial(_i - 1);
	else
		return 1;
}
/// calculation of values of derivative
template<typename X, typename Y, typename Z, std::size_t ORDERM1, std::size_t ORDER>
void calculate(std::vector<std::vector<double> >&_rI, const std::vector<double> &_rP)
{	std::cout << "X<" << ORDER - ORDERM1 << ">=" << X() << std::endl;
	std::cout << "Y<" << ORDER - ORDERM1 << ">=" << Y() << std::endl;
	std::cout << "Z<" << ORDER - ORDERM1 << ">=" << Z() << std::endl;
	_rI.at(eX).at(ORDER - ORDERM1) = X()(_rI, _rP);
	_rI.at(eY).at(ORDER - ORDERM1) = Y()(_rI, _rP);
	_rI.at(eZ).at(ORDER - ORDERM1) = Z()(_rI, _rP);
	if constexpr (ORDERM1 > 1)
		calculate<typename X::derivative::type, typename Y::derivative::type, typename Z::derivative::type, ORDERM1-1, ORDER>(_rI, _rP);
}
}
int main(int argc, char**argv)
{	if (argc != 7)
	{	std::cerr << argv[0] << ": Usage: " << argv[0] << "sigma rho beta x y z" << std::endl;
		return 1;
	}
	static constexpr std::size_t ORDER = 4;
	const std::vector<double> sP(
		{	std::atof(argv[1]),
			std::atof(argv[2]),
			std::atof(argv[3])
		}
	);
	std::vector<std::vector<double> > sI(
		{	std::vector<double>(ORDER, std::atof(argv[4])),
			std::vector<double>(ORDER, std::atof(argv[5])),
			std::vector<double>(ORDER, std::atof(argv[6]))
		}
	);
	using namespace peter;
	calculate<typename x::derivative::type, typename y::derivative::type, typename z::derivative::type, ORDER-1, ORDER>(sI, sP);
	for (auto &r : sI)
	{	for (auto &d : r)
			std::cout << d/factorial(&d - r.data()) << std::endl;
		std::cout << std::endl;
	}
}

