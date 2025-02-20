#include <iostream>
#include <cmath>
#include <vector>
#include "ctaylor/cjacobian.h"

namespace peter
{
template<typename L, typename R>
struct subtraction;
#if 1
/// the AD type used to calculate the derivatives of the calculated coefficients
/// wrt the initial conditions
/// If you do not care about this functionality, replace the #if 1 above with #if 0
typedef jacobian::cjacobian<
	boost::mp11::mp_list<
		boost::mp11::mp_size_t<0>,
		boost::mp11::mp_size_t<1>,
		boost::mp11::mp_size_t<2>
	>
> result;
template<std::size_t ENUM>
static result makeIndependent(const double _d)
{	using namespace boost::mp11;
	return jacobian::cjacobian<mp_list<mp_size_t<ENUM> > >(_d, true);
}
#else
/// calculate plain double coefficients
typedef double result;
template<std::size_t>
static result makeIndependent(const double _d)
{	return _d;
}
#endif

/// a zero value
struct zero
{	struct derivative
	{	typedef zero type;
	};
	auto operator()(const std::vector<std::vector<result> >&, const std::vector<double>&) const
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
	auto operator()(const std::vector<std::vector<result> >&_rI, const std::vector<double>&_rP) const
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
	auto operator()(const std::vector<std::vector<result> >&, const std::vector<double>&_rP) const
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
	auto operator()(const std::vector<std::vector<result> >& _rI, const std::vector<double>&_rP) const
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
/// a + (-b) == a - b
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
	auto operator()(const std::vector<std::vector<result> >& _rI, const std::vector<double>&_rP) const
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
	auto operator()(const std::vector<std::vector<result> >& _rI, const std::vector<double>&_rP) const
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
/// a*(-b) == -(a*b)
template<typename L, typename R>
struct multiplication<L, negate<R> >
{	typedef typename negate<
		typename multiplication<L, R>::type
	>::type type;
};
/// (-a)*b == -(a*b)
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
	auto operator()(const std::vector<std::vector<result> >&_rI, const std::vector<double>&) const
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
constexpr std::size_t factorial(const std::size_t _i)
{	if (_i)
		return _i*factorial(_i - 1);
	else
		return 1;
}
/// calculation of values of derivative
template<typename X, typename Y, typename Z, std::size_t ORDER, std::size_t MAX_ORDER>
void calculate(std::vector<std::vector<result> >&_rI, const std::vector<double> &_rP)
{	std::cout << "X<" << ORDER << ">=" << X() << std::endl;
	std::cout << "Y<" << ORDER << ">=" << Y() << std::endl;
	std::cout << "Z<" << ORDER << ">=" << Z() << std::endl;
	_rI.at(eX).at(ORDER) = X()(_rI, _rP);
	_rI.at(eY).at(ORDER) = Y()(_rI, _rP);
	_rI.at(eZ).at(ORDER) = Z()(_rI, _rP);
	if constexpr(ORDER < MAX_ORDER)
		calculate<typename X::derivative::type, typename Y::derivative::type, typename Z::derivative::type, ORDER + 1, MAX_ORDER>(_rI, _rP);
}
}
int main(int argc, char**argv)
{	if (argc != 7)
	{	std::cerr << argv[0] << ": Usage: " << argv[0] << " sigma rho beta x y z" << std::endl;
		std::cerr << argv[0] << ": Example: " << argv[0] << " 10 28 2.66666666 0.9 0 0" << std::endl;
		return 1;
	}
	const std::vector<double> sP(
		{	std::atof(argv[1]),
			std::atof(argv[2]),
			std::atof(argv[3])
		}
	);
	//using namespace boost::mp11;
	using namespace peter;
	static constexpr std::size_t MAX_ORDER = 3;
	std::vector<std::vector<result> > sI(
		{	std::vector<result>(MAX_ORDER + 1, makeIndependent<0>(std::atof(argv[4]))),
			std::vector<result>(MAX_ORDER + 1, makeIndependent<1>(std::atof(argv[5]))),
			std::vector<result>(MAX_ORDER + 1, makeIndependent<2>(std::atof(argv[6])))
		}
	);
	calculate<typename x::derivative::type, typename y::derivative::type, typename z::derivative::type, 1, MAX_ORDER>(sI, sP);
	for (auto &r : sI)
	{	for (auto &d : r)
			std::cout << d/factorial(&d - r.data()) << std::endl;
		std::cout << std::endl;
	}
}

