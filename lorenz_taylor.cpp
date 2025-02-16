#include <iostream>
#include <cmath>
#include <vector>


namespace peter
{
struct zero
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<std::vector<double> >&, const std::vector<double>&) const
	{	return 0.0;
	}
};
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
};
template<>
struct negate<zero>
{	typedef zero type;
};
enum enumParameter
{	eSigma,
	eRho,
	eBeta
};
template<enumParameter VAR>
struct parameter
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<std::vector<double> >&, const std::vector<double>&_rP) const
	{	return _rP.at(VAR);
	}
};
template<typename L, typename R>
struct multiplication;
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
};
template<typename L>
struct addition<L, zero>
{	typedef L type;
};
template<typename L>
struct addition<zero, L>
{	typedef L type;
};
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
};
template<typename L>
struct subtraction<L, zero>
{	typedef L type;
};
template<typename L>
struct subtraction<zero, L>
{	typedef typename negate<L>::type type;
};
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
};
template<typename L>
struct multiplication<L, zero>
{	typedef zero type;
};
template<typename L>
struct multiplication<zero, L>
{	typedef zero type;
};
template<typename T, std::size_t ORDER>
struct derive
{	typedef typename derive<T, ORDER-1>::type::derivative::type type;
};
template<typename T>
struct derive<T, 0>
{	typedef T type;
};
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
template<enumIndependent VAR, unsigned char ORDER>
struct X
{	struct derivative
	{	typedef X<VAR, ORDER+1> type;
	};
	double operator()(const std::vector<std::vector<double> >&_rI, const std::vector<double>&) const
	{	return _rI.at(VAR).at(ORDER);
	}
};
struct x:X<eX, 0>
{	struct derivative
	{	typedef typename multiplication<
			sigma,
			typename subtraction<X<eY, 0>, X<eX, 0>>::type
		>::type type;
	};
};
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
std::size_t factorial(const std::size_t _i)
{	if (_i)
		return _i*factorial(_i - 1);
	else
		return 1;
}
template<typename X, typename Y, typename Z, std::size_t ORDER>
void call(std::vector<std::vector<double> >&_rI, const std::vector<double> &_rP)
{	if constexpr (ORDER > 0)
		call<X, Y, Z, ORDER-1>(_rI, _rP);
	_rI.at(eX).at(ORDER) = typename derive<X, ORDER>::type()(_rI, _rP);
	_rI.at(eY).at(ORDER) = typename derive<Y, ORDER>::type()(_rI, _rP);
	_rI.at(eZ).at(ORDER) = typename derive<Z, ORDER>::type()(_rI, _rP);
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
	call<x, y, z, ORDER-1>(sI, sP);
	for (auto &r : sI)
		for (auto &d : r)
			std::cout << d/factorial(&d - r.data()) << std::endl;
}

