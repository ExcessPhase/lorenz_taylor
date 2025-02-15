#include <iostream>
#include <cmath>
#include <vector>


namespace peter
{
#if 0
template<unsigned char VAR, unsigned int ORDER>
struct independentVariable
{	struct derivative
	{	typedef independentVariable<VAR, ORDER+1> type;
	};
	double operator()(const std::vector<std::vector<double> >& _rI, const std::vector<double>&) const
	{	return _rI.at(VAR).at(ORDER);
	}
};
#endif
struct zero
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<double>&, const std::vector<double>&) const
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
	double operator()(const std::vector<double>&_rI, const std::vector<double>&_rP) const
	{	return -T()(_rI, _rP);
	}
	typedef negate type;
};
template<>
struct negate<zero>
{	typedef zero type;
};
template<unsigned char VAR>
struct parameter
{	struct derivative
	{	typedef zero type;
	};
	double operator()(const std::vector<double>&, const std::vector<double>&_rP) const
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
	double operator()(const std::vector<double>& _rI, const std::vector<double>&_rP) const
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
	double operator()(const std::vector<double>& _rI, const std::vector<double>&_rP) const
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
	double operator()(const std::vector<double>& _rI, const std::vector<double>&_rP) const
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
enum enumParameter
{	eSigma,
	eRho,
	eBeta
};
typedef parameter<eSigma> sigma;
typedef parameter<eRho> rho;
typedef parameter<eBeta> beta;
#if 0
dx/dt = sigma * (y - x)

dy/dt = x * (rho - z) - y

dz/dt = x * y - beta * z
#endif
struct x;
struct y;
struct z;


enum enumIndependent
{	eX,
	eY,
	eZ
};
struct z
{	struct derivative
	{	typedef subtraction<
			multiplication<
				x,
				y
			>,
			multiplication<
				beta,
				z
			>
		> type;
	};
	double operator()(const std::vector<double>& _rI, const std::vector<double>&) const
	{	return _rI.at(eZ);
	}
};
struct y
{	struct derivative
	{	typedef subtraction<
			multiplication<
				x,
				subtraction<rho, z>
			>,
			y
		> type;
	};
	double operator()(const std::vector<double>& _rI, const std::vector<double>&) const
	{	return _rI.at(eY);
	}
};
struct x
{	struct derivative
	{	typedef multiplication<
			sigma,
			subtraction<y, x>
		> type;
	};
	double operator()(const std::vector<double>& _rI, const std::vector<double>&) const
	{	return _rI.at(eX);
	}
};
std::size_t factorial(const std::size_t _i)
{	if (_i)
		return _i*factorial(_i - 1);
	else
		return 1;
}
template<typename X, std::size_t ORDER>
void call(void)
{	if constexpr (ORDER > 0)
		call<X, ORDER-1>();
	const std::vector<double> sP({10.0, 28.0, 8.0/3.0});
	const std::vector<double> sI({0.9, 0.0, 0.0});
	std::cout << typename derive<X, ORDER>::type()(sI, sP)/factorial(ORDER) << std::endl;
}
}
int main(int argc, char**argv)
{	using namespace peter;
	call<peter::x, 3>();
	std::cout << std::endl;
	call<peter::y, 3>();
	std::cout << std::endl;
	call<peter::z, 3>();
}

